#####################################################################
## Script to compute differential methylation at the feature level ##
#####################################################################
args <- commandArgs()

help <- function() {
    cat("difmet.R :
Calculates differential methylation and outputs a table with the information\n")
    cat("Usage: \n")
    cat("--meta       : Sample metadata with QC information            [required]\n")
    cat("--met        : path to methylation files                      [required]\n")
    cat("--anno       : path to annotation files                       [required]\n")
    cat("--genemeta   : Gene metadata file                             [required]\n")
    cat("--groups     : File with grouping from transcriptome analysis [required]\n")
    cat("--mincells   : Min cells per feature in each group            [required]\n")
    cat("--out        : File for output                                [required]\n")
    cat("\n")
    q()
}

io <- list()
opts <- list()

if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$sample.metadata   <- sub( '--meta=', '', args[grep('--meta', args)] )
    io$data.dir          <- sub( '--met=', '', args[grep('--met', args)] )
    io$annos_dir         <- sub( '--anno=', '', args[grep('--anno',args)] )
    io$gene.metadata     <- sub( '--genemeta=', '', args[grep('--genemeta', args)] )
    io$groups            <- sub( '--groups=', '', args[grep('--groups', args)] )
    io$outfile           <- sub( '--out=', '', args[grep('--out', args)] )
    opts$min.cells       <- sub( '--mincells=', '', args[grep('--mincells', args)] )
}



library(data.table)
library(purrr)
library(ggplot2)
library(parallel)

### functions ###
fread_gz = function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") return(fread(cmd = paste("zcat", filename), ...))
  fread(filename, ...)
}
fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")

## TEST INPUT ##
## Define I/O ##
#io <- list()
#io$basedir <- "/Volumes/Data/data/hisham/"
#io$sample.metadata <- "../scNMT_NOMeWorkFlow/tables//sample_stats_qcPass.txt" 
#io$data.dir <- "../scNMT_NOMeWorkFlow/data/met"
#io$annos_dir  <- "../scNMT_NOMeWorkFlow/data/anno"
#io$gene.metadata <- "../scRNA_SMARTseq2/data/gene_metadata.tsv"
#io$groups <- "data/groups.tsv"
#io$outfile <- "tables/difmet.tsv"

## Define options ##

# Define stage and lineage
opts$groupA <- "M7C"
opts$groupB <- "M7E"
opts$groupC <- "TDC"
opts$groupD <- "TDE"

# Overlap genomic features with nearby genes?
opts$OverlapWithGenes <- TRUE
opts$gene_window <- 5e4       # window length for the overlap

# Subset top most variable sites
opts$number_features <- 5000

# Filter by coverage
opts$min.CpGs <- 1            # Minimum number of CpG per feature in each cell
opts$min.cells <- as.numeric(opts$min.cells)  # Minimum number of cells per feature in each group

# Regress out global methylation rate
opts$regress.mean <- FALSE

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"
if (opts$regress.mean) {
  warning("Binomial test does not work when regressing out the global methylation, switching to t.test")
  opts$statistical.test <- "t.test"
}

# Minimum differential methylation (%) for statistical significance
opts$min.diff <- 5

# Multiple testing correction
opts$threshold_fdr <- 0.10


###############
## Load data ##
###############

# Load sample metadata and grouping
groups <- fread(io$groups)
sample_metadata <- fread(io$sample.metadata)
#sample_metadata$id <- gsub("(sc_[A-H][0-9]+)_.*","\\1", sample_metadata$id)

sample_metadata$id_rna <- gsub("BSM7E6", "M7E6A", sample_metadata$sample)
sample_metadata$id_rna <- sub("T_", "TD", sample_metadata$id_rna)
sample_metadata$id_rna <- sub('_S.*', '', sample_metadata$id_rna)

#sample_metadata$id_rna <- sub("_","",sample_metadata$sample)

groups$id_rna <- sub('_S.*', '', groups$id_rna)

#groups$sample <- sub("D","_",groups$id_rna)

sample_metadata <- sample_metadata %>% 
  .[context == "CG" & pass_metQC == TRUE & pass_CHHQC == TRUE & pass_CHGQC == TRUE] %>%
  merge(groups, by = "id_rna")

  
opts$cells <- sample_metadata[, id_rna]

# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>%
  .[, chr := gsub("chr", "", chr) %>% paste0("chr", .)] %>%
  .[, c("start", "end") := .(start - opts$gene_window, end + opts$gene_window)] %>%
  setkey(chr, start, end)

# Load genomic context metadata 
feature_metadata <- dir(io$annos_dir, pattern = ".bed$", full = TRUE) %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
  .[, chr := gsub("chr", "", chr) %>% paste0("chr", .)]

# run overlap
feature_metadata <- feature_metadata[, isgene := grepl("ENS", id) %>% as.numeric] %>% 
  split(by = "isgene", keep.by = FALSE)

feature_metadata$`0` <- setkey(feature_metadata$`0`, chr, start, end) %>%
  foverlaps(gene_metadata)

feature_metadata$`1` <- feature_metadata$`1`[, ens_id := id] %>% 
  merge(gene_metadata, by = "ens_id")

feature_metadata$`1` <- feature_metadata$`1`[,2:8]
colnames(feature_metadata$`1`) <- c("chr", "start", "end", "strand", "id", "anno", "gene")

feature_metadata$`0` <- feature_metadata$`0`[,c("chr", "i.start", "i.end", "i.strand", "id", "anno", "gene")]
colnames(feature_metadata$`0`) <- c("chr", "start", "end", "strand", "id", "anno", "gene")

feature_metadata <- map(feature_metadata, ~.[, .(id, anno, gene, chr, start, end, strand)]) %>%
  rbindlist()


# Load methylation data
data <- dir(io$data.dir, pattern = ".tsv.gz", full = TRUE) %>% 
  map(fread_gz) %>% 
  rbindlist()

############################
## Parse methylation data ##
############################



# Define the two exclusive groups
sample_metadata[group == opts$groupA, groupABC := "A"]
sample_metadata[group == opts$groupB, groupABC := "B"]
sample_metadata[group == opts$groupC, groupABC := "C"]
sample_metadata[group == opts$groupD, groupABC := "D"]

# Merge methylation data and sample metadata
data <- merge(data, sample_metadata[, .(sample, group = groupABC)], by = "sample")

data[is.na(data)] <- 0

# Convert beta value to M value
if (opts$statistical.test == "t.test")
  data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#############################
## Filter methylation data ##
#############################

# Filter features by coverage
data <- data[N>=opts$min.CpGs]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by=c("id","anno")] %>% 
  .[Ngroup==2] %>% 
  .[,Ngroup:=NULL]

# Filter features by minimum number of cells per group
remove_n_sites <- data %>% split(.$anno) %>% map(~ .[,.(N=min(.N)), by=c("id","group")] %>% .[N<opts$min.cells,id])

data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[!id %in% remove_n_sites[[y]]]) %>% rbindlist

# Filter by variance
if (!is.na(opts$number_features)) {
  keep_hv_sites <- data %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% setorder(-var)  %>% head(n=opts$number_features) %>% .$id)
  data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist
}



###############################################
## Associate the genomic features with genes ##
###############################################

if (opts$OverlapWithGenes==TRUE) {
  data <- merge(data, feature_metadata, by = c("id", "anno"), allow.cartesian = TRUE)} else {
  data[,gene:="NA"]
}


#######################################
## Differential methylation analysis ##
#######################################

difmetacc <- function(comparison) {
# Binomial assumption: test of equal proportions using Fisher exact test
    variable_list <- LETTERS[1:comparison]
    for (i in variable_list) {
        for (j in variable_list) {
            if (j > i) {
                A <- i
                B <- j
                if (opts$statistical.test == "binomial") {
                    diff <- data[, .(
                        A_met=sum(.SD[group==A,N*(rate/100)]), A_unmet=sum(.SD[group==A,N*(1-rate/100)]),
                        B_met=sum(.SD[group==B,N*(rate/100)]), B_unmet=sum(.SD[group==B,N*(1-rate/100)])), by = c("id","anno", "gene","chr","start","end","strand")] %>%
                        .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id","anno", "gene","chr","start","end","strand")] %>%
                        .[,c(paste0("rate",A),paste0("rate",B)):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))]
    # T-test under normality assumption
                } else if (opts$statistical.test == "t.test") {
                    diff <- data[, .(
                        N_A = .SD[group==A,.N], N_B = .SD[group==B,.N],
                        rateA = mean(.SD[group==A,rate]), rateB = mean(.SD[group==B,rate]),
                        p.value = t.test(x=.SD[group==B,m], y=.SD[group==A,m], var.equal=FALSE)[["p.value"]]), by = c("id","anno", "gene","chr","start","end","strand")]
                }

# Multiple testing correction and define significant hits
                diff %>%
                    .[,diff:=diff[,13]-diff[,14]] %>%
                    .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr")), by="anno"] %>%
                    .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
                    .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>%
                    .[, comparison := paste0(A, "_vs_", B)] %>%
                    setorderv("padj_fdr")
                diff[,9] <- round(diff[,9],2)
                diff[,10] <- round(diff[,10],2)
                diff[,11] <- round(diff[,11],2)
                fwrite_tsv(diff, paste0(io$outfile, A, "_vs_", B, ".tsv"))
            }
        }
    }
}

nclusters <- length(unique(data$group))

mclapply(nclusters, difmetacc, mc.cores=nclusters*2)
