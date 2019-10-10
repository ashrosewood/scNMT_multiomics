args <- commandArgs()

help <- function(){
    cat("correlations_met.R :
Correlates methylation and RNA\n")
    cat("Usage: \n")
    cat("--anno        : path to annotation files                               [required]\n")
    cat("--meta        : stats table with QC info                               [required]\n")
    cat("--SO          : Seurat object with RNA information                     [required]\n")
    cat("--accdir      : path to annotated methylation files                    [required]\n")
    cat("--genefile    : gene metadata file                                     [required]\n")
    cat("--plotdir      : output directory                                      [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$anno_dir   <- sub( '--anno=', '', args[grep('--anno', args)] )
    io$rna_sce    <- sub( '--SO=', '', args[grep('--SO', args)] )
    io$meta_data  <- sub( '--meta=', '', args[grep('--meta', args)] )
    io$plot_dir   <- sub( '--plotdir=', '', args[grep('--plotdir', args)] )
    io$met_dir    <- sub( '--accdir=', '', args[grep('--accdir', args)] )
    io$gene_file  <- sub( '--genefile=', '', args[grep('--genefile', args)] )
}


library(Seurat)
library(scater)
library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)

### TEST INPUT ###
#io$meta_data <- "../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt"
#io$met_dir <- "data/met"
#io$rna_sce <- "../scNMT_transcriptomeMapping/data/SeuratObject.rds"
#io$anno_dir <- "../test_git/scNMT_NOMeWorkFlow/data/anno"
#io$gene_file <- "../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv"
#io$plot_dir <- "plots/cor_met"


opts$anno_regex <- ""
opts$gene_overlap_dist <- 1e5 # overlap annoations with genes within xx bp
opts$min_weight_met <- 1
opts$min_cells_met <- 10 # loci must have observations in this many cells
opts$min_cells_rna <- 5 # genes must have exp>0 in this many cells
opts$anno_regex <- "CGI_promoter|_ER_peaks|H3K27ac_peaks|body|Repressed|Enhancer|CTCF"
opts$filt_rna_var <- 0.5 # select the top xx fraction by variance
opts$filt_met_var <- 0.5 # select the top xx fraction by variance

opts$p_cutoff <- 0.1


### functions ###
fread_gz <- function(path, ...){fread(cmd = paste("zcat", path), ...)}

fread_gz <- function(path, ...){fread(x, ...)}

### load metadata and select cells ####

meta <- fread(io$meta_data) %>%
  .[pass_accQC == TRUE & pass_metQC == TRUE]

### load rna and format as data.table ###
rna <- readRDS(io$rna_sce)
rna <- rna@assays$RNA@data
colnames(rna) <- str_extract(colnames(rna),"_[A-Z0-9]+_")
colnames(rna) <- gsub("([A-Z])0([0-9])","\\1\\2",colnames(rna))

meta$id_rna <- str_extract(meta$id,"_[A-Z0-9]+_")
rna <- as.data.frame(rna)
rna$ens_id <- rownames(rna)
#rna <- rna[,unique(meta$id_rna)]
#rna <- setDT(rna, keep.rownames = "ens_id")
#rna <- melt(rna, id.vars = "ens_id", value.name = "exp", variable.name = "id_rna")
#rna <- merge(rna, meta[, .(id_rna,sample)], by = "id_rna", allow.cartesian=TRUE)

rna <- rna %>%
  .[, unique(meta[, id_rna])] %>%
  setDT(keep.rownames = "ens_id") %>%
  melt(id.vars = "ens_id", value.name = "exp", variable.name = "id_rna") %>%
  merge(meta[, .(id_rna, sample)], by = "id_rna",allow.cartesian=TRUE)

### filter rna ###

rna <- rna[, keep := sum(exp>0) > opts$min_cells_rna, ens_id] %>%
  .[keep == TRUE]

filt_genes_var <- rna[, var(exp), ens_id] %>%
  .[order(-rank(V1))] %>%
  .[1: (.N * opts$filt_rna_var), ens_id]

rna <- rna[ens_id %in% filt_genes_var]

rna$gene <- rna$ens_id

### load met data ###

met <- dir(io$met_dir, pattern = ".gz$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  # map(fread_gz) %>%
  map(fread) %>%
  rbindlist()


### annotate with gene name and keep only loci with nearby genes ###

genes <- fread(io$gene_file) %>%
  .[, c("start", "end", "chr") := .(start - opts$gene_overlap_dist, 
                                    end + opts$gene_overlap_dist,
                                    gsub("chr", "", chr))] %>%
  setkey(chr, start, end)

 
anno <- dir(io$anno_dir, full = TRUE, pattern = ".bed$") %>%
  .[grep(opts$anno_regex, .)] %>%
  map(fread) %>%
  rbindlist() %>% setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
  .[, .(chr = gsub("chr", "", chr),
        start = start,
        end = end,
        strand = strand,
        id = id,
        anno = anno
  )] %>% .[!chr %in% c("MT", "Y")] %>%
  .[, gene_id := grepl("ENSG", id)] %>%
  split(by = "gene_id", keep.by = FALSE)


anno[["TRUE"]] <- anno[["TRUE"]] %>%
  .[, ens_id := id] %>%
  merge(genes[, .(ens_id, gene)], by = "ens_id")

anno[["FALSE"]] <- setkey(anno[["FALSE"]], chr, start, end) %>%
  foverlaps(genes, nomatch = 0L)

anno <- map(anno, ~.[, .(id, anno, gene, ens_id)]) %>%
  rbindlist()


met <- merge(met, anno, by = c("anno", "id"), allow.cartesian = TRUE) # note some loci have >1 gene -> cartesian join

### join met with rna ###

metrna <- merge(met, rna, by = c("sample", "gene"))



### filter met data ###

metrna <- metrna[N >= opts$min_weight_met] %>%
  .[, keep := .N > opts$min_cells_met, .(id, anno, gene, ens_id.x)] %>%
  .[keep == TRUE]

filt_met_var <- metrna[, var(rate), .(id, anno, gene, ens_id.x)] %>%
  unique(by = c("id", "anno", "V1")) %>% # need to filter out multiple instances due to some loci having >1 gene
  split(by = "anno") %>%
  map(~.[order(-rank(V1))][1: (.N * opts$filt_met_var), .(id, anno)]) %>%
  rbindlist() %>%
  setkey(id, anno)

metrna <- metrna[filt_met_var, on = c("id", "anno")]

### compute correlation coeff ###

compute_cor <- function(exp, rate, N){
  wtd.cor(exp, rate, N) %>%
    as.list() %>%
    set_names(c("r", "std_err", "t", "p"))
}

cors <- metrna[, compute_cor(exp, rate, N), .(anno, id, gene, ens_id.x)] %>%
  .[, padj := p.adjust(p, method = "fdr"), .(anno)] %>%
  .[, logpadj := -log10(padj)] %>%
  .[, sig := padj < opts$p_cutoff]

# labs <- paste0(c("q < ", "q >= "), opts$p_cutoff)
labs <- c("NS", "Significant")

p <- ggplot(cors, aes(r, logpadj, colour = sig, label = gene)) +
  geom_point() +
  facet_wrap(~anno) +
  scale_colour_manual(values = c("grey", "navy"), labels = labs, name = NULL) +
  geom_hline(yintercept = -log10(opts$p_cutoff), colour = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = mean(r)), colour = "blue") +
  #geom_text_repel(data = cors[sig == TRUE]) +
  labs(x = c("Weighted Pearson R"), y = "-log10 q-value") +
  guides(label = FALSE) +
  theme_bw()
p

if(!exists(io$plot_dir)){
  dir.create(io$plot_dir)
}

save_plot(paste(io$plot_dir, "met_rna_correlations.pdf", sep="/"), p, base_width = 12, base_height = 12)

fwrite(cors, paste(io$plot_dir, "met_rna_correlations.tsv", sep="/"), sep="\t", 
       col.names = T, row.names = F, quote=F, na="NA")
