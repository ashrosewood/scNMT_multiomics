args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--met      : methylation directory                          [required]\n")
    cat("--acc      : accessibility directory                        [required]\n")
    cat("--rna      : Seurat object with RNA info                    [required]\n")
    cat("--qcinfo   : Sample metadata file with QC info              [required]\n")
    cat("--genemeta : Gene metadata file                             [required]\n")
    cat("--anno     : Directory with annnotation bed files           [required]\n")
    cat("--groups   : File containing the transcriptome cluster info [required]\n")
    cat("--outdir   : output directory                               [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$met_dir         <- sub( '--met=', '', args[grep('--met', args)] )
    io$acc_dir         <- sub( '--acc=', '', args[grep('--acc', args)] )
    io$rna_file        <- sub( '--rna=', '', args[grep('--rna', args)] )
    io$sample.metadata <- sub( '--qcinfo=', '', args[grep('--qcinfo', args)] )
    io$annos_dir       <- sub( '--anno=', '', args[grep('--anno', args)] )
    io$gene.metadata   <- sub( '--genemeta=', '', args[grep('--genemeta', args)] )
    io$groups          <- sub( '--groups=', '', args[grep('--groups', args)] )
    io$outfile         <- sub( '--outdir=', '', args[grep('--outdir', args)] )
}



library(data.table)
library(purrr)
library(Seurat)
library(ggplot2)
library(cowplot)

### functions ###
fread_gz = function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") return(fread(cmd = paste("zcat <", filename), ...))
  fread(filename, ...)
}

### TEST INPUT ###

## Define I/O ##
#io <- list()
#io$sample.metadata <- "../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt" 
#io$annos_dir  <- "data/anno"
#io$gene.metadata <- "../scNMT_transcriptomeMapping/data/gene_metadata.tsv"
#io$groups <- "data/groups.tsv"
#io$met_dir <- "../scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
#io$acc_dir <- "../scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised"
#io$rna_file <- "../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds"
#io$outfile <- "plots/zoom"
io$rna_de <- "../scNMT_transcriptomeMapping/data/seurat/post_DEgenes_CCreduced.tsv"
io$met_de <- "tables/difmetA_vs_B.tsv"
io$acc_de <- "tables/difaccA_vs_B.tsv"

## Options ##
opts <- list()
opts$gene <- "PDZK1"
#opts$gene <- "COASY"
opts$window <- 10000
opts$slide <- 1000
opts$up    <- 20000
opts$down  <- 20000
opts$features <- c("promoters", 
                   "T47D_ER_peaks", "T47D_H3K27ac_peaks")


dir.create(io$outfile, recursive = TRUE)

rna_de <- read.table(io$rna_de, header=T)
met_de <- read.table(io$met_de, header=T)
acc_de <- read.table(io$acc_de, header=T)

met_de <- subset(met_de, sig == TRUE & padj_fdr < .01)
acc_de <- subset(acc_de, sig == TRUE & padj_fdr < .01)

overlap <- intersect(rna_de$gene,intersect(met_de$gene,acc_de$gene))

meta <- fread(io$sample.metadata) %>%
  .[pass_metQC == TRUE & pass_accQC == TRUE & pass_CHGQC == TRUE & pass_CHHQC == TRUE]

cells <- meta[, unique(sample)] %>%
  gsub("sc_", "", .)

gene <- fread(io$gene.metadata) %>%
  .[gene == opts$gene] %>%
  .[, c("start", 
        "end", 
        "chr") := .(start - opts$up, 
                    end + opts$down,
                    chr)]

gene_win <- gene[, .(gene = gene, chr = chr, start = seq(start, end, by = opts$slide))] %>%
  .[, end := start + opts$window] %>%
  setkey(chr, start, end)

seurat <- readRDS(io$rna_file)
#test <- subset(seurat, features = opts$gene)
rna <- data.table(sample = colnames(seurat), 
                  expr = as.vector(seurat@assays$RNA@data[opts$gene, ]))

groups <- fread(io$groups) %>%
  .[, .(sample = gsub("sc_", "", id_rna), group = as.factor(group))]

groups$sample <- sub("D","_", groups$sample)

metacc <- list(acc = paste0(io$acc_dir, "/", cells, "_GpC.gz"),
               met = paste0(io$met_dir, "/", cells, "_CpG.gz")) %>%
  map(~.[file.exists(.)]) %>%
  map(~map(., ~fread_gz(.) %>%
             
            .[, .(chr, start = pos, end = pos, rate)] %>%
            setkey(chr, start, end) %>%
            foverlaps(gene_win, nomatch=0L)) %>%
        
        map2(cells, ~.x[, sample := .y]) %>%
        rbindlist() %>% 
        .[, .(rate = mean(rate)), .(sample, chr, gene, start, end)] %>%
        merge(groups, by = "sample")
            
          )

# find overlapping annotations

anno <- paste0(io$annos_dir, "/", opts$features, ".bed") %>%
  .[file.exists(.)] %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "st", "id", "anno")) %>%
  .[, chr := gsub("chr", "", chr) %>% paste0("chr", .)] %>%
  setkey(chr, start, end) %>%
  foverlaps(gene %>% setkey(chr, start, end), nomatch = 0L) %>%
  .[, .(i.start, i.end, anno = as.factor(anno))] %>%
  unique()


# compute means, sds and cors

means <- map(metacc, ~.[, .(mean = mean(rate), sd = sd(rate), .N), .(start, end, group)])

rna$sample <- sub("D","_", rna$sample)

cors <- map2(metacc, names(metacc), ~.x[, omic := .y]) %>%
  rbindlist() %>%
  dcast(sample + start + end ~ omic, value.var = "rate") %>%
  merge(rna, by = "sample", allow.cartesian = TRUE) %>%
  .[, .(metacc = cor(met, acc, use = "pair"),
        metrna = cor(met, expr, use = "pair"),
        accrna = cor(acc, expr, use = "pair"))
    , .(start, end)] %>%
  melt(id.vars = c("start", "end"), variable.name = "comparison", value.name = "r") %>%
  .[, comparison := factor(comparison, levels = c("metrna", "metacc","accrna"))]

# make plots

mean_plot <- map2(means, names(means), 
                  ~{
                    if (.y == "acc") {
                      ylab <- "Accessibility"
                      colours <- c("dodgerblue1", "dodgerblue4", "blue")
                    } else {
                      ylab <- "Methylation"
                      colours <- c("orangered1", "orangered4", "red")
                    }
                    ylab <- ifelse(.y == "acc", "Accessibility", "Methylation")
                    ggplot(.x, aes(start, mean, colour = group, fill = group)) +
                      geom_line(size = 1) +
                      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), 
                                  alpha = 0.1,
                                  linetype = 0) +
                      # ylim(0, NA) +
                      theme_cowplot() +
                      ylab(ylab) + 
                      xlab("") +
                      scale_colour_manual(values = colours) +
                      scale_fill_manual(values = colours) +
                      theme(line = element_blank(), legend.title = element_blank()) +
                      theme(axis.text.x = element_blank())
                    
                   })

# mean_plot$met <- mean_plot$met + xlab(paste0("Genomic position (", gene[,chr], ")")) 
# mean_plot$acc <- mean_plot$acc + theme(axis.text.x = element_blank())

cor_plot <- ggplot(cors, aes(start, r, fill = comparison, colour = comparison)) +
  geom_line(size = 1) +
  theme_cowplot() +
  xlab("") +
  ylab("Correlation") +
  theme(line = element_blank(), legend.title = element_blank()) +
  theme(axis.text.x = element_blank())

gene.anno <- fread(io$gene.metadata) %>%
  .[gene == opts$gene] %>%
  .[, c("i.start", 
        "i.end", 
        "anno") := .(start, 
                    end, 
                    gene)] %>%
  .[,.(i.start, i.end, anno)]

boo <- rbind(gene.anno, anno)

anno_plot <- ggplot(boo) +
  geom_rect(aes(xmin = i.start, 
                xmax = i.end,
                # ymin = 0, ymax = 1,
                ymin=as.numeric(anno)-0.20,
                ymax=as.numeric(anno)+0.20,
                fill=anno), 
            alpha=0.5) +
  theme_cowplot() +
  xlab(paste0("Genomic position (", gene[,chr], ")")) +
  theme(line = element_blank(), legend.title = element_blank()) +
  theme(axis.text.y = element_blank())

final_plot <- cowplot::plot_grid(mean_plot$acc, 
                                 mean_plot$met, 
                                 cor_plot, 
                                 anno_plot,
                                 align = "v", 
                                 ncol = 1, 
                                 rel_heights=c(1/6, 1/3, 1/3, 1/6))


out_file <- paste0(io$outfile, "/", opts$gene, "zoom_plot", ".pdf")
save_plot(out_file, final_plot, base_height = 12, base_width = 8)
