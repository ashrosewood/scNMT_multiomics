
args <- commandArgs()

help <- function(){
    cat("plot_difmet_difacc.R :
Plots the differential methylation and differential accessibility analysis done in the prior scripts\n")
    cat("Usage: \n")
    cat("--met       : differential methylation table   [required]\n")
    cat("--acc       : differential accessibility table [required]\n")
    cat("\n")
    q()
}

io <- list()
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$difmet   <- sub( '--met=', '', args[grep('--met', args)] )
    io$difacc   <- sub( '--acc=', '', args[grep('--acc', args)] )
}

if ("nVennR" %in% rownames(installed.packages()) == FALSE) {
    install.packages("nVennR", repos="https://ftp.osuosl.org/pub/cran/")
}

library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(nVennR)


### TEST INPUT ###   
#io <- list()
#io$difmet <- "tables/difmet.tsv"
#io$difacc <- "tables/difacc.tsv"
io$data.dir <- "tables"


io$difrna <- "../scNMT_transcriptomeMapping/data/seurat/post_DEgenes.tsv"
#io$difrna <- "../scNMT_transcriptomeMapping/data/seurat/post_DEgenes_CCreduced.tsv"


rnadif <- read.table(io$difrna, header = T, sep = "\t")

metdata <- dir(io$data.dir, pattern = "difmet", full = TRUE) %>%
    map(fread) %>%
    map(~.[, .(id, gene, anno, padj_fdr, log_padj_fdr, diff, sig)])
names(metdata) <- sapply(strsplit(dir(io$data.dir, pattern = "difmet", full = TRUE), "\\/|\\.| "), "[", 2)


accdata <- dir(io$data.dir, pattern = "difacc", full = TRUE) %>%
    map(fread) %>%
    map(~.[, .(id, gene, anno, padj_fdr, log_padj_fdr, diff, sig)])
names(accdata) <- sapply(strsplit(dir(io$data.dir, pattern = "difacc", full = TRUE), "\\/|\\.| "), "[", 2)

for (i in seq(1:length(metdata))) {

    dif <- list(metdata[[i]],accdata[[i]])
    names(dif) <- c(names(metdata)[i], names(accdata)[i])

    dif_merge <- map2(copy(dif), names(dif),  ~setnames(.x, 
                                              c("padj_fdr", 
                                                "log_padj_fdr", 
                                                "diff", 
                                                "sig"), 
                                              c(paste0("p_", .y), 
                                                paste0("logp_", .y),
                                                paste0("dif_", .y), 
                                                paste0("sig_", .y)))) %>%
        purrr::reduce(merge, by = c("id", "anno", "gene"))

# volcano plots #

    volc_plot <- function(dt, nlabs = 10L, title = NULL){
  
        labs <- dt[sig==TRUE][order(padj_fdr)]
        nlabs <- min(nlabs, nrow(labs))
        labs <- labs[1:nlabs]
  
        ggplot(dt, aes(diff, log_padj_fdr)) +
            geom_point(alpha = 0.75, aes(colour = sig)) +
            scale_colour_manual(values = c("black", "red")) +
            guides(colour = FALSE) +
            theme_cowplot() +
            xlab("% difference") +
            ylab("log10 q-value") +
            geom_text_repel(data = labs, aes(label = gene)) +
            ggtitle(title)
    }

    volc <- map(dif, ~ggplot(., aes(diff, log_padj_fdr, colour = sig)) +
              geom_point() +
                  facet_wrap(~anno))

    ending <- substr(names(dif)[1], 7, 12)

    png(paste0("plots/acc_volcano_", ending, ".png"))
    print(volc[2])
    dev.off()

    png(paste0("plots/met_volcano_", ending, ".png"))
    print(volc[1])
    dev.off()

# dif Scatter #
    dif_merge$dif_met <- dif_merge[,6]
    dif_merge$dif_acc <- dif_merge[,10]
    dif_merge$sig <- dif_merge[,7]*dif_merge[,11]
    sp <- ggplot(dif_merge, aes(dif_met, dif_acc, colour = sig)) +
        geom_point(alpha = 0.2)

    png(paste0("plots/dif_scatter_", ending, ".png"))
    print(sp)
    dev.off()


# dif Venn Diagram #

    rna <- unique(rnadif$gene)
    met <- subset(dif$met, dif$met$sig == T)
    met <- unique(met$gene)
    acc <- subset(dif$acc, dif$acc$sig == T)
    acc <- unique(acc$gene)

    rnamet <- intersect(rna, met)
    accrna <- intersect(rna, acc)
    accmet <- intersect(acc, met)

    accmetrna <- intersect(rnamet, acc)

    outfile <- paste0("plots/dif_venn_", ending, ".svg")
        
    myV <- createVennObj(nSets = 3, sNames = c("RNA","methylation","accessibility"), sSizes = c(0, length(acc)-length(accmet)-length(accrna)-length(accmetrna), length(met)-length(accmet)-length(rnamet)-length(accmetrna), length(accmet)-length(accmetrna), length(rna)-length(rnamet)-length(accrna)-length(accmetrna), length(accrna)-length(accmetrna), length(rnamet)-length(accmetrna), length(accmetrna)))

    myV <- plotVenn(nVennObj = myV, borderWidth = 2, setColors = c('blue', 'red', 'green'), outFile = outfile)

}
