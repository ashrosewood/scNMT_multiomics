
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

io$difrna <- "../scNMT_transcriptomeMapping/tables/DE_genes.tsv"

dif <- map(list(met = io$difmet, acc = io$difacc), fread) %>% 
  map(~.[, .(id, gene, anno, padj_fdr, log_padj_fdr, diff, sig)]) 

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

rnadif <- read.table(io$difrna, header = T, sep = "\t")

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

png("plots/acc_volcano.png")
volc$acc
dev.off()

png("plots/met_volcano.png")
volc$met
dev.off()

# dif scatter plot #
sp <- ggplot(dif_merge, aes(dif_met, dif_acc, colour = sig_met*sig_acc)) +
  geom_point(alpha = 0.2)

png("plots/dif_scatter.png")
sp
dev.off()


# dif Venn Diagram #

rna <- unique(rnadif$gene)
met <- subset(dif$met, sig == T)
met <- unique(met$gene)
acc <- subset(dif$acc, sig == T)
acc <- unique(acc$gene)

rnamet <- intersect(rna, met)
accrna <- intersect(rna, acc)
accmet <- intersect(acc, met)

accmetrna <- intersect(rnamet, acc)

myV <- createVennObj(nSets = 3, sNames = c("RNA","methylation","accessibility"), sSizes = c(0, length(acc)-length(accmet)-length(accrna)-length(accmetrna), length(met)-length(accmet)-length(rnamet)-length(accmetrna), length(accmet)-length(accmetrna), length(rna)-length(rnamet)-length(accrna)-length(accmetrna), length(accrna)-length(accmetrna), length(rnamet)-length(accmetrna), length(accmetrna)))

myV <- plotVenn(nVennObj = myV, borderWidth = 2, setColors = c('blue', 'red', 'green'), outFile = 'plots/dif_venn.svg')
