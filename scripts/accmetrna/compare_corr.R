args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--metrna      : methylation/RNA correlations tsv                       [required]\n")
    cat("--accrna      : accessibility/RNA correlations tsv                     [required]\n")
    cat("--plotdir      : output directory                                       [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$rna_met    <- sub( '--metrna=', '', args[grep('--metrna', args)] )
    io$rna_acc    <- sub( '--accrna=', '', args[grep('--accrna', args)] )
    io$plot_dir   <- sub( '--plotdir=', '', args[grep('--plotdir', args)] )
}


library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)
library(ggrepel)


### TEST INPUT ###
#io$rna_met  <- "plots/cor_met/met_rna_correlations.tsv"
#io$rna_acc  <- "plots/cor_acc/acc_rna_correlations.tsv"
#io$plot_dir <- "plots/cor_accmetrna"

opts$p_cutoff <- 0.1

met <- fread(io$rna_met)
acc <- fread(io$rna_acc)

df <- merge(met, acc, by=c("anno", "id", "gene", "ens_id.x"), 
            suffixes=c(".met", ".acc"))

df$sig <- FALSE

df[r.met >= 0.5 & r.acc > 0.5 | r.met >= 0.5 & r.acc <= -0.5, sig:=TRUE]
df[r.met <= -0.5 & r.acc < -0.5 | r.met <= -0.5 & r.acc >= 0.5, sig:=TRUE]

labs <- c("NS", "Significant")


p <- ggplot(df, aes(r.met, r.acc, colour = sig)) +
    geom_point() +
  facet_wrap(~anno) +
  xlim(-1,1) +
  ylim(-1,1) +
                                        #  scale_colour_manual(values = c("grey", "navy"), labels = labs, name = NULL) +
      scale_colour_manual(values = c("grey", "navy")) +
  geom_hline(yintercept = 0, colour = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "blue", linetype = "dashed") +
#  geom_text_repel(data = df[sig == TRUE]) +
  labs(x = c("Met-vs-RNA Pearson R"), y = "ACC-vs-RNA Pearson R") +
  guides(label = FALSE) +
  theme_bw()
p

if(!exists(io$plot_dir)){
  dir.create(io$plot_dir)
}

save_plot(paste(io$plot_dir, "accmetrna_correlations.scatter.pdf", sep="/"), p, base_width = 6, base_height = 6)

fwrite(df
       , paste(io$plot_dir, "accmetrna_correlations.tsv", sep="/"), sep="\t" 
      , col.names = T, row.names = F, quote=F, na="NA")


acc[,Type := "acc"]
met[,Type := "met"]

merg <- rbind(acc, met)

merg <- subset(merg, sig == TRUE)

p <-
  ggplot(merg, aes(x=anno, y=r, fill=Type)) +
  geom_boxplot(aes(fill = Type), alpha=1.0, outlier.shape = NA) +
  #geom_jitter(alpha=0.5, color=c("#00136C", "#F87D42")) +
  #geom_jitter(aes(color = Type), alpha=0.5) +
  scale_color_manual(values=c("#00136C", "#F87D42"))+
 # scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"),
  #                  labels=c("CG methylation","GC accessibility"))+
  scale_fill_manual("legend", values = c("acc" = "#F87D42", "met" = "#00136C"))+
  ylab("Correlations") +
  xlab("") +
  ggtitle("Compare correlations relative to RNA") +
  coord_flip()+
  theme_bw() +
  theme(
    axis.title.y = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.title.x = element_text(colour="black", size=17, margin=margin(0,20,0,0)),
    axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.text.y = element_text(colour="black", size=11),
    axis.ticks = element_line(colour="black"),
    legend.title= element_text(size=15),
    legend.text = element_text(size=15)
  )
p

if(!exists(io$plot_dir)){
  dir.create(io$plot_dir)
}

save_plot(paste(io$plot_dir, "accmetrna_correlations_significant.boxplot.pdf", sep="/"), p, base_width = 6, base_height = 6)
