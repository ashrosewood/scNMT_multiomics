args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--model      : trained model                       [required]\n")
    cat("--data       : data matrix list                    [required]\n")
    cat("--plotdir    : output directory for plots          [required]\n")
    cat("--clusters   : number of desired clusters          [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$model      <- sub( '--model=', '', args[grep('--model', args)] )
    io$data       <- sub( '--data=', '', args[grep('--data', args)] )
    io$plotdir   <- sub( '--plotdir=', '', args[grep('--plotdir', args)] )
    opts$clusters <- sub( '--clusters=', '', args[grep('--clusters', args)] )
}

library(MultiAssayExperiment)
library(MOFA2)
library(ggplot2)
library(rhdf5)
library(stringr)
library(umap)
library(Rtsne)
library(irlba)
library(data.table)
library(png)
library(Seurat)
library(purrr)
library(ggpubr)

### TEST INPUT ###
#io$data <- "data/all_matrix_list.rds"
#io$model <- "data/hdf5/model_trained.hdf5"
#io$plotdir <- "plots"
#opts$clusters <- 3

my_data <- readRDS(io$data)

for (i in seq(1:length(my_data))) {
    my_data[[i]] <- my_data[[i]][(which(rowSums(my_data[[i]], na.rm = T) != 0)),]
}

good_cells <- intersect(intersect(colnames(my_data$rna),colnames(my_data$met_body)), colnames(my_data$acc_body))

# Create MOFAobject
MOFAobject <- create_mofa(list("RNA" = my_data$rna[,good_cells], "met_body" = my_data$met_body[,good_cells], "met_Enhancer" = my_data$met_Enhancer[,good_cells], "met_H3K27ac" = my_data$met_MCF7_H3K27ac_peaks[,good_cells], "acc_body" = my_data$acc_body[,good_cells], "acc_Enhancer" = my_data$acc_Enhancer[,good_cells], "acc_H3K27ac" = my_data$acc_MCF7_H3K27ac_peaks[,good_cells]))

MOFAobject <- load_model(io$model)

sample_metadata <- data.frame(sample = samples_names(MOFAobject)[[1]],
                              condition = substr(samples_names(MOFAobject)[[1]],1,3)
                              )

samples_metadata(MOFAobject) <- sample_metadata

plot_data_overview(MOFAobject)

plot_factor_cor(MOFAobject, cluster = T)

png(paste0(io$plotdir,"/full_variance.png"))
plot_variance_explained(MOFAobject, x="factor", y="view", legend = T)
dev.off()

plot_top_weights(MOFAobject, view = "RNA", factor = 1, nfeatures = 10, scale = T, abs=T)
plot_top_weights(MOFAobject, view = "acc_body", factor = 1, nfeatures = 10, scale = T, abs=T)

plot_factor(MOFAobject, factor = 1, scale = T, add_violin = T, dodge = T, dot_size = 1, legend = T)

plot_data_heatmap(MOFAobject, view = "RNA", factor = 1, features = 20, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE)

plot_data_scatter(MOFAobject, view = "RNA", factor = 1, features = 5, add_lm = T, color_by = "condition")

set.seed(42)

model <- run_tsne(MOFAobject)

plot_dimred(model, method = "TSNE", color_by = "condition")

################### old MOFA ####################################

'''
variance_grid <- function (object, cluster = TRUE, ...) 
{
    R2_list <- calculateVarianceExplained(object, ...)
    fvar_m <- R2_list$R2Total
    fvar_mk <- R2_list$R2PerFactor
    fvar_mk_df <- reshape2::melt(fvar_mk, varnames = c("factor", 
        "view"))
    fvar_mk_df$factor <- factor(fvar_mk_df$factor)
    if (cluster & ncol(fvar_mk) > 1) {
        hc <- hclust(dist(t(fvar_mk)))
        fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk)[hc$order])
    }
    hm <- ggplot(fvar_mk_df, aes_string(x = "view", y = "factor")) + 
        geom_tile(aes_string(fill = "value"), color = "black") + 
        guides(fill = guide_colorbar("R2")) + scale_fill_gradientn(colors = c("gray97", 
        "darkblue"), guide = "colorbar") + ylab("Latent factor") + 
        theme(plot.title = element_text(size = 17, hjust = 0.5), 
            axis.title.x = element_blank(), axis.text.x = element_text(size = 11, 
                angle = 60, hjust = 1, vjust = 1, color = "black"), 
            axis.text.y = element_text(size = 12, color = "black"), 
            axis.title.y = element_text(size = 15), axis.line = element_blank(), 
            axis.ticks = element_blank(), panel.background = element_blank())
    hm <- hm + ggtitle("Variance explained per factor") + 
        guides(fill = guide_colorbar("R2"))
    hm
}

png(paste0(io$plotdir,"/factor_variance.png"))
variance_grid(MOFAobject)
dev.off()

png(paste0(io$plotdir,"/topweights_F1.png"))
plotTopWeights(
  object = MOFAobject,
  view = "rna", 
  factor = 1, 
  nfeatures = 10
)
dev.off()

png(paste0(io$plotdir,"/topweights_F2.png"))
plotTopWeights(
  object = MOFAobject,
  view = "rna", 
  factor = 2, 
  nfeatures = 10
)
dev.off()


png(paste0(io$plotdir,"/weights_HM_F1.png"))
my_factor1 <- sort(getFactors(MOFAobject,"LF1")[,1])
my_order_samples <- names(my_factor1)
my_df <- data.frame(
  row.names = my_order_samples,
#  culture = getCovariates(MOFAobject, "culture")[my_order_samples],
  factor = my_factor1
)

plotDataHeatmap(
  object = MOFAobject, 
  view = "rna", 
  factor = "LF1", 
  features = 20, 
  transpose = TRUE, 
  show_colnames=FALSE, show_rownames=TRUE, # pheatmap options
  cluster_cols = FALSE, annotation_col=my_df # pheatmap options
)
dev.off()

png(paste0(io$plotdir,"/weights_HM_F2.png"))
my_factor1 <- sort(getFactors(MOFAobject,"LF2")[,1])
my_order_samples <- names(my_factor1)
my_df <- data.frame(
  row.names = my_order_samples,
#  culture = getCovariates(MOFAobject, "culture")[my_order_samples],
  factor = my_factor1
)

plotDataHeatmap(
  object = MOFAobject, 
  view = "rna", 
  factor = "LF2", 
  features = 20, 
  transpose = TRUE, 
  show_colnames=FALSE, show_rownames=TRUE, # pheatmap options
  cluster_cols = FALSE, annotation_col=my_df # pheatmap options
)
dev.off()

clusters <- clusterSamples(MOFAobject, k=opts$clusters)

theme_pub <- function() {
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.2)),
    axis.title = element_text(size=rel(1.2), color="black"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

factors <- seq(1:2)
Z <- getFactors(MOFAobject) %>% .[,factors]

# Scale Z by the variance explained
# r2 <- apply(calculateVarianceExplained(model, factors=factors)$R2PerFactor,1,sum)
# Z <- Z %*% diag(r2*100)
algorithms <- c("umap")
umap.defaults$n_neighbors <- 25
umap.defaults$min_dist <- 0.55

for (algorithm in algorithms) {

  set.seed(42)
  if (algorithm=="tsne") {
    tsne <- Rtsne(Z, check_duplicates=FALSE, pca=FALSE, theta=0.5, dims=2)
    Z.out <- tsne$Y
  } else if (algorithm=="umap") {
    umap.out <- umap(Z, config = umap.defaults)
    Z.out <- umap.out$layout
  }

  # Flip a factor
  Z.out[,2] <- -Z.out[,2]

  to.plot <- Z.out %>% as.data.table %>% .[,sample:=rownames(Z)]
  to.plot$pclusters <- clusters

  p1 <- ggplot(to.plot, aes(x=V1, y=V2, color=pclusters)) +
    geom_point(alpha=0.7, size=2.0) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    theme_pub() + theme(legend.position = "none")

  png(paste0(io$plotdir,"/MOFA_UMAP.png"))
  print(p1)
  dev.off()

  # Save coordinates
  # fwrite(to.plot[,c("sample","V1","V2")], sprintf("%s/%s_coordinates.txt",io$outdir,algorithm))
}

'''
