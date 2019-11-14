args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--met      : methylation directory                       [required]\n")
    cat("--acc      : accessibility directory                     [required]\n")
    cat("--rna      : Seurat object with RNA info                 [required]\n")
    cat("--qcinfo   : Sample metadata file with QC info           [required]\n")
    cat("--genemeta : Gene metadata file                          [required]\n")
    cat("--anno     : Directory with annnotation bed files        [required]\n")
    cat("--outdir   : output directory                            [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$met.dir         <- sub( '--met=', '', args[grep('--met', args)] )
    io$acc.dir         <- sub( '--acc=', '', args[grep('--acc', args)] )
    io$rna.file        <- sub( '--rna=', '', args[grep('--rna', args)] )
    io$sample.metadata <- sub( '--qcinfo=', '', args[grep('--qcinfo', args)] )
    io$annos_dir       <- sub( '--anno=', '', args[grep('--anno', args)] )
    io$gene_metadata   <- sub( '--genemeta=', '', args[grep('--genemeta', args)] )
    io$outdir          <- sub( '--outdir=', '', args[grep('--outdir', args)] )
}



library(data.table)
library(dplyr)
library(Seurat)
library(purrr)
library(SummarizedExperiment)


opts$met.annoname <- c(
    "CGI_promoters1000",
    "MCF7_ChromHMM_Enhancer",
    "MCF7_H3K27ac_peaks",
    "nonCGI_promoters1000"
)

opts$acc.annoname <- c(
    "CGI_promoters1000",
    "MCF7_ChromHMM_Enhancer",
    "MCF7_H3K27ac_peaks",
    "nonCGI_promoters1000"
)

opts$met.annos <- c(
    "CGI_promoter",
    "Enhancer",
    "MCF7_H3K27ac_peaks",
    "nonCGI_promoter"
)

opts$acc.annos <- c(
    "CGI_promoter",
    "Enhancer",
    "MCF7_H3K27ac_peaks",
    "nonCGI_promoter"
)

####### TEST INPUT #########
#io$met.dir  <- "../scNMT_NOMeWorkFlow/data/met"
#io$acc.dir  <- "../scNMT_NOMeWorkFlow/data/acc"
#io$rna.file <- "../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds"

#io$sample.metadata <- "../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt"
#io$annos_dir       <- "data/anno"
#io$gene_metadata   <- "../scNMT_transcriptomeMapping/data/gene_metadata.tsv"

#io$outdir <- "data"

#io$met.stats <- "/met/stats/samples/sample_stats.txt"
#io$acc.stats <- "/acc/stats/samples/sample_stats.txt"

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

tmp <- fread(io$sample.metadata)
opts$met_cells <- tmp %>% .[pass_metQC==T, id]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id]

opts$overlapGenes  <- TRUE
opts$gene_window  <- 5e4

opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature
opts$acc_nfeatures <- 1000

opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

###############
## Load data ##
###############

# Load Methylation data
met_dt <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[sample%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id","anno","rate","Nmet","N","sample", "var"))

# Load Accessibility data
acc_dt <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[sample%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id","anno","rate","Nmet","N","sample", "var"))

# Load RNA data
sce <- readRDS(io$rna.file) # %>% .[,opts$rna_cells]

# Load annotation metadata
feature_metadata <- lapply(unique(c(opts$met.annoname, opts$acc.annoname)), function(i) 
  fread(sprintf("%s/%s.bed",io$annos_dir,i), stringsAsFactors=T)[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

# Load gene metadata 
gene_metadata <- fread(io$gene_metadata,stringsAsFactors=T) %>% 
  setnames(c("ens_id","gene"),c("id","gene")) %>% 
  .[,chr:=as.factor(sub("chr","",chr))]

temp <- as.data.frame(gene_metadata)
rownames(temp) <- temp$id

################
## Parse data ##
################

# Parse gene and feature metadata
feature_metadata_filt.met <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% met_dt[anno==y,id]] ) %>%
  rbindlist
feature_metadata_filt.acc <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% acc_dt[anno==y,id]] ) %>%
  rbindlist

gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","gene")] %>% 
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)

## Parse RNA expression data ##

# Convert to data.table
rna_dt <- as.matrix(sce@assays[["RNA"]]@data) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "gene") %>%
  merge(temp %>% tibble::rownames_to_column("ens_id") %>% .[,c("gene","ens_id")], by = "gene")
# rna_dt[,c("id_rna","gene","ens_id"):=list(as.factor(id_rna),as.factor(gene),as.factor(ens_id))]

## Parse accessibility data ##

acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


##############################
## Merge data with metadata ##
##############################
sample_metadata <- fread(io$sample.metadata) %>%
  .[id%in%opts$met_cells | id %in% opts$acc_cells ]


#acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
#met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()
#rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage","stage_lineage")], by="id_rna") %>% droplevels()


###########################################################
## Associate the genomic features with overlapping genes ##
###########################################################

# Methylation
if (opts$overlapGenes) {
  met_list <- list()
  for (i in unique(met_dt$anno)){
    
    # Subset corresponding anno
    met_tmp <- met_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSG", unique(met_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.met[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      met_list[[i]] <- merge(met_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSG", unique(met_tmp$id)) == TRUE)) {
      met_list[[i]] <- merge(met_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  met_dt <- rbindlist(met_list)
  rm(met_list, met_tmp,feature_metadata_tmp,ov)
} else {
  met_dt[,gene:="NA"]
}

# Accessibility
if (opts$overlapGenes) {
  acc_list <- list()
  for (i in unique(acc_dt$anno)){
    
    # Subset corresponding anno
    acc_tmp <- acc_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSG", unique(acc_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.acc[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      acc_list[[i]] <- merge(acc_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSG", unique(acc_tmp$id)) == TRUE)) {
      acc_list[[i]] <- merge(acc_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  acc_dt <- rbindlist(acc_list)
  rm(acc_list, acc_tmp,feature_metadata_tmp,ov)
} else {
  acc_dt[,gene:="NA"]
}

#############################
## Add Gene names to dt    ##
#############################

#met_dt$gene <- gene_metadata[met_dt$id,3]
#acc_dt <- merge(acc_dt, gene_metadata, by = "id")
#rna_dt <- merge(rna_dt, gene_metadata, by = "id")


#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt$var[is.na(met_dt$var)] <- 0
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by  minimum number of cells
met_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by variance
met_dt <- met_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0.1] %>% .[,var:=NULL] %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt$var[is.na(acc_dt$var)] <- 0
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by variance
acc_dt <- acc_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0.1] %>% .[,var:=NULL] %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Remove lowly expressed genes
rna_dt <- rna_dt[,mean:=mean(expr),by="ens_id"] %>% .[mean>0.1] %>% .[,mean:=NULL]

# Remove genes with constant expression levels
rna_dt <- rna_dt[,var:=var(expr),by="ens_id"] %>% .[var>0.1] %>% .[,var:=NULL]

# Filter genes with low cellular detection rate and sites with low coverage across samples
rna_dt <- rna_dt[,cdr:=sum(expr>0)/length(unique(rna_dt$id_rna)), by="ens_id"] %>% .[cdr>=opts$rna_min.cdr] %>% .[,cdr:=NULL]

# Extract top N highly variable genes
rna_dt <- rna_dt[,var:=var(expr), by="ens_id"] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

############################
## Regress out covariates ##
############################

# RNA: number of expressed genes
foo <- data.table(id_rna=colnames(sce), covariate=sce$nFeature_RNA/nrow(sce))
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]], by=c("gene")] %>%
  .[,covariate:=NULL]


# Methylation: differences in mean methylation rate
met.stats <- fread(io$sample.metadata)
met.stats <- subset(met.stats, context == "GC")
foo <- met_dt[,.(covariate=mean(m)),by=c("id")]
foo <- met.stats %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("sample","mean")]
met_dt <- met_dt %>% merge(foo, by="sample") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

# Accessibility: differences in global accessibility rate
acc.stats <- fread(io$sample.metadata)
acc.stats <- subset(acc.stats, context == "CG")
foo <- acc.stats %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("sample","mean")]
acc_dt <- acc_dt %>% merge(foo, by="sample") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

#####################################
## Select highly variable features ##
#####################################

# RNA: Extract top N highly variable genes
keep_hv_genes <- rna_dt[,.(var=var(expr)), by="ens_id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id
rna_dt <- rna_dt[ens_id%in%as.character(keep_hv_genes)] %>% droplevels()

# Accessibility: Extract top N most variable features
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

# Methylation: Extract top N most variable features
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

met_dt[,var(m),by=c("id","anno")]

#############################
## Create joint data.frame ##
#############################

data1 <- rna_dt %>% .[,c("id_rna","ens_id","expr")] %>%  
  setnames(c("sample","feature","value")) %>% .[,c("feature_group","sample_group"):=list("RNA","MCF7")]
data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("met_",feature), paste0("met_",feature_group), "MCF7")]
data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group), "MCF7")]

temp2 <- gsub("(B)(10)1(_[A-Z0-9]+)","\\1C\\2\\3", data1$sample)

data1$sample <- temp2

all_rep <- Reduce(intersect,list(data1$sample,data2$sample,data3$sample))

test1 <- subset(data1, sample %in% all_rep)
test2 <- subset(data2, sample %in% all_rep)
test3 <- subset(data3, sample %in% all_rep)

data <- rbind(test1,test2,test3)

outfile <- paste0(io$outdir,"/data.txt")
fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
system(sprintf("pigz -f %s",outfile))


###########################
## Create input matrices ##
###########################

met_cells <- as.character(unique(data2$sample))
rna_cells <- as.character(unique(data1$sample))
acc_cells <- as.character(unique(data3$sample))

rna_dt$sample <- temp2

met_anno_list <- c("CGI_promoter", "Enhancer", "MCF7_H3K27ac_peaks", "nonCGI_promoter")
acc_anno_list <- c("CGI_promoter", "Enhancer",  "MCF7_H3K27ac_peaks", "nonCGI_promoter")

rna_matrix <- rna_dt[,c("gene","expr","sample")] %>%
  .[,c("sample","gene"):=list(as.character(sample),as.character(gene))] %>%
  .[,sample:=factor(sample,levels=Reduce(intersect,list(rna_cells,acc_cells,met_cells)))] %>%
  dcast(sample~gene, value.var="expr", drop=F, fun.aggregate=sum) %>% matrix.please() %>% t

met_matrix_list <- list()
for (n in unique(met_dt$anno)) {
  met_matrix_list[[paste("met",n,sep="_")]] <- met_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(intersect,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F, fun.aggregate=sum) %>% matrix.please() %>% t
  
  cat(sprintf("%s methylation matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(met_matrix_list[[paste("met",n,sep="_")]]), ncol(met_matrix_list[[paste("met",n,sep="_")]]),
              100*mean(is.na(met_matrix_list[[paste("met",n,sep="_")]]))))
}

cat("\n")

acc_matrix_list <- list()
for (n in unique(acc_dt$anno)) {
  acc_matrix_list[[paste("acc",n,sep="_")]] <- acc_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(intersect,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F, fun.aggregate=sum) %>% matrix.please() %>% t
  
  cat(sprintf("%s accessibility matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(acc_matrix_list[[paste("acc",n,sep="_")]]), ncol(acc_matrix_list[[paste("acc",n,sep="_")]]),
              100*mean(is.na(acc_matrix_list[[paste("acc",n,sep="_")]]))))
}

all_matrix_list <- c(rna=list(rna_matrix),met_matrix_list,acc_matrix_list)

saveRDS(all_matrix_list, "data/all_matrix_list.rds")

###### SMALLER DATASET ###########
'''
met_matrix_list <- list()
for (n in unique(met_anno_list)) {
  met_matrix_list[[paste("met",n,sep="_")]] <- met_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(intersect,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F, fun.aggregate=sum) %>% matrix.please() %>% t
  
  cat(sprintf("%s methylation matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(met_matrix_list[[paste("met",n,sep="_")]]), ncol(met_matrix_list[[paste("met",n,sep="_")]]),
              100*mean(is.na(met_matrix_list[[paste("met",n,sep="_")]]))))
}

cat("\n")

acc_matrix_list <- list()
for (n in unique(acc_anno_list)) {
  acc_matrix_list[[paste("acc",n,sep="_")]] <- acc_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(intersect,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F, fun.aggregate=sum) %>% matrix.please() %>% t
  
  cat(sprintf("%s accessibility matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(acc_matrix_list[[paste("acc",n,sep="_")]]), ncol(acc_matrix_list[[paste("acc",n,sep="_")]]),
              100*mean(is.na(acc_matrix_list[[paste("acc",n,sep="_")]]))))
}

all_matrix_list <- c(rna=list(rna_matrix),met_matrix_list,acc_matrix_list)

saveRDS(all_matrix_list, "data/all_matrix_list.rds")
'''
