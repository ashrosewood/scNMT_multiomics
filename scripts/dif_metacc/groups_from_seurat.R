args <- commandArgs()

help <- function(){
    cat("groups_from_seurat.R :
Grabs groups from Seurat object\n")
    cat("Usage: \n")
    cat("--SO        : Seurat object from transcriptome analysis   [required]\n")
    cat("--out       : filename for output                         [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$seurat   <- sub( '--SO=', '', args[grep('--SO', args)] )
    io$out      <- sub( '--out=', '', args[grep('--out', args)] )
}


library(Seurat)
library(data.table)
library(purrr)

fwrite_tsv <- partial(fwrite, sep = "\t", na="NA")

### TEST INPUT ###
#io <- list()
#io$seurat <- "../scRNA_SMARTseq2/data/seurat/SeuratObject.rds"
#io$out <- "data/groups.tsv"

dir.create(dirname(io$out))

seurat <- readRDS(io$seurat)
Idents(seurat) <- "origin"
groups <- Idents(seurat) %>% 
  as.data.frame() %>% 
  setDT(keep.rownames = "id_rna") %>%
  setnames(".", "group", skip_absent=T)  


groups[, id := strsplit(id_rna, "_") %>% 
         map_chr(2) %>% 
         gsub("([A-H])0", "\\1", .) %>% 
         paste0("sc_", .)]

groups$group <- substr(groups$group, 1, 3)

fwrite_tsv(groups, io$out)
