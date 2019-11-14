args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--name        : name for output model                  [required]\n")
    cat("--outdir      : output directory                       [required]\n")
    cat("--data        : data matrix list                       [required]\n")
    cat("--scale       : TRUE or FALSE                          [required]\n")
    cat("--factor      : number of latent factors in model      [default = 10]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$data     <- sub( '--data=', '', args[grep('--data', args)] )
    io$trial    <- sub( '--name=', '', args[grep('--name', args)] )
    io$outdir   <- sub( '--outdir=', '', args[grep('--outdir', args)] )
    opts$scale  <- sub( '--scale=', '', args[grep('--scale', args)] )
    num_factor  <- sub( '--factor=', '', args[grep('--factor', args)] )
}

if (identical(num_factor,character(0))){
   opts$factor <- 10
}else{
   opts$factor <- as.numeric(num_factor)
}

suppressMessages(library(MOFA))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(scater))
suppressMessages(library(reticulate))
suppressMessages(library(argparse))


#io$outdir <- "data"
#io$trial <- "test"
#io$data <- "data/all_matrix_list.rds"
#opts$scale <- FALSE
#opts$factor <- 10

# print(py_config())
# print(.libPaths())

my_data <- readRDS(io$data)

my_data[[3]] <- my_data[[3]][-996,]

# Create MOFAobject
MOFAobject <- createMOFAobject(my_data)

# Data processing options
DataOptions <- getDefaultDataOptions()
DataOptions$scaleViews <- opts$scale

# Model options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- opts$factor

# Training options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 5000
TrainOptions$DropFactorThreshold <- 0.01

# Prepare MOFAobject for training
MOFAmodel <- prepareMOFA(MOFAobject, 
  DataOptions = DataOptions, 
  ModelOptions = ModelOptions, 
  TrainOptions = TrainOptions
)

# Train the model
outfile <- sprintf("%s/hdf5/model_%s.hdf5",io$outdir,io$trial)
model <- runMOFA(MOFAmodel, outfile)










