# Install required packages
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
if (!"usefulLDfunctions" %in% installed.packages()) {
  devtools::install_github("lldelisle/usefulLDfunctions")
}
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("digest")

# introduce variable that will be used
wd <- "/data/home/ldelisle/scriptsForLopezDelisleZakanyBochatonEtAl2024"
# all files are relative to wd
# rds files will be put in a RDS folder
# rds are only created once, if pipeline changes, delete manually all rds
RDSfolder <- "../mountDuboule/Lucille/Bochaton/RDS"
# metadata sample is a csv file which has to have a column named
# "Sample.ID" which will have a unique number allowing you to define which experiment you wanna merge
# It also needs a "RDS" column which will be used to load RDS files. Do not add the extension .RDS in this column
# finally a column "Directory" needs to be setup
# if 2 CellPlex or regular 10x experiments were done label them 1, 2, ...
# in summary the columns needed are: Sample.ID, RDS, Directory
# all column will be added as metadata to seurat object so you can make any additional column
metadata.file <- "scRNAseq/metadata.csv"
separator.in.metadata.file <- ","
# samples.used is defining the order for plotting

samples.used <- 1:2
desiredName.for.RDS <- "combined.WT.mutant"

# if you want to remove some genes from analysis, please add them in this variable
manual.removal.genes <- NULL

# to remove Hox genes from PCA calculation (and thus from UMAP and clustering) uncomment this
# hoxa.b.d.genes <-paste0("Hox", rep(letters[c(1,2,4)], each = 13), rep(1:13, 3))
# manual.removal.genes <- hoxa.b.d.genes

seurat.all.sorted.subset <- list()

# set working directory
setwd(wd)
current.version <- "2023.11.06"

# load metadata of scRNAseq samples that will be analyzed
metadata.sample <- read.table(file = metadata.file, sep = separator.in.metadata.file,
                              header = TRUE, check.names = FALSE)
all.sorted <- metadata.sample$RDS[match(samples.used, metadata.sample$Sample.ID)]
if (any(is.na(all.sorted)) || any(all.sorted == "")) {
  stop("One of your selected sample has no value for the column RDS.")
}
hash <- digest(object = paste(sort(c(all.sorted, current.version)), collapse = ","), algo = "crc32", serialize = FALSE)
nameRDS <- paste0(desiredName.for.RDS, "_", hash, ".RDS")
if (!file.exists(file.path("RDS", "merged", nameRDS))) {
  for (my.sample in all.sorted) {
    cat(paste0("loading ",my.sample, "\n"))
    if (!file.exists(file.path(RDSfolder, "single.sample", paste0(my.sample, ".RDS")))) {
      stop("The RDS file does not exists, run step1 or fix the metadata csv")
    }
    seurat.all.sorted.subset[[my.sample]] <- readRDS(file.path(RDSfolder, "single.sample", paste0(my.sample, ".RDS")))
    all.meta <- metadata.sample[metadata.sample$RDS == my.sample, ]
    #add metadata from the metadata.sample for each samples
    for (current.metadata in colnames(all.meta)) {
      seurat.all.sorted.subset[[my.sample]][[current.metadata]] <- all.meta[[current.metadata]]
    }
  }
  cat("Merging dataset...\n")
  combined.seurat <- merge(seurat.all.sorted.subset[[all.sorted[1]]],
                           y = subsetByNamesOrIndices(seurat.all.sorted.subset,
                                                      all.sorted[2:length(all.sorted)]),
                           add.cell.ids = all.sorted,
                           project = "aggregate.all")
  cat("done \n")
  #set order for ploting seurat object in the same order as the sample.used
  
  for (current.metadata in colnames(metadata.sample)) {
    combined.seurat[[current.metadata]] <- factor(combined.seurat[[]][[current.metadata]],
                                                  levels = unique(metadata.sample[match(samples.used,
                                                                                        metadata.sample$Sample.ID),
                                                                                  current.metadata]))
  }
  # Normalize Data, find variable features, scale and regress cell cycle score
  combined.seurat <- NormalizeData(combined.seurat,
                                   verbose = FALSE)
  combined.seurat <- FindVariableFeatures(combined.seurat,
                                          verbose = FALSE, nfeatures = 2000)
  cat("Scaling merged object..\n")
  combined.seurat <- ScaleData(combined.seurat,
                               verbose = TRUE,
                               vars.to.regress = c("percent.mt", "S.Score",
                                                   "G2M.Score"))
  cat("Done\n")
  # checking number of 10x well to see whether batch correction is useful
  if (length(unique(metadata.sample$Directory[match(samples.used, metadata.sample$Sample.ID)])) > 1) {
    cat("more than one 10x well detected, restricting variable features to limit batch effect",
        "\nvariable features will include gene with expression between 0.05 and 0.8 percentile\n")
    average.expression <- rowMeans(AverageExpression(combined.seurat, group.by = "Directory")[[1]])
    genes.to.exclude <- names(average.expression[average.expression < quantile(average.expression,
                                                                               0.05) | average.expression > quantile(average.expression,
                                                                                                                     0.8)])
    
    
    used.var.gene <- setdiff(VariableFeatures(combined.seurat),
                             c(genes.to.exclude, manual.removal.genes))
    
  } else {
    cat("Only experiments from a single 10x well is detected,\n 
      all variable genes are used as there shouldn't be batch effect\n")
    #run PCA analysis
    used.var.gene <- setdiff(VariableFeatures(combined.seurat),
                             c(manual.removal.genes))
  }
  cat("Performing PCA analysis..\n")
  combined.seurat <- RunPCA(combined.seurat, npcs = 50, verbose = FALSE,
                            features = used.var.gene)
  cat("Done\n")
  cat("Saving RDS file...\n")
  dir.create(file.path(RDSfolder,"merged"), showWarnings = F, recursive = T)
  saveRDS(combined.seurat, file.path(RDSfolder,"merged", nameRDS))
  cat("Done\n")
} else {
  cat("RDS of merged datasets already exist, proceed directly with step3\n")
}
cat(paste0("the seurat object has been saved to: ", nameRDS,"\n"))

