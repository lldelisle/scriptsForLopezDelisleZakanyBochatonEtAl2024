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
safelyLoadAPackageInCRANorBioconductor("demuxmix")
safelyLoadAPackageInCRANorBioconductor("ggplot2")

# introduce variable that will be used
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
wd <- "/data/home/ldelisle/scriptsForLopezDelisleZakanyBochatonEtAl2024"
# all files are relative to wd
# rds files will be put in a RDS folder in the working directory
# rds are only created once, if pipeline changes, delete manually all rds
raw.files <- "../mountDuboule/Lucille/Bochaton/"
directoryGex <- file.path(raw.files, "GEX")
directoryCMO <- file.path(raw.files, "CMO")
RDSfolder <- "../mountDuboule/Lucille/Bochaton/RDS"
output.directory <- "scRNAseq"
output.metrics <- file.path(output.directory, "metrics/")

# metadata sample is a csv file which has to have
# a column named "Directory" which indicates the directory with cell ranger format of CMO and GEX.
# It also need to have a column "Sample.ID" which will have a unique number
# This allows you to define which experiment you wanna run
# It needs a Multiplexing column which value can be "CellPlex" or NA
# Finally it needs a RDS column which will give the names to the RDS files of the non CellPlex samples
# in summary 4 columns are required : "Sample.ID", "RDS", "Multiplexing" and "Directory"
metadata.file <- "scRNAseq/metadata.csv"
separator.in.metadata.file <- ","
# samples.used is defining the order for plotting
samples.used <- c(1:2)
min.cmo.umi <- 5

# Define library directory
# if needed, update on getwd()/.Renviron

#define functions
create.seurat.from.folder <- function(sample.name, directoryGex) {
  ## create a seurat object from a directory containing files in the cell ranger 10x format
  # and compute the percentage of mitochondria to later do filtering
  cat("Read GEX matrix...")
  all.sorted.raw.data <- Read10X(file.path(directoryGex, sample.name))
  cat("\n")
  seurat.object <- CreateSeuratObject(counts = all.sorted.raw.data, project = sample.name, min.cells = 3,
                                      min.features = 200)
  seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^mt-")
  return(seurat.object)
}

CellCycleScore.from.seuratObject <- function(seurat.object, output.metrics, s.genes, g2m.genes) {
  ## subset seurat object based on percent.mt, between 0.05% and 8% and nCount_RNA, between 0.4 and 2.5 of mean RNA
  # and compute cell cycle score
  # save metrics pre-subset
  dir.create(output.metrics, showWarnings = F, recursive = T)
  ggsave(filename = paste0(output.metrics, Project(seurat.object), ".RNA.QC.pre-subset.pdf"),
         VlnPlot(seurat.object,
                 features = c("nFeature_RNA",
                              "nCount_RNA",
                              "percent.mt"),
                 ncol = 3), device = "pdf")
  mean.rna <- mean(seurat.object$nCount_RNA)
  seurat.object.subset <- subset(seurat.object,
                                 subset = nCount_RNA > 0.4 * mean.rna &
                                   nCount_RNA < 2.5 * mean.rna &
                                   percent.mt < 8 & percent.mt > 0.05)
  ggsave(filename = paste0(output.metrics, Project(seurat.object.subset), ".RNA.QC.post-subset.pdf"),
         VlnPlot(seurat.object.subset,
                 features = c("nFeature_RNA",
                              "nCount_RNA",
                              "percent.mt"),
                 ncol = 3), device = "pdf")
  seurat.object.subset <- NormalizeData(seurat.object.subset, verbose = FALSE)
  seurat.object.subset <- FindVariableFeatures(seurat.object.subset,
                                               nfeatures = 3000,
                                               selection.method = "vst", verbose = FALSE)
  seurat.object.subset <- CellCycleScoring(seurat.object.subset,
                                           s.features = s.genes,
                                           g2m.features = g2m.genes,
                                           set.ident = TRUE, verbose = FALSE)
  return(seurat.object.subset)
}


# set working directory
setwd(wd)

#load metadata of scRNAseq samples that will be analyzed
mandatory.columns <- c("Sample.ID", "RDS", "Multiplexing", "Directory")
metadata.sample <- read.table(file = metadata.file, sep = separator.in.metadata.file,
                              header = TRUE, check.names = FALSE)
if (ncol(metadata.sample) < 4) {
  stop("Less than 4 columns detected into your metadata csv, check you correctly set the separator")
}
if (any(!mandatory.columns %in% colnames(metadata.sample))) {
  stop("Your metadata csv does not have the mandatory columns")
}
if (anyDuplicated(na.omit(metadata.sample$RDS))) {
  stop("Some rows have identical values in the RDS column, this is not possible")
}
# use names of each samples match is to match the number to the Sample.ID column of metadata file
all.directory <- unique(metadata.sample$Directory[match(samples.used, metadata.sample$Sample.ID)])

# this loop is done on each 10x well, for CellPlex samples,
# it annotates CellPlex cells to their respective sample (as defined in the CMO directory)
# and keep the "negative" annotation
# then, it subsets samples to remove low quality cells and putatively un-resolved doublets
# then, it computes the cellcyclescore of each cells
# and finally in the case of cellPlex split the object between samples and save their respective rds
# Here we set the list RDS variable which will give the names of RDS files to put in the RDS column of metadata
list.RDS <- ""
for (mydirectory in all.directory) {
  # Subset the metadata.sample
  all.meta <- metadata.sample[metadata.sample$Directory == mydirectory, ]
  # check if sample is CellPlex or regular 10x
  if (is.na(all.meta[1, "Multiplexing"])) {
    if (nrow(all.meta) != 1) {
      stop("directory Name present multiple time")
    }
    if (is.na(all.meta$RDS)) {
      stop("RDS value for a non CellPlex sample needs to be set")
    }
    if (!file.exists(file.path(RDSfolder, "single.sample", paste0(all.meta$RDS, ".RDS")))) {
      cat(paste0(mydirectory, " RDS is not found will proceed with this new sample\n"))
      cat("sample is not cellPlex\n")
      # Create Seurat object
      seurat.object <- create.seurat.from.folder(mydirectory, directoryGex)
      # Subset and compute cell cycle score
      seurat.object.subset <- CellCycleScore.from.seuratObject(seurat.object, output.metrics, s.genes, g2m.genes)
      # Save RDS
      dir.create(file.path(RDSfolder, "single.sample"), showWarnings = FALSE, recursive = TRUE)
      saveRDS(seurat.object.subset, file = file.path(RDSfolder, "single.sample", paste0(all.meta$RDS, ".RDS")))
      rds.file <- all.meta$RDS
      list.RDS <- paste0(list.RDS, "\n", rds.file)
    } else {
      cat(paste0(mydirectory, " RDS is found will skip this sample\n"))
    }
  } else if (all.meta[1, "Multiplexing"] == "CellPlex") {
    # if sample is 10x, perform annotation of cells to demultiplex per sample
    if (!file.exists(file.path(RDSfolder, "cellplex", paste0(mydirectory, ".RDS")))) {
      cat(paste0(mydirectory, " RDS is not found will proceed for new sample\n"))
      cat("sample is cellPlex, will perform demultiplexing\n")
      cat("Read CMO matrix...")
      cmo.umi.mtx <- Read10X(file.path(directoryCMO, mydirectory),
                             gene.column = 1)
      cat("\n done\n")
      # First remove the 'unmapped' row
      cmo.umi.mtx <- cmo.umi.mtx[setdiff(rownames(cmo.umi.mtx), "unmapped"), ]
      # This translation (see https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-)
      # Is part of the galaxy workflow
      # translation <- read.delim(file.path("C:/Users/mayran/Desktop/lab/R/2020.analyses/tests/",
      #                           "translation_3M-february-2018.txt.gz"),
      #                           header = F)
      # colnames(cmo.umi.mtx) <- translation$V1[match(colnames(cmo.umi.mtx), translation$V2)]
      # Filter cellular barcodes with less than min.cmo.umi umis assigned to sample barcodes
      cmo.umi.mtx.subset <- cmo.umi.mtx[, colSums(cmo.umi.mtx) > min.cmo.umi]
      # Create Seurat object
      seurat.object <- create.seurat.from.folder(mydirectory, directoryGex)
      # Filter the cmo matrix to only keep cellular barcodes which are in the Seurat object
      cmo.umi.mtx.subset.in.seurat <- cmo.umi.mtx.subset[, intersect(colnames(seurat.object),
                                                                     colnames(cmo.umi.mtx.subset))]
      # Demultiplex with demuxmix using the nCounts_RNA
      cat("performing demultiplexing...\n")
      dmm2 <- demuxmix(as.matrix(cmo.umi.mtx.subset.in.seurat),
                       rna = seurat.object[[]][colnames(cmo.umi.mtx.subset.in.seurat), "nCount_RNA"])
      cat("done\n")
      classification2 <- dmmClassify(dmm2)
      # Save the demuxmix result to tsv
      write.table(classification2, file = file.path(raw.files, "CMO", paste0(mydirectory, "classification2.tsv")))

      # add classifications on seurat object
      seurat.object$CMOtype <- classification2[colnames(seurat.object), "Type"]
      seurat.object$allCMO <- classification2[colnames(seurat.object), "HTO"]
      seurat.object$Prob.CMO <- classification2[colnames(seurat.object), "Prob"]
            # For cellular barcodes present in the Seurat object but not in the CMO demultiplexing result
      # We set the CMOtype to 'noCMO'
      seurat.object$CMOtype[is.na(seurat.object$CMOtype)] <- "noCMO"
      # sample.CMO will contain: noCMO, negative, multiplet, uncertain, or the name of the sample
      seurat.object$sample.CMO <- seurat.object$CMOtype
      seurat.object$sample.CMO[seurat.object$sample.CMO == "singlet"] <- seurat.object$allCMO[seurat.object$sample.CMO == "singlet"]
      # By default CITE-seq-count add the sequence of the sample barcode,
      # for example: WT_120h_rep1-CCGTCGTCCAAGCAT
      # So we remove the last '-' followed by the sequence
      seurat.object$sample.CMO <- gsub("-[ATCG]*$", "", seurat.object$sample.CMO)
      # display number of cells in each category pre-cleaning
      print(table(seurat.object$sample.CMO, exclude = NULL))
      print(table(seurat.object$CMOtype, exclude = NULL))
      # switches identity to keep only the singlet and cells with no sample assigned
      Idents(seurat.object) <- "CMOtype"
      seurat.object.clean <- subset(seurat.object, idents = c("singlet", "negative"))
      Idents(seurat.object.clean) <- "sample.CMO"
      print(table(seurat.object.clean$sample.CMO, exclude = NULL))
      print(table(seurat.object.clean$CMOtype, exclude = NULL))

      #subset the CellPlex object base on nCount and percent.mt and compute cell cycle scoring
      seurat.object.clean.subset <- CellCycleScore.from.seuratObject(seurat.object.clean, output.metrics, s.genes, g2m.genes)
      #save the RDS of the full Cellplex sample
      cat(paste0("saving RDS of ", mydirectory, "\n"))
      dir.create(file.path(RDSfolder, "cellplex"), showWarnings = FALSE, recursive = TRUE)
      saveRDS(seurat.object.clean.subset, file = file.path(RDSfolder, "cellplex", paste0(mydirectory, ".RDS")))
      #save the RDS of the individual cellplex sample
      split.seurat.by.sample <- SplitObject(seurat.object.clean.subset, split.by = "sample.CMO")
      if (nrow(all.meta) != length(split.seurat.by.sample) + 1) {
        warning("The number of row matching this directory does not match the number of samples in CMO")
      }
      dir.create(file.path(RDSfolder, "single.sample"), showWarnings = FALSE, recursive = TRUE)
      for (sample in names(split.seurat.by.sample)) {
        cat(paste0("saving RDS of ", sample, " from ", mydirectory, "\n"))
        if (!paste0(mydirectory, "_", sample) %in% all.meta$RDS) {
          warning(paste("The RDS file", paste0(mydirectory, "_", sample),
                        "will be generated but is not in the RDS column of metadata csv.",
                        "Fill the metadata csv before the next step."))
        }
        saveRDS(split.seurat.by.sample[[sample]], file.path(RDSfolder, "single.sample", paste0(mydirectory, "_", sample, ".RDS")))
        rds.file <- paste0(mydirectory, "_", sample)
        list.RDS <- paste0(list.RDS, "\n", rds.file)
      }
    } else {
      cat(paste0(mydirectory, " RDS is found will skip this sample\n"))
    }
  } else {
    stop("Category multiplexing not found")
  }
}
cat(paste0("\nHere is the list of rds files which needs to be accordingly present in your metadata file.\n",
           "Ensure to place CellPlex samples in the correct order as it is processed in a random order",
           list.RDS, "\n\n"))
