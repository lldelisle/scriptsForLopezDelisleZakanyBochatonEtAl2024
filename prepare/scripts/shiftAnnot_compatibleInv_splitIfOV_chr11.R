## This script takes as input a br file and an annotation file (bed, gtf...)
## It will create where the br file is, one folder per mutant genome with the shifted annotations inside.
## If an annotation overlap an insertion/deletion, it will be split in 2 annotations with the same names.
## Idem if it overlap an inversion.

###############################################
## The script assumes the chromosome is chr2 but you can change it here:
chrToShift <- "chr11"
###############################################

options(scipen = 999)
options(stringsAsFactors = F)
library(tools)

if (length(commandArgs(TRUE)) == 0) {
  script.basename <- getSrcDirectory(function(x) {
    x
  })
  cat("Choose the br file.\n")
  pathForBr <- file.choose()
  cat("Choose the file to convert.\n")
  fileToConvert <- file.choose()
  cat("Which is the number of the column with the chromosme name.\n")
  colChr <- as.numeric(readLines(con = stdin(), n = 1))
  cat("Which is the number of the column with the start position.\n")
  colStart <- as.numeric(readLines(con = stdin(), n = 1))
  cat("Which is the number of the column with the end position.\n")
  colEnd <- as.numeric(readLines(con = stdin(), n = 1))
  cat("Which is the number of the column with the strand information (put 0 if there is no).\n")
  colStrand <- as.numeric(readLines(con = stdin(), n = 1))
  cat("Is your annotation file BED-like (0-based half-open): put T or GTF-like (1-based closed): put F\n")
  isBED <- as.logical(readLines(con = stdin(), n = 1))
  cat("Which is the name you want to add to each output.\n")
  outputName <- readLines(con = stdin(), n = 1)
} else {
  if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help") {
    cat("Usage: Rscript shiftAnnot_compatibleInv_splitOV.R pathForBr fileToConvert colCr colStart colEnd colStrand BEDlike outputName \n")
    stop()
  }
  pathForBr <- commandArgs(TRUE)[1]
  fileToConvert <- commandArgs(TRUE)[2]
  colChr <- as.numeric(commandArgs(TRUE)[3])
  colStart <- as.numeric(commandArgs(TRUE)[4])
  colEnd <- as.numeric(commandArgs(TRUE)[5])
  colStrand <- as.numeric(commandArgs(TRUE)[6])
  isBED <- as.logical(commandArgs(TRUE)[7])
  outputName <- commandArgs(TRUE)[8]
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
}
################################################################################
other.name <- paste(sep = "/", script.basename, "shiftAnnotFunctions_compatibleInv.R")
source(other.name)
cat("Loading input file...")
if (file_ext(fileToConvert) == "gz") {
  firstChars <- substr(readLines(gzfile(fileToConvert), n = 1), 1, 5)
  if (firstChars == "track") {
    annotationData <- read.table(gzfile(fileToConvert), sep = "\t", stringsAsFactors = F, comment.char = "#", quote = "", skip = 1)
  } else {
    annotationData <- read.table(gzfile(fileToConvert), sep = "\t", stringsAsFactors = F, comment.char = "#", quote = "")
  }
} else {
  firstChars <- substr(readLines(fileToConvert, n = 1), 1, 5)
  if (firstChars == "track") {
    annotationData <- read.table(fileToConvert, sep = "\t", stringsAsFactors = F, comment.char = "#", quote = "", skip = 1)
  } else {
    annotationData <- read.table(fileToConvert, sep = "\t", stringsAsFactors = F, comment.char = "#", quote = "")
  }
}
cat("loaded.\n")
outputFolder <- dirname(pathForBr)
# Change to UCSC:
if (!(grepl("chr", annotationData[1, colChr]))) {
  annotationData[, colChr] <- paste0(rep("chr", nrow(annotationData)), annotationData[, colChr])
  annotationData[annotationData[, colChr] == "chrMT", colChr] <- "chrM"
  todelete <- grep("^chr(GL|JH)", annotationData[, colChr])
  annotationData <- annotationData[-todelete, ]
}
# Change to 1-based if needed:
if (isBED) {
  annotationData[, colStart] <- annotationData[, colStart] + 1
}
brdf <- read.delim(pathForBr, stringsAsFactors = F)
brdf[brdf == ""] <- NA
for (i in 1:nrow(brdf)) {
  cat("#####################\n")
  print(brdf[i, ])
  cat("####################\n")
  genome <- brdf$genome[i]
  shiftedAnnotations <- shiftDFFromBR(annotationData, genome, brdf, colChr, colStart, colEnd, colStrand, verbose = T, chromoWithTg = chrToShift, splitIfOverlap = T)
  # Put back to BED if needed:
  if (isBED) {
    shiftedAnnotations[, colStart] <- shiftedAnnotations[, colStart] - 1
  }
  cat("\n\n\n\n")
  dir.create(paste0(outputFolder, "/", genome), showWarnings = F)
  write.table(shiftedAnnotations, paste0(outputFolder, "/", genome, "/", genome, "_vSplit_", outputName), sep = "\t", row.names = F, col.names = F, quote = F)
}
