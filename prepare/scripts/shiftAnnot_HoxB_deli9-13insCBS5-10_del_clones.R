## This script takes as input an annotation file (bed, gtf...)
## A cassette description file with the results of blat
## A BR file
## I assume that there is no inversion and that the cassette is at least 100bp
## It will shift the annotations.
## If an annotation overlap an insertion/deletion, it will be split in 2 annotations with the same names.
## Idem if it overlap an inversion.

options(scipen = 999)
options(stringsAsFactors = F)
# if (!"devtools" %in% installed.packages()){
#   install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
# }
# devtools::install_github("lldelisle/usefulLDfunctions")
# library(usefulLDfunctions)
library(tools)

if (length(commandArgs(TRUE)) == 0) {
  script.basename <- getSrcDirectory(function(x) {
    x
  })
  cat("Choose the file to convert.\n")
  fileToConvert <- file.choose()
  cat("Choose the cassette description file.\n")
  fileCassette <- file.choose()
  cat("Choose the br file.\n")
  fileBR <- file.choose()
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
  cat(paste0("Full output path(Input is ", fileToConvert, ")\n"))
  outputName <- readLines(con = stdin(), n = 1)
} else {
  if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help") {
    cat("Usage: Rscript shiftAnnot_HoxB_deli9-13insCBS5-10_del_clones.R fileToConvert CassetteDesc Br colChr colStart colEnd colStrand isBED fullOutputPath \n")
    stop()
  }
  fileToConvert <- commandArgs(TRUE)[1]
  fileCassette <- commandArgs(TRUE)[2]
  fileBR <- commandArgs(TRUE)[3]
  colChr <- as.numeric(commandArgs(TRUE)[4])
  colStart <- as.numeric(commandArgs(TRUE)[5])
  colEnd <- as.numeric(commandArgs(TRUE)[6])
  colStrand <- as.numeric(commandArgs(TRUE)[7])
  isBED <- as.logical(commandArgs(TRUE)[8])
  outputName <- commandArgs(TRUE)[9]
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
}
################################################################################
# Source the functions which are in the other file:
other.name <- file.path(script.basename, "shiftAnnotFunctions_compatibleInv.R")
source(other.name)

cat("Loading cassette file..")
cassette.df <- read.delim(fileCassette)
if (!all(colnames(cassette.df) == c("start_c", "end_c", "start_chr11", "end_chr11"))) {
  stop("The cassette file is not valid")
}

cat("Loading BR file..")
brdf_deli913insCTCF_clone <- read.delim(fileBR)
if (!all(colnames(brdf_deli913insCTCF_clone)[1:4] == c("genome", "br1", "br2", "tg1"))) {
  stop("The BR file is not valid")
}
if (nrow(brdf_deli913insCTCF_clone) != 1) {
  stop("Only BR with only one row are accepted")
}

cat("Loading input file...")
# annotationData <- usefulLDfunctions:::.readFileFromConditionOnNcols(fileToConvert,
#                                                                     paste0(">=", max(colChr, colStart, colEnd, colStrand)),
#                                                                     keepQuote=T)
annotationData <- .readFileFromConditionOnNcols(fileToConvert,
  paste0(">=", max(colChr, colStart, colEnd, colStrand)),
  keepQuote = T
)
cat("loaded.\n")

# Convert to UCSC format
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
# I assume that there is no inversion and that the cassette is at least 100bp
i <- 1
shift <- 0
while (brdf_deli913insCTCF_clone[1, 1 + i * 3] < 100) {
  shift <- shift - (brdf_deli913insCTCF_clone[1, 1 + i * 3 - 1] - brdf_deli913insCTCF_clone[1, 1 + i * 3 - 2] + 1)
  i <- i + 1
}
cassette.begin <- brdf_deli913insCTCF_clone[1, 1 + i * 3 - 2] + shift
# The coordinate of CTCF cassette (1) corresponds to
# cassette.begin in the mutant genome coordinate
annot.CTCFcassette.list <- apply(
  cassette.df,
  1,
  function(v) {
    temp.df <- getDFfromCooPotInv(annotationData,
      colChr, colStart,
      colEnd, colStrand,
      chrToGet = "chr11",
      startToGet = v["start_chr11"],
      endToGet = v["end_chr11"]
    )
    if (nrow(temp.df) == 0) {
      return(NULL)
    }
    return(shiftDFOfBP(temp.df,
      colStart,
      colEnd,
      bpToShift = cassette.begin - 2 + v["start_c"]
    ))
  }
)


if (all(sapply(annot.CTCFcassette.list, is.null))) {
  annot.CTCFcassette <- NULL
} else {
  annot.CTCFcassette <- do.call(rbind, annot.CTCFcassette.list)
}

# Shift with br files
brdf_delB <- read.delim(text = "genome\tbr1\tbr2\ttg1\ndelB\t96193822\t96368393\t0")

cat("Shifting in chr11:\n")
shiftedAnnotations_deli913insCTCF_clone <-
  shiftDFFromBR(annotationData,
    genome = brdf_deli913insCTCF_clone$genome[1], brdf_deli913insCTCF_clone, colChr, colStart, colEnd,
    colStrand, verbose = T, chromoWithTg = "chr11", splitIfOverlap = T
  )

# Construct chr11_delB
annotations_chr11 <- annotationData[annotationData[, colChr] == "chr11", ]

cat("Shifting in chr11_delB:\n")
shiftedAnnotations_chr11 <-
  shiftDFFromBR(annotations_chr11,
    genome = "delB", brdf_delB, colChr, colStart, colEnd,
    colStrand, verbose = T, chromoWithTg = "chr11", splitIfOverlap = T
  )
if (nrow(shiftedAnnotations_chr11) > 1) {
  shiftedAnnotations_chr11[, colChr] <- "chr11_delB"
}

shiftedAnnotations <- rbind(
  shiftedAnnotations_deli913insCTCF_clone,
  shiftedAnnotations_chr11,
  annot.CTCFcassette
)
# Put back to BED if needed:
if (isBED) {
  shiftedAnnotations[, colStart] <- shiftedAnnotations[, colStart] - 1
}

dir.create(dirname(outputName), showWarnings = F)
write.table(shiftedAnnotations,
  file = outputName,
  sep = "\t", row.names = F, col.names = F, quote = F
)
