options(stringsAsFactors = F)
# This function invert the strand (+ -> - and - -> +, anything else stay the same)
changeStrand <- function(s) {
  if (s == "+") {
    return("-")
  } else if (s == "-") {
    return("+")
  } else {
    return(s)
  }
}

# This function will shift the annotationData according to the
# br1=position of start deletion/inversion,
# br2=position of end deletion/inversion and
# length.transgene which is either the size of the transgene
# or the inv keyword,
# to be compatible with any annotationData you need to specify
# which column is for chromosome, start, end
# and if you have an inversion, strand.
# You can choose to split the annotation
# if it totally overlap the 2 breakpoints
# in the case of a transgene insertion.
# We assume here that all coordinates are 1-base included
shiftDF <- function(annotationData, br1, br2, length.transgene,
                    colChr, colStart, colEnd, colStrand = 0, verbose = T,
                    chromoWithTg = "chr2", splitIfOverlap = F) {
  annotationDF <- annotationData
  # shiftAnnot is the new annotation which will be returned.
  shiftAnnot <- annotationDF
  length.deletion <- br2 - br1 + 1
  # First, look for the annotations between the two breakpoints.
  deletedOrInverted <- which(annotationDF[, colChr] == chromoWithTg &
    annotationDF[, colStart] >= br1 &
    annotationDF[, colEnd] <= br2)

  if (length(deletedOrInverted) > 0) {
    nbdel <- length(deletedOrInverted)
    if (verbose) {
      cat(paste("deleting/inverting", nbdel, "annotations.\n"))
      cat(paste("First ones:\n"))
      print(head(annotationDF[deletedOrInverted, ]))
      cat("\n")
    }
    if (length.transgene == "inv") {
      # If it is an inversion, you need to change start, end and strand
      shiftAnnot[deletedOrInverted, colStart] <- br2 - (annotationDF[deletedOrInverted, colEnd] - br1)
      shiftAnnot[deletedOrInverted, colEnd] <- br2 - (annotationDF[deletedOrInverted, colStart] - br1)
      if (colStrand != 0) {
        shiftAnnot[deletedOrInverted, colStrand] <- sapply(annotationDF[deletedOrInverted, colStrand], changeStrand)
      }
    } else {
      # Else you just need to remove them
      shiftAnnot <- shiftAnnot[-deletedOrInverted, ]
      annotationDF <- annotationDF[-deletedOrInverted, ]
    }
  }

  # Second, look for the annotations which begin before the first breakpoint and end after the second one.
  ov.all <- which(annotationDF[, colChr] == chromoWithTg & annotationDF[, colStart] < br1 & annotationDF[, colEnd] > br2)

  if (length(ov.all) > 0) {
    nbov <- length(ov.all)
    if (verbose) {
      cat(paste(nbov, "annotations cross the whole boundary.\n"))
      cat(paste("First ones:\n"))
      print(head(annotationDF[ov.all, ]))
      cat("\n")
    }
    if (length.transgene != "inv") {
      # If it is an insertion/deletion
      if (splitIfOverlap) {
        # If you want to split the annotation
        # The initial annotation is shorten on the right
        shiftAnnot[ov.all, colEnd] <- br1 - 1
        # Another annotation is copied and start and end are changed to match the new coordinates.
        newAnnot <- annotationData[ov.all, ]
        newAnnot[, colStart] <- br2 + 1 - length.deletion + as.numeric(length.transgene)
        newAnnot[, colEnd] <- newAnnot[, colEnd] - length.deletion + as.numeric(length.transgene)
        shiftAnnot <- rbind(shiftAnnot, newAnnot)
      } else {
        # If you do not want to split, you just change the right to match the insertion/deletion.
        shiftAnnot[ov.all, colEnd] <- shiftAnnot[ov.all, colEnd] - length.deletion + as.numeric(length.transgene)
      }
    } else {
      # If it is an inversion, only if you have a strand information you change the initial annotation.
      if (colStrand != 0) {
        # Then the annotation will be present 3 times
        # First only shorten on the right
        shiftAnnot[ov.all, colEnd] <- br1 - 1
        # In the middle it is inverted and correspond to the inversion coordinates
        newAnnotMiddle <- annotationData[ov.all, ]
        newAnnotMiddle[, colStart] <- br1
        newAnnotMiddle[, colEnd] <- br2
        newAnnotMiddle[, colStrand] <- sapply(newAnnotMiddle[, colStrand], changeStrand)
        shiftAnnot <- rbind(shiftAnnot, newAnnotMiddle)
        # Finally on the right it is shorten on the left
        newAnnotRight <- annotationData[ov.all, ]
        newAnnotRight[, colStart] <- br2 + 1
        shiftAnnot <- rbind(shiftAnnot, newAnnotRight)
      }
    }
  }

  # Third look at annotations which overlap the breakpoint 2.
  ov.right <- which(annotationDF[, colChr] == chromoWithTg & annotationDF[, colStart] <= br2 & annotationDF[, colEnd] > br2 & annotationDF[, colStart] >= br1)

  if (length(ov.right) > 0) {
    nbov <- length(ov.right)
    if (verbose) {
      cat(paste(nbov, " annotations overlap the right-hand deletion boundary.\n"))
      cat(paste("First ones:\n"))
      print(head(annotationDF[ov.right, ]))
      cat("\n")
    }
    # If it is an inversion
    if (length.transgene == "inv") {
      # I split into 2
      newAnnot <- annotationDF[ov.right, ]
      shiftAnnot[ov.right, colStart] <- br2 + 1
      newAnnot[, colStart] <- br1
      newAnnot[, colEnd] <- br1 + (br2 - annotationDF[ov.right, colStart])
      if (colStrand != 0) {
        newAnnot[, colStrand] <- sapply(newAnnot[, colStrand], changeStrand)
      }
      shiftAnnot <- rbind(shiftAnnot, newAnnot)
    } else {
      # If it is an insertion/deletion I change the left coordinate to match the end of the insertion
      shiftAnnot[ov.right, colStart] <- br2 + 1 - length.deletion + as.numeric(length.transgene)
      # I change the right coordinate to make it compatible with the insertion/deletion size.
      shiftAnnot[ov.right, colEnd] <- shiftAnnot[ov.right, colEnd] - length.deletion + as.numeric(length.transgene)
    }
  }

  # Roughly the same for the annotations overlapping the first breakpoint.
  ov.left <- which(annotationDF[, colChr] == chromoWithTg & annotationDF[, colStart] < br1 & annotationDF[, colEnd] >= br1 & annotationDF[, colEnd] <= br2)
  if (length(ov.left) > 0) {
    nbov <- length(ov.left)
    if (verbose) {
      cat(paste(nbov, " annotations overlap the left-hand deletion boundary.\n"))
      cat(paste("First ones:\n"))
      print(head(annotationDF[ov.left, ]))
      cat("\n")
    }
    if (length.transgene == "inv") {
      # I split into 2
      newAnnot <- annotationDF[ov.left, ]
      shiftAnnot[ov.left, colEnd] <- br1 - 1
      newAnnot[, colEnd] <- br2
      newAnnot[, colStart] <- br2 - (annotationDF[ov.left, colEnd] - br1)
      if (colStrand != 0) {
        newAnnot[, colStrand] <- sapply(newAnnot[, colStrand], changeStrand)
      }
      shiftAnnot <- rbind(shiftAnnot, newAnnot)
    } else {
      shiftAnnot[ov.left, colEnd] <- br1 - 1
    }
  }

  # Finally if it is not an inversion, all downstream annotation need to change the coordiates.

  if (length.transgene != "inv") {
    downstream <- which(annotationDF[, colChr] == chromoWithTg & annotationDF[, colStart] > br2)

    if (length(downstream) > 0) {
      nbdown <- length(downstream)
      if (verbose) {
        cat(paste(nbdown, " after the deletedOrInverted region.\n"))
      }
      shiftAnnot[downstream, colStart] <- annotationDF[downstream, colStart] - length.deletion + as.numeric(length.transgene)
      shiftAnnot[downstream, colEnd] <- annotationDF[downstream, colEnd] - length.deletion + as.numeric(length.transgene)
    }
  }
  return(shiftAnnot)
}

# This function will shift sequencially the annotations using the information in the br dataframe for the given genome.
shiftDFFromBR <- function(annotationData, genome, br, colChr, colStart, colEnd, colStrand = 0, verbose = T, chromoWithTg = "chr2", splitIfOverlap = F) {
  if (genome == "mm10" | tolower(genome) == "wt") {
    return(annotationData)
  } else if (genome %in% br$genome) {
    brs <- br[br$genome == genome, ]
    shiftedAnnot <- annotationData
    i <- 1
    # While the start position of the breakpoint is defined in the table:
    while ((1 + 3 * (i - 1) + 1) < ncol(brs)) {
      if (is.na(brs[1, 1 + 3 * (i - 1) + 1])) {
        # If it is NA that means that there is no more step.
        i <- Inf
      } else {
        # Else we shift the annotations:
        shiftedAnnot <- shiftDF(
          shiftedAnnot, brs[1, 1 + 3 * (i - 1) + 1], brs[1, 1 + 3 * (i - 1) + 2], brs[1, 1 + 3 * (i - 1) + 3],
          colChr, colStart, colEnd, colStrand, verbose, chromoWithTg, splitIfOverlap
        )
        # If it was an insertion/deletion, the next breakpoints need to be updated:
        if (brs[1, 1 + 3 * (i - 1) + 3] != "inv") {
          brs[1, c(seq(1 + 3 * (i - 1) + 1, ncol(brs), 3), seq(1 + 3 * (i - 1) + 2, ncol(brs), 3))] <- # Starts and ends
            brs[1, c(seq(1 + 3 * (i - 1) + 1, ncol(brs), 3), seq(1 + 3 * (i - 1) + 2, ncol(brs), 3))] - # Initial Starts and ends -
            (brs[1, 1 + 3 * (i - 1) + 2] - brs[1, 1 + 3 * (i - 1) + 1] + 1) + as.numeric(brs[1, 1 + 3 * (i - 1) + 3]) # Current end - current start + 1 - current length of the transgene
        }
        i <- i + 1
      }
    }
    return(shiftedAnnot)
  }
}

# This function will get in the annotationData
# The annotations on chrToGet that overlaps the interval
# startToGet-endToGet
# All annotations will be shifted so that the startToGet will be 1
getDFfromCoo <- function(annotationData, colChr, colStart, colEnd,
                         chrToGet, startToGet, endToGet) {
  # getAnnot is the new annotation which will be returned.
  getAnnot <- annotationData[annotationData[, colChr] == chrToGet &
    annotationData[, colStart] <= endToGet &
    annotationData[, colEnd] >= startToGet, ]
  if (nrow(getAnnot) == 0) {
    return(getAnnot)
  }
  # I cut around the coordinates to get
  getAnnot[getAnnot[, colStart] <= startToGet, colStart] <- startToGet
  getAnnot[getAnnot[, colEnd] >= endToGet, colEnd] <- endToGet
  # I shift to get startToGet at 1
  getAnnot[, colStart] <- getAnnot[, colStart] - (startToGet - 1)
  getAnnot[, colEnd] <- getAnnot[, colEnd] - (startToGet - 1)

  return(getAnnot)
}

shiftDFOfBP <- function(annotationData, colStart, colEnd, bpToShift) {
  shiftedAnnot <- annotationData
  shiftedAnnot[, colStart] <- shiftedAnnot[, colStart] + bpToShift
  shiftedAnnot[, colEnd] <- shiftedAnnot[, colEnd] + bpToShift
  return(shiftedAnnot)
}

# From https://github.com/lldelisle/usefulLDfunctions/blob/0fcac5da1fd783e409c958644443a33905685262/R/myBasicFunctions.R#LL81C1-L129C1


#' Put the content of a tab separated file (with or without header gzip or not) in a dataframe
#' From the first line where the number of fields follow the condition cond
#'
#' @param fn the name of the file (tab delimited file with optionnally headers). Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param cond the condition that the number of columns should follow (for example "==4" or ">=3")
#' @param keepQuote logical whether the quotes should be kept useful for gtf files (default is F)
#' @return The dataframe containing the values of the file \code{fn} (the header is removed).
#' @importFrom utils read.delim
.readFileFromConditionOnNcols <- function(fn, cond, keepQuote = F) {
  # Require packages base, utils
  # i will be the first line (excluding commented lines) with data (no header)
  i <- 1
  while (TRUE) {
    # header is a data.frame containing the i-th line
    # (excluding commented lines)
    header <- tryCatch(
      utils::read.delim(gzfile(fn),
        nrows = 1,
        h = F, skip = (i - 1),
        comment.char = "#"
      ),
      error = function(e) {
        NULL
      }
    )
    if (is.null(header)) {
      return(NULL)
    }
    if (!eval(parse(text = paste0("ncol(header)", cond)))) {
      # if the number of columns does not fill the condition cond
      # the i-th line (excluding comments) is a header
      i <- i + 1
    } else {
      # if the number of columns fills the condition cond
      # the i-th line (excluding comments) is the first one with data
      break
    }
  }
  # change the quote char if keepQuote
  if (keepQuote) {
    quoteChar <- ""
  } else {
    # I use the default of read.delim
    quoteChar <- "\""
  }
  # return the data frame from the i-th line (excluding comments)
  return(utils::read.delim(gzfile(fn),
    h = F,
    skip = (i - 1),
    comment.char = "#",
    quote = quoteChar
  ))
}

# This function will get in the annotationData
# The annotations on chrToGet that overlaps the interval
# startToGet-endToGet
# If startToGet <= endToGet
# All annotations will be shifted so that the startToGet will be 1
# If startToGet > endToGet
# The orientation of annotations will be changed and
# The annotations will be shifted so that the endToGet will be 1
getDFfromCooPotInv <- function(annotationData, colChr, colStart, colEnd, colStrand,
                               chrToGet, startToGet, endToGet) {
  if (startToGet <= endToGet) {
    return(getDFfromCoo(annotationData, colChr, colStart, colEnd, chrToGet, startToGet, endToGet))
  }
  # We assume startToGet > endToGet
  # First we get the annotations in the same orientation:
  temp.annots <- getDFfromCoo(annotationData, colChr, colStart, colEnd, chrToGet, endToGet, startToGet)
  getAnnot <- temp.annots
  if (nrow(getAnnot) == 0) {
    return(getAnnot)
  }
  # We change strand:
  getAnnot[, colStrand] <- sapply(temp.annots[, colStrand], changeStrand)
  # We change the start:
  getAnnot[, colStart] <- (startToGet - endToGet  + 1) - temp.annots[, colEnd] + 1
  # We change the end:
  getAnnot[, colEnd] <- (startToGet - endToGet + 1) - temp.annots[, colStart] + 1
  return(getAnnot)
}
