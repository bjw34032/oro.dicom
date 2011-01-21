##
## Copyright (c) 2010, Brandon Whitcher
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

getOrientation <- function(xyz, delta=0.0001) {
  oX <- ifelse(xyz[1] < 0, "R", "L")
  oY <- ifelse(xyz[2] < 0, "A", "P")
  oZ <- ifelse(xyz[3] < 0, "F", "H")
  aX <- abs(xyz[1])
  aY <- abs(xyz[2])
  aZ <- abs(xyz[3])

  orientation <- NULL
  for (i in 1:3) {
    if (aX > delta && aX > aY && aX > aZ) {
      orientation <- paste(orientation, oX, sep="")
      aX <- 0
    } else {
      if (aY > delta && aY > aX && aY > aZ) {
        orientation <- paste(orientation, oY, sep="")
        aY <- 0
      } else {
        if (aZ > delta && aZ > aX && aZ > aY) {
          orientation <- paste(orientation, oZ, sep="")
          aZ <- 0
        }
      }
    }
  }
  return(orientation)
}

swapDimension <- function(img, dcm, digits=4) {
  imagePositionPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImagePositionPatient", FALSE), 3)
  if (nrow(imagePositionPatient) != nsli(img)) {
    imagePositionPatient <- attributes(img)$ipp
  }
  imageOrientationPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImageOrientationPatient", FALSE), 6)
  ## Ensure all rows of imageOrientationPatient are identical!
  pixelSpacing <-
    header2matrix(extractHeader(dcm$hdr, "PixelSpacing", FALSE), 2)
  ## Ensure all rows of pixelSpacing are identical!
  sliceThickness <- extractHeader(dcm$hdr, "SliceThickness")
  pixdim <- c(unique(pixelSpacing), unique(sliceThickness))
  iop.signif <- signif(imageOrientationPatient, digits)
  first.row <- getOrientation(unique(iop.signif)[1:3])
  first.col <- getOrientation(unique(iop.signif)[4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient (hope for the best).")
  }
  X <- nrow(img)
  Y <- ncol(img)
  Z <- dim(img)[3]
  W <- dim(img)[4]
  ld <- as.numeric(length(dim(img)))
  ## AXIAL
  if (is.axial(imageOrientationPatient)) {
    if (unlist(strsplit(first.row,""))[1] %in% c("A","P")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row,""))[1] == "R") {
      img <- switch(as.character(ld), "3" = img[X:1,,], "4" = img[X:1,,,],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col,""))[1] == "A") {
      img <- switch(as.character(ld), "3" = img[,Y:1,], "4" = img[,Y:1,,])
    }
    ## The z-axis is increasing toward the HEAD of the patient.
    z.index <- order(imagePositionPatient[,3])
    if (ld == 3) {
      img <- img[,,z.index]
    }
    if (ld == 4) {
      img <- img[,,Z:1,]
    }
  }
  ## CORONAL
  if (is.coronal(imageOrientationPatient)) {
    if (unlist(strsplit(first.row,""))[1] %in% c("H","F")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row,""))[1] == "R") {
      img <- switch(as.character(ld), "3" = img[X:1,,], "4" = img[X:1,,,],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col,""))[1] == "H") {
      img <- switch(as.character(ld), "3" = img[,Y:1,], "4" = img[,Y:1,,])
    }
    ## The y-axis is increasing to the posterior side of the patient.
    z.index <- order(imagePositionPatient[,2])
    if (ld == 3) {
      img <- img[,,z.index]
      index <- c(1,3,2)
    }
    if (ld == 4) {
      img <- img[,,Z:1,]
      index <- c(1,3,2,4)
    }
    img <- aperm(img, index) # re-organize orthogonal views
    pixdim <- pixdim[index[1:3]]
  }
  ## SAGITTAL
  if (is.sagittal(imageOrientationPatient)) {
    if (unlist(strsplit(first.row,""))[1] %in% c("H","F")) {
      if (ld == 3) {
        index <- c(2,1,3)
      }
      if (ld == 4) {
        index <- c(2,1,3,4)
      }
      img <- aperm(img, index)
      pixdim <- pixdim[index[1:3]]
    }
    if (unlist(strsplit(first.row,""))[1] == "P") {
      img <- switch(as.character(ld), "3" = img[X:1,,], "4" = img[X:1,,,],
                    stop("Dimension \"DIM\" incorrectly specified."))
    }
    if (unlist(strsplit(first.col,""))[1] == "H") {
      img <- switch(as.character(ld), "3" = img[,Y:1,], "4" = img[,Y:1,,])
    }
    ## The x-axis is increasing to the left hand side of the patient.
    z.index <- order(imagePositionPatient[,1])
    if (ld == 3) {
      img <- img[,,z.index]
      index <- c(3,1,2)
    }
    if (ld == 4) {
      img <- img[,,Z:1,]
      index <- c(3,1,2,4)
    }
    img <- aperm(img, index) # re-organize orthogonal views
    pixdim <- pixdim[index[1:3]]
  }
  imagePositionPatient <- imagePositionPatient[z.index,]
  if (any(is.na(imagePositionPatient))) {
    stop("Missing values are present in ImagePositionPatient.")
  }
  attr(img,"ipp") <- imagePositionPatient
  attr(img,"pixdim") <- pixdim
  return(img)
}

is.axial <- function(imageOrientationPatient, axial=c("L","R","A","P")) {
  first.row <- getOrientation(imageOrientationPatient[1,1:3])
  first.col <- getOrientation(imageOrientationPatient[1,4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% axial &
         unlist(strsplit(first.col, ""))[1] %in% axial)
}

is.coronal <- function(imageOrientationPatient,
                       coronal=c("L","R","H","F")) {
  first.row <- getOrientation(imageOrientationPatient[1,1:3])
  first.col <- getOrientation(imageOrientationPatient[1,4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% coronal &
         unlist(strsplit(first.col, ""))[1] %in% coronal)
}

is.sagittal <- function(imageOrientationPatient,
                        sagittal=c("A","P","H","F")) {
  first.row <- getOrientation(imageOrientationPatient[1,1:3])
  first.col <- getOrientation(imageOrientationPatient[1,4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    warning("Oblique acquisition in ImageOrientationPatient.")
  }
  return(unlist(strsplit(first.row, ""))[1] %in% sagittal &
         unlist(strsplit(first.col, ""))[1] %in% sagittal)
}
