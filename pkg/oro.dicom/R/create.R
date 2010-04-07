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

create3D <- function(dcm, mode="double", transpose=TRUE, pixelData=TRUE,
                     mosaic=FALSE, mosaicXY=NULL, ...) {
  if (pixelData) {
    if (is.null(dcm$hdr)) {
      stop("DICOM \"hdr\" information is not present.")
    }
    if (is.null(dcm$img)) {
      stop("DICOM \"img\" information is not present.")
    }
  } else {
    if (is.null(dcm$img)) {
      dcm <- list(hdr=dcm, img=NULL) # Only a list of headers as input
    }
  }
  X <- unique(extractHeader(dcm$hdr, "Rows"))
  if (length(X) != 1) {
    stop("Row lengths are not identical.")
  }
  Y <- unique(extractHeader(dcm$hdr, "Columns"))
  if (length(Y) != 1) {
    stop("Column lengths are not identical.")
  }
  if (mosaic) {
    if (is.null(dim(dcm$img)))
      stop("Multiple MOSAIC files detected, please use create4D()")
    if (is.null(mosaicXY)) {
      acquisitionMatrix <-
        header2matrix(extractHeader(dcm$hdr, "AcquisitionMatrix", FALSE), 4)
      x <- acquisitionMatrix[1,1]
      y <- acquisitionMatrix[1,4]
      if (is.na(x) || is.na(y))
        stop("Missing AcquisitionMatrix, please specify \"mosaicXY\"")
    } else {
      x <- mosaicXY[1]
      y <- mosaicXY[2]
    }
    z <- (X/x) * (Y/y)
    img <- array(0, c(x,y,z))
    k <- 1
    for (i in (X/x):1) {
      for (j in 1:(Y/y)) {
        img[,,k] <- dcm$img[((i-1)*x)+1:x, ((j-1)*y)+1:y]
        k <- k+1
      }
    }
    storage.mode(img) <- mode
  } else {
    ## Check if the DICOM list has length > 1
    Z <- ifelse(is.null(dim(dcm$img)), length(dcm$hdr), 1)
    img <- array(0, c(X,Y,Z))
    storage.mode(img) <- mode
    sliceLocation <- extractHeader(dcm$hdr, "SliceLocation")
    if (any(is.na(sliceLocation))) {
      stop("Missing values are present in SliceLocation.")
    }
    if (pixelData) {
      for (z in 1:Z) {
        z.order <- order(sliceLocation)[z]
        img[,,z] <- dcm$img[[z.order]]
      }
    } else {
      for (z in 1:Z) {
        z.order <- order(sliceLocation)[z]
        img[,,z] <- dicomInfo(names(dcm$hdr)[z.order])$img
      }
    }
  }
  sliceLocation <<- sliceLocation[order(sliceLocation)]
  if (transpose) {
    img <- aperm(img, c(2,1,3))
  }
#  patientPosition <- unique(extractHeader(dcm$hdr, "PatientPosition", FALSE))
#  if (length(patientPosition) != 1) {
#    stop("PatientPosition(s) are not identical.")
#  }
#  if (! mosaic) {
#    if (patientPosition == "FFS") {
#      sliceLocation <- sliceLocation[order(sliceLocation)]
#      img <- img[,,Z:1]
#      sliceLocation <<- rev(sliceLocation)
#    } else {
#  sliceLocation <<- sliceLocation[order(sliceLocation)]
#    }
#  }
  return(img)
}

create4D <- function(dcm, mode="double", transpose=TRUE, pixelData=TRUE,
                     mosaic=FALSE, mosaicXY=NULL, W=NULL, ...) {
  if (pixelData) {
    if (is.null(dcm$hdr)) {
      stop("DICOM \"hdr\" information is not present.")
    }
    if (is.null(dcm$img)) {
      stop("DICOM \"img\" information is not present.")
    }
  } else {
    if (is.null(dcm$img)) {
      dcm <- list(hdr=dcm, img=NULL) # Only a list of headers as input
    }
  }
  X <- unique(extractHeader(dcm$hdr, "Rows"))
  if (length(X) != 1) {
    stop("Row lengths are not identical.")
  }
  Y <- unique(extractHeader(dcm$hdr, "Columns"))
  if (length(Y) != 1) {
    stop("Column lengths are not identical.")
  }
  if (mosaic) {
    if (is.null(mosaicXY)) {
      acquisitionMatrix <-
        header2matrix(extractHeader(dcm$hdr, "AcquisitionMatrix", FALSE), 4)
      x <- acquisitionMatrix[1,1]
      y <- acquisitionMatrix[1,4]
      if (is.na(x) || x == 0 || is.na(y) || y == 0)
        stop("Missing AcquisitionMatrix, please specify \"mosaicXY\"")
    } else {
      x <- mosaicXY[1]
      y <- mosaicXY[2]
    }
    z <- (X/x) * (Y/y)
    W <- length(dcm$hdr)
    img <- array(0, c(x,y,z,W))
    for (w in 1:W) {
      k <- 1
      for (i in (X/x):1) {
        for (j in 1:(Y/y)) {
          img[,,k,w] <- dcm$img[[w]][((i-1)*x)+1:x, ((j-1)*y)+1:y]
          k <- k+1
        }
      }
    }
    storage.mode(img) <- mode
  } else {
    ## Check if the DICOM list has length > 1
    Z <- ifelse(is.null(dim(dcm$img)), length(dcm$hdr), 1)
    img <- array(0, c(X,Y,Z,W))
    storage.mode(img) <- mode
    sliceLocation <- extractHeader(dcm$hdr, "SliceLocation")
    if (any(is.na(sliceLocation))) {
      stop("Missing values are present in SliceLocation.")
    }
    if (pixelData) {
      for (z in 1:Z) {
        z.order <- order(sliceLocation)[z]
        img[,,z] <- dcm$img[[z.order]]
      }
    } else {
      for (z in 1:Z) {
        z.order <- order(sliceLocation)[z]
        img[,,z] <- dicomInfo(names(dcm$hdr)[z.order])$img
      }
    }
  }
  if (transpose) {
    img <- aperm(img, c(2,1,3,4))
  }
  patientPosition <- unique(extractHeader(dcm$hdr, "PatientPosition", FALSE))
  if (length(patientPosition) != 1) {
    stop("PatientPosition(s) are not identical.")
  }
  if (! mosaic) {
    if (patientPosition == "FFS") {
      sliceLocation <- sliceLocation[order(sliceLocation)]
      img <- img[,,Z:1,]
      sliceLocation <<- rev(sliceLocation)
    } else {
      sliceLocation <<- sliceLocation[order(sliceLocation)]
    }
  }
  return(img)
}
