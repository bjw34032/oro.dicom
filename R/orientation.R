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

swapDimension <- function(img, dcm) {
  patientPosition <- unique(extractHeader(dcm$hdr, "PatientPosition", FALSE))
  if (length(patientPosition) != 1) {
    stop("PatientPosition(s) are not identical.")
  }
  imageOrientationPatient <-
    header2matrix(extractHeader(dcm$hdr, "ImageOrientationPatient", FALSE), 6)
  ## Ensure all rows of imageOrientationPatient are identical!
  first.row <- getOrientation(imageOrientationPatient[1,1:3])
  first.col <- getOrientation(imageOrientationPatient[1,4:6])
  if (nchar(first.row) > 1 || nchar(first.col) > 1) {
    stop("Oblique acquisition in ImageOrientationPatient.")
  }
  X <- nrow(img)
  Y <- ncol(img)
  Z <- dim(img)[3]
  axial <- c("L","R","A","P")
  if (first.row %in% axial && first.col %in% axial) {
    cat("## Axial acquisition", fill=TRUE)
    if (first.row %in% c("A","P")) {
      img <- aperm(img, c(2,1,3))
    }
    if (first.row == "R") {
      img <- img[X:1,,]
    }
    if (first.col == "A") {
      img <- img[,Y:1,]
    }
    if ((patientPosition %in% c("FFS","FFP") &&
         sign(diff(sliceLocation))[1] < 0) ||
        (patientPosition %in% c("HFS","HFP") &&
         sign(diff(sliceLocation))[1] > 0)) {
      img <- img[,,Z:1]
      sliceLocation <<- rev(sliceLocation)
    }
  }
  coronal <- c("L","R","H","F")
  if (first.row %in% coronal && first.col %in% coronal) {
    cat("## Coronal acquisition", fill=TRUE)
    if (first.row %in% c("H","F")) {
      img <- aperm(img, c(2,1,3))
    }
    if (first.row == "R") {
      img <- img[X:1,,]
    }
    if (first.col == "H") {
      img <- img[,Y:1,]
    }
    if ((patientPosition %in% c("HFS","FFS") &&
         sign(diff(sliceLocation))[1] < 0) ||
        (patientPosition == c("HFP","FFP") &&
         sign(diff(sliceLocation))[1] > 0)) {
      img <- img[,,Z:1]
      sliceLocation <<- rev(sliceLocation)
    }
    img <- aperm(img, c(1,3,2)) # re-organize orthogonal views
  }
  sagittal <- c("A","P","H","F")
  if (first.row %in% sagittal && first.col %in% sagittal) {
    cat("## Sagittal acquisition", fill=TRUE)
      if (first.row %in% c("H","F")) {
      img <- aperm(img, c(2,1,3))
    }
    if (first.row == "P") {
      img <- img[X:1,,]
    }
    if (first.col == "H") {
      img <- img[,Y:1,]
    }
    img <- aperm(img, c(2,3,1)) # re-organize orthogonal views
  }
  ## sliceLocation <- extractHeader(dcm$hdr, "SliceLocation")
  if (any(is.na(sliceLocation))) {
    stop("Missing values are present in SliceLocation.")
  }
  return(img)
}
