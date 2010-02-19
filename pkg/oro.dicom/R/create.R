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

create3D <- function(dcm, X, Y, Z, sliceLocation, patientPosition="HFS",
                     mode="double", transpose=TRUE, pixelData=TRUE,
                     path=NULL) {
  if (pixelData) {
    if (is.null(dcm$hdr)) {
      stop("DICOM \"hdr\" information is not present.")
    }
    if (is.null(dcm$img)) {
      stop("DICOM \"img\" information is not present.")
    }
  } else {
    if (is.null(dcm$hdr)) {
      dcm <- list(hdr=dcm, img=NULL) # Only a list of headers as input
    }
  }
  img <- array(0, c(X,Y,Z))
  storage.mode(img) <- mode
  if (pixelData) {
    for (z in 1:Z) {
      z.order <- order(sliceLocation)[z]
      img[,,z] <- dcm$img[[z.order]]
    }
  } else {
    for (z in 1:Z) {
      z.order <- order(sliceLocation)[z]
      img[,,z] <- dicomInfo(file.path(path, names(dcm$hdr)[[z.order]]))$img
    }
  }
  if (transpose) {
    img <- aperm(img, c(2,1,3))
  }
  sliceLocation <<- sliceLocation[order(sliceLocation)]
  if (patientPosition == "FFS") {
    img <- img[,,Z:1]
    sliceLocation <<- rev(sliceLocation)
  }
  return(img)
}

