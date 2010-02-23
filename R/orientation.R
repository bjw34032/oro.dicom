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

getOrientation <- function(xyz) {
  oX <- ifelse(xyz[1] < 0, "R", "L")
  oY <- ifelse(xyz[2] < 0, "A", "P")
  oZ <- ifelse(xyz[3] < 0, "F", "H")
  aX <- abs(xyz[1])
  aY <- abs(xyz[2])
  aZ <- abs(xyz[3])

  orientation <- NULL
  epsilon <- 0.0001
  for (i in 1:3) {
    if (aX > epsilon && aX > aY && aX > aZ) {
      orientation <- paste(orientation, oX, sep="")
      aX <- 0
    } else {
      if (aY > epsilon && aY > aX && aY > aZ) {
        orientation <- paste(orientation, oY, sep="")
        aY <- 0
      } else {
        if (aZ > epsilon && aZ > aX && aZ > aY) {
          orientation <- paste(orientation, oZ, sep="")
          aZ <- 0
        }
      }
    }
  }
  return(orientation)
}

