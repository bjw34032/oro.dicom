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

dicom.VRfields <- function() {
  VR <- matrix(c("AE", "Application Entity", 16, 0,
                 "AS", "Age String", 4, 1, 
                 "AT", "Attribute Tag", 4, 1, 
                 "CS", "Code String", 16, 0, 
                 "DA", "Date", 8, 1, 
                 "DS", "Decimal String", 16, 0, 
                 "DT", "Date Time", 26, 0, 
                 "FL", "Floating Point Single", 4, 1, 
                 "FD", "Floating Point Double", 8, 1, 
                 "IS", "Integer String", 12, 0, 
                 "LO", "Long Strong", 64, 0, 
                 "LT", "Long Text", 10240, 0, 
                 "OB", "Other Byte String", 0, 0, 
                 "OW", "Other Word String", 0, 0, 
                 "PN", "Person Name", 64, 0, 
                 "SH", "Short String", 16, 0, 
                 "SL", "Signed Long", 4, 1, 
                 "SQ", "Sequence of Items", 0, 0, 
                 "SS", "Signed Short", 2, 1, 
                 "ST", "Short Text", 1024, 0, 
                 "TM", "Time", 16, 0, 
                 "UI", "Unique Identifier UID", 64, 0, 
                 "UL", "Unsigned Long", 4, 1, 
                 "UN", "Unknown", 0, 0, 
                 "US", "Unsigned Short", 2, 1, 
                 "UT", "Unlimited Text", 0, 0),
               ncol=4, byrow=TRUE)
  dimnames(VR) <- list(NULL, c("code", "name", "bytes", "fixed"))
  VR <- data.frame(VR, stringsAsFactors=FALSE)
  VR$bytes <- as.numeric(VR$bytes)
  VR$fixed <- as.numeric(VR$fixed)
  return(VR)
}

dicom.fields <- function() {
  data("dicom.fields.data")
  dicom.fields.data
}
