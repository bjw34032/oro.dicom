##
## Copyright (c) 2010-2014, Brandon Whitcher
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

dicom2analyze <- function(dcm, datatype=4, reslice=TRUE, DIM=3,
                          descrip="SeriesDescription", ...) {
    switch(as.character(DIM),
           "2" = { 
               dcmList <- list(hdr=list(dcm$hdr), img=list(dcm$img))
               img <- create3D(dcmList, ...)
               },
           "3" = { img <- create3D(dcm, ...) },
           "4" = { img <- create4D(dcm, ...) },
           stop("Dimension parameter \"DIM\" incorrectly specified."))
    if (DIM %in% 3:4 && reslice) {
        img <- swapDimension(img, dcm)
    }
    aim <- oro.nifti::anlz(img, datatype=datatype)
    if (is.null(attr(img,"pixdim"))) {
        ## (x,y) pixel dimensions
        aim@"pixdim"[2:3] <- as.numeric(unlist(strsplit(extractHeader(dcm$hdr, "PixelSpacing", FALSE)[1], " ")))
        ## z pixel dimensions
        aim@"pixdim"[4] <- ifelse(aim@"dim_"[1] > 2,
                                  extractHeader(dcm$hdr, "SliceThickness")[1],
                                  1)
    } else {
        aim@"pixdim"[2:4] <- attr(img,"pixdim")
    }
    ## description
    for (i in 1:length(descrip))
        if (i == 1) {
            descrip.string <- extractHeader(dcm$hdr, descrip[i], FALSE)[1]
        } else {
            descrip.string <- paste(descrip.string,
                                    extractHeader(dcm$hdr, descrip[i], FALSE)[1],
                                    sep="; ")
        }
    if (nchar(descrip.string) > 80) {
        warning("Description is greater than 80 characters and has been truncated")
        aim@"descrip" <- substring(descrip.string, 1, 80)
    } else {
        aim@"descrip" <- descrip.string
    }
    ## scannum
    aim@"scannum" <- unlist(substring(extractHeader(dcm$hdr, "StudyID")[1], 1, 10))
    ## patient_id
    aim@"patient_id" <- substring(extractHeader(dcm$hdr, "PatientID")[1], 1, 10)
    ## exp_date
    aim@"exp_date" <- substring(extractHeader(dcm$hdr, "StudyDate")[1], 1, 10)
    ## exp_time
    aim@"exp_time" <- substring(extractHeader(dcm$hdr, "StudyTime")[1], 1, 10)
    return(aim)
}

dicom2nifti <- function(dcm, datatype=4, units=c("mm","sec"), rescale=FALSE,
                        reslice=TRUE, qform=TRUE, sform=TRUE, DIM=3,
                        descrip="SeriesDescription", aux.file=NULL, ...) {
    switch(as.character(DIM),
           "2" = { 
               dcmList <- list(hdr=list(dcm$hdr), img=list(dcm$img))
               img <- create3D(dcmList, ...)
               },
           "3" = { img <- create3D(dcm, ...) },
           "4" = { img <- create4D(dcm, ...) },
           stop("Dimension parameter \"DIM\" incorrectly specified."))
    if (DIM %in% 3:4 && reslice) {
        img <- swapDimension(img, dcm)
    }
    nim <- oro.nifti::nifti(img, datatype=datatype)
    if (is.null(attr(img, "pixdim"))) {
        ## (x,y) pixel dimensions
        pixelSpacing <- extractHeader(dcm$hdr, "PixelSpacing", FALSE)
        nim@"pixdim"[2:3] <- header2matrix(pixelSpacing, 2)[1,]
        ## z pixel dimensions
        nim@"pixdim"[4] <- ifelse(nim@"dim_"[1] > 2,
                                  extractHeader(dcm$hdr, "SliceThickness")[1], 1)
    } else {
        nim@"pixdim"[2:4] <- attr(img,"pixdim")
    }
    ## description
    if (! is.null(descrip)) {
        for (i in 1:length(descrip)) {
            if (i == 1) {
                descrip.string <- extractHeader(dcm$hdr, descrip[i], FALSE)[1]
            } else {
                descrip.string <- paste(descrip.string,
                                        extractHeader(dcm$hdr, descrip[i], FALSE)[1],
                                        sep="; ")
            }
        }
        if (nchar(descrip.string) > 80) {
            warning("Description is greater than 80 characters and will be truncated")
        }
        nim@"descrip" <- descrip.string
    }
    ## aux_file
    if (! is.null(aux.file)) {
        if (nchar(descrip.string) > 24) {
            warning("aux_file is greater than 24 characters and will be truncated")
        }
        nim@"aux_file" <- aux.file
    }
    ## units
    if (length(units) == 2) {
        nim@"xyzt_units" <- oro.nifti::space.time2xyzt(units[1], units[2])
    } else {
        stop("units must be a length = 2 vector")
    }
    ## qform
    if (qform) { # Basic LAS convention corresponds to the xform matrix
        nim@"qform_code" <- 2
        nim@"quatern_b" <- 0
        nim@"quatern_c" <- 1
        nim@"quatern_d" <- 0
        nim@"pixdim"[1] <- -1.0 # qfac
    }
    ## sform
    if (sform) { # Basic LAS convention corresponds to the xform matrix
        nim@"sform_code" <- 2
        nim@"srow_x" <- pixdim(nim)[2] * c(-1,0,0,0)
        nim@"srow_y" <- pixdim(nim)[3] * c(0,1,0,0)
        nim@"srow_z" <- pixdim(nim)[4] * c(0,0,1,0)
    }
    ## rescale
    if (rescale) {
        nim@"scl_slope" <- extractHeader(dcm$hdr, "RescaleSlope")[1]
        nim@"scl_inter" <- extractHeader(dcm$hdr, "RescaleIntercept")[1]
    }
    return(nim)
}
