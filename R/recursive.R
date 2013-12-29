##
## Copyright (c) 2010-2013, Brandon Whitcher
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

rereadDICOMFile <- function(fname, endian="little", flipud=TRUE, DICM=TRUE,
                            skipSequence=FALSE, readPixelData=TRUE, warn=-1,
                            debug=FALSE) {
    ## Warnings?
    oldwarn <- getOption("warn")
    options(warn = warn)
    ##
    fsize <- file.info(fname)$size
    fraw <- readBin(fname, "raw", n=as.integer(fsize), endian=endian)
    skip128 <- fraw[1:128]
    if (debug) {
        cat("#", "First 128 bytes of DICOM header =", fill=TRUE)
        print(skip128)
    }
    if (DICM) {
        if (rawToChar(fraw[129:132]) != "DICM") {
            stop("DICM != DICM")
        }
    }
    dicomHeader <- sequence <- NULL
    seq.txt <- ""
    ## Call parseDICOMHeader() to recursively parse the DICOM header
    dcm <- parseDICOMHeader(fraw[133:fsize], seq.txt, endian=endian, verbose=debug)
    hdr <- as.data.frame(dcm$header, stringsAsFactors=FALSE)
    row.names(hdr) <- NULL
    names(hdr) <- c("group", "element", "name", "code", "length", "value", "sequence")
    ##
    if (dcm$pixel.data && readPixelData) {
        if (debug) {
            cat("##### Reading PixelData (7FE0,0010) #####", fill=TRUE)
        }
        img <- parsePixelData(fraw[(132 + dcm$data.seek + 1):fsize], hdr, endian, flipud)
    } else { 
        if (dcm$spectroscopy.data && readPixelData) {
            if (debug) {
                cat("##### Reading SpectroscopyData (5600,0020) #####", fill=TRUE)
            }
            img <- parseSpectroscopyData(fraw[(132 + dcm$data.seek + 1):fsize], hdr, endian)
        } else {
            img <- NULL
        }
    }
    ## Warnings?
    options(warn = oldwarn)
    ##
    list(hdr = hdr, img = img)
}

.rawToCharWithEmbeddedNuls <- function(str.raw, to="UTF-8") {
    iconv(rawToChar(str.raw[str.raw != as.raw(0)]), to=to)
}

parseDICOMHeader <- function(rawString, sq.txt="", endian="little", 
                             verbose=FALSE) {
    ##
    ## "The default DICOM Transfer Syntax, which shall be supported by
    ## all AEs, uses Little Endian encoding and is specified in Annex
    ## A.1." (PS 3.5-2004, page 38)
    ##
    ## PS 3.5-2004, Sect 7.1.2: Data Element Structure with Explicit VR
    ## Explicit VRs store VR as text chars in 2 bytes.
    ## VRs of OB, OW, SQ, UN, UT have VR chars, then 0x0000, then 32 bit VL:
    ##
    ## +-----------------------------------------------------------+
    ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 |
    ## +----+----+----+----+----+----+----+----+----+----+----+----+
    ## |<Group-->|<Element>|<VR----->|<0x0000->|<Length----------->|<Value->
    ##
    ## Other Explicit VRs have VR chars, then 16 bit VL:
    ##
    ## +---------------------------------------+
    ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
    ## +----+----+----+----+----+----+----+----+
    ## |<Group-->|<Element>|<VR----->|<Length->|<Value->
    ##
    ## Implicit VRs have no VR field, then 32 bit VL:
    ##
    ## +---------------------------------------+
    ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
    ## +----+----+----+----+----+----+----+----+
    ## |<Group-->|<Element>|<Length----------->|<Value->
    ##
    rawToHex <- function(bytes) {
        toupper(paste(rev(bytes), collapse=""))
    }
    is.item <- function (group, element) {
        group == "FFFE" && element %in% c("E000","E00D","E0DD")
    }
    strseek <- dseek <- 0
    dicomHeader <- NULL
    pixelData <- spectroscopyData <- FALSE
    while (! (pixelData || spectroscopyData) && strseek < length(rawString)) { ##
        rm(group, element, dictionaryIndex, dic, rawValue, VR, vr, value, length)
        group <- rawToHex(rawString[strseek + 1:2])
        element <- rawToHex(rawString[strseek + 3:4])
        if (! any(dictionaryIndex <- group == dicom.dic$group & element == dicom.dic$element)) {
            ## Private tag = Unknown
            dic <- data.frame(group = group, element = element, code = "UN", 
                              offset = 1, name = "Unknown", stringsAsFactors=FALSE)
        } else {
            dic <- dicom.dic[dictionaryIndex,]
        }
        if (verbose) {
            cat("#", group, element, dic$name, dic$code, sep="\t")
        }
        b56 <- .rawToCharWithEmbeddedNuls(rawString[strseek + 5:6])
        b78 <- readBin(rawString[strseek + 7:8], "integer", size=2, endian=endian)
        b47 <- readBin(rawString[strseek + 5:8], "integer", size=4, endian=endian)
        strseek <- strseek + 8
        if (b56 %in% c("OB","OW","SQ","UN","UT") && ! is.item(group, element)) {
            ## Explicit VR
            length <- readBin(rawString[strseek + 1:4], "integer", size=4, endian=endian)
            strseek <- strseek + 4
            if (b56 != "SQ") {
                rawValue <- rawString[strseek + 1:length]
            }
            vr <- b56
        } else {
            if (b56 %in% c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                           "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
                && ! is.item(group, element)) {
                ## Explicit VR
                length <- b78
                if (dic$code != "SQ") {
                    rawValue <- rawString[strseek + 1:length]
                }
                vr <- b56
            } else {
                ## Implicit VR
                length <- b47
                vr <- NULL
                if (is.item(group, element)) {
                    length <- 0
                    rawValue <- raw(0)
                } else {
                    if (dic$code != "SQ") {
                        rawValue <- rawString[strseek + 1:length]
                    } else {
                        vr <- "SQ"
                    }
                }
            }
        }
        if (! is.null(vr)) {
            VR <- dicom.VR[dicom.VR$code == vr, ]
        } else {
            VR <- dicom.VR[dicom.VR$code == "UN", ]
        }
        if (verbose) {
            cat("", VR$code, length, sep="\t")
        }
        if (group == "7FE0" && element == "0010" && sq.txt == "") { 
            ## PixelData
            value <- "PixelData"
            pixelData <- TRUE
            dseek <- strseek
        } else {
            if (group == "5600" && element == "0020" && sq.txt == "") { 
                ## SpectroscopyData
                value <- "SpectroscopyData"
                spectroscopyData <- TRUE
                dseek <- strseek + 4 # HACK: not sure why I need to skip an extra four bytes
            } else {
                if (dic$code != VR$code) { 
                    ## Should this be a stop()?
                    warning("DICOM dictionary does not match VR code")
                }
                value <- switch(VR$code,
                                UL = ,
                                US = readBin(rawValue, "integer", n=length/VR$bytes, 
                                             size=VR$bytes, signed=FALSE, endian=endian),
                                SL = ,
                                SS = readBin(rawValue, "integer", n=length/VR$bytes, 
                                             size=VR$bytes, signed=TRUE, endian=endian),
                                FD = ,
                                FL = readBin(rawValue, "numeric", n=length/VR$bytes, 
                                             size=VR$bytes, signed=TRUE, endian=endian),
                                OB = ,
                                OW = .rawToCharWithEmbeddedNuls(rawValue),
                                SQ = "Sequence",
                                {
                                    if (length > 0) {
                                        tmpString <- .rawToCharWithEmbeddedNuls(rawValue)
                                        tmpString <- sub(" +$", "", tmpString)     # remove white space at end
                                        tmpString <- gsub("[\\/]", " ", tmpString) # remove "/"s
                                        tmpString <- gsub("[\\^]", " ", tmpString) # remove "^"s
                                    } else {
                                        tmpString <- ""
                                    }
                                    tmpString
                                })
            }
        }
        if (verbose) {
            cat("", value, sq.txt, sep="\t", fill=TRUE)
        }
        dicomHeaderRow <- c(group, element, dic$name, dic$code, length, value, sq.txt)
        dicomHeader <- rbind(dicomHeader, dicomHeaderRow)
        if (group == "FFFE" && element == "E0DD") { # SequenceDelimitationItem
            dseek <- strseek
            break
        }
        if (VR$code == "SQ") {
            groupElement <- paste("(", group, ",", element, ")", sep="")
            if (length > 0) {
                ## Pass length of bytes provided explicitly by the sequence tag
                dcm <- parseDICOMHeader(rawString[strseek + 1:length], 
                                        paste(sq.txt, groupElement), 
                                        verbose=verbose)
            } else {
                ## Pass remaining bytes and look for SequenceDelimitationItem tag
                dcm <- parseDICOMHeader(rawString[(strseek + 1):length(rawString)], 
                                        paste(sq.txt, groupElement),
                                        verbose=verbose)
                length <- dcm$data.seek
            }
            dicomHeader <- rbind(dicomHeader, dcm$header)
        }
        strseek <- strseek + ifelse(length >= 0, length, 0)
    } ##
    list(header = dicomHeader, pixel.data = pixelData, data.seek = dseek,
         spectroscopy.data = spectroscopyData)
}

parsePixelData <- function(rawString, hdr, endian="little", flipupdown=TRUE) {
    rows <- as.numeric(with(hdr, value[name == "Rows" & sequence == ""]))
    columns <- as.numeric(with(hdr, value[name == "Columns" & sequence == ""]))
    bytes <- as.numeric(with(hdr, value[name == "BitsAllocated" & sequence == ""])) / 8
    length <- as.numeric(with(hdr, length[name == "PixelData" & sequence == ""]))
    if (length <= 0) {
        stop("Number of bytes in PixelData not specified")
    }
    imageData <- readBin(rawString[1:length], "integer", n=length, size=bytes,
                         endian=endian)
    if (length == rows * columns * bytes) { # 2D PixelData
        imageData <-  t(matrix(imageData[1:(columns * rows)], columns, rows))
        if (flipupdown) {
            imageData <- imageData[rows:1, ]
        }
    } else { # 3D PixelData
        k <- length / rows / columns / bytes
        if (k == trunc(k)) {
            warning("3D DICOM file detected!")
            imageData <- array(imageData[1:(columns * rows * k)], c(columns, rows, k))
            imageData <- aperm(imageData, c(2,1,3))
            if (flipupdown) {
                imageData <- imageData[rows:1, , ]
            }
        } else {
            stop("Number of bytes in PixelData does not match dimensions of image")
        }
    }
    return(imageData)
}

parseSpectroscopyData <- function(rawString, hdr, endian="little") {
    numberOfFrames <- as.numeric(with(hdr, value[name == "NumberOfFrames" & sequence == ""]))
    rows <- as.numeric(with(hdr, value[name == "Rows" & sequence == ""]))
    columns <- as.numeric(with(hdr, value[name == "Columns" & sequence == ""]))
    dataPointRows <- as.numeric(with(hdr, value[name == "DataPointRows" & sequence == ""]))
    dataPointColumns <- as.numeric(with(hdr, value[name == "DataPointColumns" & sequence == ""]))
    dataRepresentation <- with(hdr, value[name == "DataRepresentation" & sequence == ""])
    nComponents <- ifelse(dataRepresentation == "COMPLEX", 2, 1)
    whichComponent <- ifelse(nComponents > 1 && dataRepresentation == "IMAGINARY", 2, 1)
    valuesPerFrame <- columns * rows * dataPointRows * dataPointColumns
    bytes <- as.numeric(with(hdr, value[name == "BitsAllocated" & sequence == ""])) / 8
    length <- valuesPerFrame * numberOfFrames
    imageData <- readBin(rawString[1:length], "numeric", n=length, size=4, endian=endian)
    odd <- seq(1, length, by=2)
    even <- seq(2, length, by=2)
    complex(real=imageData[odd], imaginary=imageData[even])
}
