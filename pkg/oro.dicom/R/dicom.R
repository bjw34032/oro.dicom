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

dicomInfo <- function(fname, endian="little", flipud=TRUE, skip128=TRUE,
                      DICM=TRUE, pixelData=TRUE, warn=-1, debug=FALSE) {
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

  ## Sub-routines
  
  unsigned.header <- function(VR, implicit, fid, endian) {
    ## "Unsigned Long" and "Unsigned Short"
    length <- ifelse(implicit,
                     readBin(fid, integer(), size=4, endian=endian),
                     readBin(fid, integer(), size=2, endian=endian))
    value <- readBin(fid, integer(), n=length/VR$bytes,
                     size=VR$bytes, signed=FALSE, endian=endian)
    list(length=length, value=paste(value, collapse=" "))
  }
  
  signed.header <- function(VR, implicit, fid, endian) {
    ## "Signed Long" and "Signed Short"
    length <- ifelse(implicit,
                     readBin(fid, integer(), size=4, endian=endian),
                     readBin(fid, integer(), size=2, endian=endian))
    value <- readBin(fid, integer(), n=length/VR$bytes,
                     size=VR$bytes, signed=TRUE, endian=endian)
    list(length=length, value=paste(value, collapse=" "))
  }
  
  floating.header <- function(VR, implicit, fid, endian) {
    ## "Floating Point Single" and "Floating Point Double"
    length <- ifelse(implicit,
                     readBin(fid, integer(), size=4, endian=endian),
                     readBin(fid, integer(), size=2, endian=endian))
    value <- readBin(fid, numeric(), n=length/VR$bytes,
                     size=VR$bytes, signed=TRUE, endian=endian)
    list(length=length, value=value)
  }
  
  other.header <- function(fid, implicit, endian) {
    ## "OtherByteString" or "OtherWordString"
    if (implicit) {
      length <- readBin(fid, integer(), size=4, endian=endian)
    } else {
      skip <- readBin(fid, integer(), size=2, endian=endian)
      length <- readBin(fid, integer(), size=4, endian=endian)
    }
    seek(fid, where=seek(fid) + length) # skip over this field
    list(length=length, value <- "skipped")
  }
  
  sequence.header <- function(group, element, fid, endian) {
    ## "Sequence of Items" with bytes = 0
    skip <- readBin(fid, integer(), size=2, endian=endian)
    length <- readBin(fid, integer(), size=4, endian=endian)
    if (is.null(SQ) && length < 0) {
      SQ <<- paste("(", group, ",", element, ")", sep="")
    }
    seek(fid, where=seek(fid) + max(0,length)) # skip over all sequence items
    list(length=length, value="sequence")
  }
  
  null.header <- function(VR, group, element, fid, endian) {
    if (VR$bytes > 0) {
      if (! implicit) {
        length <- readBin(fid, integer(), size=2, endian=endian)
      }
      ## Trim trailing white space
      value <- sub(" +$", "",
                   iconv(rawToChar(readBin(fid, "raw", length)), to="UTF-8"))
      ## Replace all "\\"s with " "
      value <- gsub("[\\]", " ", value)
      ## Remove all non {a-zA-Z} characters with white space
      ## value <- gsub("[^{a-zA-Z}]", " ", value)
      if (VR$code == "UI") {
        value <- sub("\\0", "", value) # Remove trailing \0
      }
    } else {
      if (! implicit) {
        skip <- readBin(fid, integer(), size=2, endian=endian)
        length <- readBin(fid, integer(), size=4, endian=endian)
      }
      if (length >= 0 &&
          (VR$code != "SQ" && (group != "FFFE" &&
             (! element %in% c("E000","E00D","E0DD"))))) {
        skip <- readBin(fid, integer(), length, size=1, endian=endian)
        value <- "skip"
      } else {
        value <- "nothing matched"
      }
    }
    list(length=length, value=value)
  }
  
  unknown.header <- function(VR, implicit, fid, endian) {
    if (VR$bytes > 0) {
      length <- ifelse(implicit,
                       readBin(fid, integer(), size=4, endian=endian),
                       readBin(fid, integer(), size=2, endian=endian))
      value <- iconv(rawToChar(readBin(fid, "raw", length)), to="UTF-8")
      value <- sub(" +$", "", value) # remove white space at end
      value <- gsub("[\\]", " ", value) # remove "\\"s
    } else {
      if (implicit) {
        length <- readBin(fid, integer(), size=4, endian=endian)
      } else {
        skip <- readBin(fid, integer(), size=2, endian=endian)
        length <- readBin(fid, integer(), size=4, endian=endian)
      }
      skip <- readBin(fid, integer(), max(0,length), size=1, endian=endian)
      value <- "skip"
    }
    list(length=length, value=value)
  }

  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Some shortcuts...
  data("dicom.dic")
  dcm.group <- dicom.dic$group # [, "group"]
  dcm.element <- dicom.dic$element # [, "element"]
  data("dicom.VR")
  VRcode <- dicom.VR$code # [, "code"]
  ## Open connection
  fid <- file(fname, "rb")
  ## First 128 bytes are not used
  if (skip128) {
    seek(fid, where=128)
  }
  ## Next four bytes spell "DICM"
  if (DICM) {
    if (readChar(fid, n=4) != "DICM") {
      stop("DICM != DICM")
    }
  }
  hdr <- NULL
  pixel.data <- FALSE
  SQ <- NULL
  while (! pixel.data) {
    implicit <- FALSE
    group <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    element <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    pixel.data <- ifelse(group == "7FE0" & element == "0010", TRUE, FALSE)
    index <- which(dcm.group %in% group & dcm.element %in% element)
    vrstr <- iconv(rawToChar(readBin(fid, "raw", 2)), to="UTF-8") # readChar(fid, n=2)
    if (debug && length(index) != 1) {
      warning(sprintf("DICOM tag (%s,%s) is not in current dictionary.",
                      group, element))
    }
    if (! vrstr %in% VRcode) {
      ## Implicit VR (see diagram above)
      implicit <- TRUE
      seek(fid, where=seek(fid) - 2) # go back two bytes
      vrstr <- dicom.dic$code[index] # VR string from (group,element)
    }
    if (any(VRindex <- VRcode %in% vrstr)) {
      VR <- dicom.VR[VRcode %in% vrstr, ]
    } else {
      VR <- dicom.VR[VRcode %in% "UN", ]
    }
    if (pixel.data && pixelData) {
      if (debug) {
        cat("##### Reading PixelData (7FE0,0010)", fill=TRUE)
      }
      if (implicit) {
        length <- readBin(fid, integer(), size=4, endian=endian)
      } else {
        skip <- readBin(fid, integer(), size=2, endian=endian)
        length <- readBin(fid, integer(), size=4, endian=endian)
      }
      BitsAllocated <- which(hdr[, 3] %in% "BitsAllocated")[1]
      bytes <- as.numeric(hdr[BitsAllocated, 6]) / 8
      ## Assuming only integer() data are being provided
      img <- readBin(fid, integer(), length, size=bytes, endian=endian)
      out <- list(length=length, value=NULL)
    } else {
      out <- switch(VR$code,
                    UL = unsigned.header(VR, implicit, fid, endian),
                    US = unsigned.header(VR, implicit, fid, endian),
                    SL = signed.header(VR, implicit, fid, endian),
                    SS = signed.header(VR, implicit, fid, endian),
                    FS = floating.header(VR, implicit, fid, endian),
                    FD = floating.header(VR, implicit, fid, endian),
                    OB = other.header(fid, implicit, endian),
                    OW = other.header(fid, implicit, endian),
                    SQ = sequence.header(group, element, fid, endian),
                    unknown.header(VR, implicit, fid, endian))
    }
    name <- ifelse(any(index), dicom.dic$name[index], "Unknown")
    hdr <- rbind(hdr, c(group, element, name, VR$code, out$length,
                        out$value, ifelse(is.null(SQ), "", SQ)))
    if (debug) {
      cat("", seek(fid), group, element, name, VR$code, out$length,
          out$value, ifelse(is.null(SQ), "", SQ), sep="\t", fill=TRUE)
    }
    if (name == "SequenceDelimitationItem") {
      SQ <- NULL
    }
  }
  close(fid)

  hdr <- as.data.frame(hdr)
  names(hdr) <- c("group", "element", "name", "code", "length", "value")
  hdr$name <- as.character(hdr$name)
  hdr$length <- as.numeric(hdr$length)
  hdr$value <- as.character(hdr$value)

  if (pixelData) {
    nr <- as.numeric(hdr$value[hdr$name %in% "Rows"])
    nc <- as.numeric(hdr$value[hdr$name %in% "Columns"])
    img <- t(matrix(img[1:(nc*nr)], nc, nr))
    if (flipud) {
      img <- img[nr:1,]
    }
  } else {
    img <- NULL
  }

  ## Warnings?
  options(warn=oldwarn)
  list(hdr=hdr, img=img)
}

dicomSeparate <- function(path, verbose=FALSE, counter=100,
                          recursive=TRUE, exclude=NULL, ...) {
  if (recursive) {
    filenames <- list.files(path, full.names=TRUE, recursive=TRUE)
  } else {
    filenames <- list.files(path, full.names=TRUE)
  }
  if (! is.null(exclude)) {
    filenames <- grep(exclude, filenames, value=TRUE, invert=TRUE)
  }
  nfiles <- length(filenames)
  headers <- images <- vector("list", nfiles)
  names(images) <- names(headers) <- filenames
  for (i in 1:nfiles) {
    if (verbose && (i %% counter == 0)) {
      cat("  ", i, "files processed...", fill=TRUE)
    }
    dcm <- dicomInfo(filenames[i], ...)
    images[[i]] <- dcm$img
    headers[[i]] <- dcm$hdr
  }
  list(hdr=headers, img=images)
}

dicom2analyze <- function(img, hdr, descrip="SeriesDescription", ...) {
  require("oro.nifti")
  aim <- oro.nifti::anlz(img, ...)
  ## (x,y) pixel dimensions
  aim@"pixdim"[2:3] <- as.numeric(unlist(strsplit(extractHeader(hdr, "PixelSpacing", FALSE)[1], " ")))
  ## z pixel dimensions
  aim@"pixdim"[4] <- ifelse(aim@"dim_"[1] > 2,
                            extractHeader(hdr, "SliceThickness")[1],
                            1)
  ## description
  for (i in 1:length(descrip))
    if (i == 1) {
      descrip.string <- extractHeader(hdr, descrip[i], FALSE)[1]
    } else {
      descrip.string <- paste(descrip.string,
                              extractHeader(hdr, descrip[i], FALSE)[1],
                              sep="; ")
    }
  if (nchar(descrip.string) > 80) {
    warning("Description is greater than 80 characters and has been truncated")
    aim@"descrip" <- substring(descrip.string, 1, 80)
  } else {
    aim@"descrip" <- descrip.string
  }
  ## originator
  aim$"originator" <- substring(extractHeader(hdr, "RequestingPhysician")[1],
                                1, 10)
  ## scannum
  aim@"scannum" <- substring(extractHeader(hdr, "StudyID")[1], 1, 10)
  ## patient_id
  aim@"patient_id" <- substring(extractHeader(hdr, "PatientID")[1], 1, 10)
  ## exp_date
  aim@"exp_date" <- substring(extractHeader(hdr, "StudyDate")[1], 1, 10)
  ## exp_time
  aim@"exp_time" <- substring(extractHeader(hdr, "StudyTime")[1], 1, 10)
  return(aim)
}

dicom2nifti <- function(img, hdr, units=c("mm","sec"), rescale=FALSE, 
                        descrip="SeriesDescription", ...) {
  require("oro.nifti")
  nim <- oro.nifti::nifti(img, ...)
  ## (x,y) pixel dimensions
  nim@"pixdim"[2:3] <- as.numeric(unlist(strsplit(extractHeader(hdr, "PixelSpacing", FALSE)[1], " ")))
  ## z pixel dimensions
  nim@"pixdim"[4] <- ifelse(nim@"dim_"[1] > 2,
                            extractHeader(hdr, "SliceThickness")[1],
                            1)
  ## description
  for (i in 1:length(descrip))
    if (i == 1) {
      descrip.string <- extractHeader(hdr, descrip[i], FALSE)[1]
    } else {
      descrip.string <- paste(descrip.string,
                              extractHeader(hdr, descrip[i], FALSE)[1],
                              sep="; ")
    }
  if (nchar(descrip.string) > 80)
    warning("Description is greater than 80 characters and will be truncated")
  nim@"descrip" <- descrip.string
  ## units
  if (length(units) == 2) {
    nim@"xyzt_units" <- space.time2xyzt(units[1], units[2])
  } else {
    stop("units must be a length=2 vector")
  }
  ## rescaling (more of a CT thing?)
  if (rescale) {
    nim@"scl_slope" <- extractHeader(hdr, "RescaleSlope")[1]
    nim@"scl_inter" <- extractHeader(hdr, "RescaleIntercept")[1]
  }
  return(nim)
}

