##
## Copyright (c) 2010-2011, Brandon Whitcher
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

## Sub-routines
.readCharWithEmbeddedNuls <- function(fid, n, to="UTF-8") {
  txt <- readBin(fid, "raw", n)
  iconv(rawToChar(txt[txt != as.raw(0)]), to=to)
}

.unsigned.header <- function(VR, implicit, fid, endian) {
  ## "Unsigned Long" and "Unsigned Short"
  length <- ifelse(implicit,
                   readBin(fid, integer(), size=4, endian=endian),
                   readBin(fid, integer(), size=2, endian=endian))
  value <- readBin(fid, integer(), n=length/VR$bytes,
                   size=VR$bytes, signed=FALSE, endian=endian)
  list(length=length, value=paste(value, collapse=" "))
}

.signed.header <- function(VR, implicit, fid, endian) {
  ## "Signed Long" and "Signed Short"
  length <- ifelse(implicit,
                   readBin(fid, integer(), size=4, endian=endian),
                   readBin(fid, integer(), size=2, endian=endian))
  value <- readBin(fid, integer(), n=length/VR$bytes,
                   size=VR$bytes, signed=TRUE, endian=endian)
  list(length=length, value=paste(value, collapse=" "))
}

.floating.header <- function(VR, implicit, fid, endian) {
  ## "Floating Point Single" and "Floating Point Double"
  length <- ifelse(implicit,
                   readBin(fid, integer(), size=4, endian=endian),
                   readBin(fid, integer(), size=2, endian=endian))
  value <- readBin(fid, numeric(), n=length/VR$bytes,
                   size=VR$bytes, signed=TRUE, endian=endian)
  list(length=length, value=value)
}

.other.header <- function(fid, implicit, endian) {
  ## "OtherByteString" or "OtherWordString"
  if (implicit) {
    length <- readBin(fid, integer(), size=4, endian=endian)
  } else {
    skip <- readBin(fid, integer(), size=2, endian=endian)
    length <- readBin(fid, integer(), size=4, endian=endian)
  }
  seek(fid, where=seek(fid) + length) # skip over this field
  list(length=length, value="skipped")
}

.sequence.header <- function(group, element, fid, implicit, endian,
                            skipSQ, SQ, EOS) {
  ## "Sequence of Items" with bytes = 0
  if (implicit) {
    length <- readBin(fid, integer(), size=4, endian=endian)
  } else {
    skip <- readBin(fid, integer(), size=2, endian=endian)
    length <- readBin(fid, integer(), size=4, endian=endian)
  }
  ## skip <- readBin(fid, integer(), size=2, endian=endian)
  ## length <- readBin(fid, integer(), size=4, endian=endian)
  if (length < 0) {
    ## Append (group,element) doublets for nested SequenceItem tags
    ## SQ <- paste(SQ, "(", group, ",", element, ")", sep="")
    SQ <- c(SQ, paste("(", group, ",", element, ")", sep=""))
  } else {
    if (skipSQ) {
      ## skip over all sequence items
      seek(fid, where=seek(fid) + max(0,length))
    } else {
      ## SQ <- paste(SQ, "(", group, ",", element, ")", sep="")
      SQ <- c(SQ, paste("(", group, ",", element, ")", sep=""))
      EOS <- c(EOS, seek(fid) + max(0,length))
    }
  }
  list(length=length, value="sequence", SQ=SQ, EOS=EOS)
}

.unknown.header <- function(VR, implicit, fid, endian, skipSQ) {
  ## Unknown header!
  if (VR$bytes > 0) {
    length <- ifelse(implicit,
                     readBin(fid, integer(), size=4, endian=endian),
                     readBin(fid, integer(), size=2, endian=endian))
    value <- .readCharWithEmbeddedNuls(fid, length)
    value <- sub(" +$", "", value) # remove white space at end
    value <- gsub("[\\]", " ", value) # remove "\\"s
  } else {
    if (implicit) {
      length <- readBin(fid, integer(), size=4, endian=endian)
    } else {
      skip <- readBin(fid, integer(), size=2, endian=endian)
      length <- readBin(fid, integer(), size=4, endian=endian)
    }
    if (skipSQ) {
      skip <- readBin(fid, integer(), max(0,length), size=1, endian=endian)
      value <- "skip"
    } else {
      value <- ""
    }
  }
  list(length=length, value=value)
}

.pixeldata.header <- function(hdr, implicit, fid, endian) {
  if (implicit) {
    length <- readBin(fid, integer(), size=4, endian=endian)
  } else {
    skip <- readBin(fid, integer(), size=2, endian=endian)
    length <- readBin(fid, integer(), size=4, endian=endian)
  }
  M <- which(hdr[, 3] == "Rows")[1]
  M <- as.numeric(hdr[M, 6])
  N <- which(hdr[, 3] == "Columns")[1]
  N <- as.numeric(hdr[N, 6])
  ## If the length is not provided, calculate from DICOM headers
  if (length < 0) {
    length <- M * N
  }
  bitsAllocated <- which(hdr[, 3] == "BitsAllocated" & nchar(hdr[, 7]) == 0)[1]
  bytes <- as.numeric(hdr[bitsAllocated, 6]) / 8
  ## Assuming only integer() data are being provided
  img <- readBin(fid, integer(), length, size=bytes, endian=endian)
  list(img=img, length=length, value="")
}

readDICOMFile <- function(fname, endian="little", flipud=TRUE, skip128=TRUE,
                          DICM=TRUE, skipSequence=TRUE, pixelData=TRUE,
                          warn=-1, debug=FALSE) {
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

  ## Warnings?
  oldwarn <- getOption("warn")
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
    if (readChar(fid, nchars=4) != "DICM") {
      stop("DICM != DICM")
    }
  }
  SQ <- EOS <- NULL
  hdr <- NULL
  pixel.data <- FALSE
  file.size <- file.info(fname)$size
  while (!pixel.data && !(seek(fid) >= file.size)) {
    seek.old <- seek(fid)
    implicit <- FALSE
    group <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    element <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    index <- which(dcm.group == group & dcm.element == element)
    name <- ifelse(any(index), dicom.dic$name[index], "Missing")
    vrstr <- .readCharWithEmbeddedNuls(fid, n=2)
    if (debug && length(index) < 1) {
      warning(sprintf("DICOM tag (%s,%s) is not in current dictionary.",
                      group, element))
    }
    if (! vrstr %in% VRcode) {
      ## Implicit VR (see diagram above)
      implicit <- TRUE
      seek(fid, where=seek(fid) - 2) # go back two bytes
      vrstr <- dicom.dic$code[index] # VR string from (group,element)
    }
    if (any(VRindex <- VRcode == vrstr)) {
      VR <- dicom.VR[VRindex, ]
    } else {
      VR <- dicom.VR[VRcode == "UN", ]
    }
    ## Check for PixelData (group,element) and _NOT_ a SequenceItem
    pixel.data <- (group == "7FE0" && element == "0010" &&
                   (is.null(SQ) || nchar(SQ) == 0))
    if (pixel.data && pixelData) {
      ## Read in the image data
      if (debug) {
        cat("##### Reading PixelData (7FE0,0010)", fill=TRUE)
      }
      out <- .pixeldata.header(hdr, implicit, fid, endian)
    } else {
      if (!is.null(SQ) && skipSequence &&
          (group == "FFFE" && element == "E000")) {
        out <- list(length=4, value="item")
        seek(fid, where=seek(fid) + out$length)
      } else {
        switch(VR$code,
               UL = ,
               US = { out <- .unsigned.header(VR, implicit, fid, endian) },
               SL = ,
               SS = { out <- .signed.header(VR, implicit, fid, endian) },
               FS = ,
               FD = { out <- .floating.header(VR, implicit, fid, endian) },
               OB = ,
               OW = { out <- .other.header(fid, implicit, endian) },
               SQ = {
                 out <- .sequence.header(group, element, fid, implicit,
                                         endian, skipSequence, SQ, EOS)
                 SQ <- out$SQ
                 EOS <- out$EOS
               },
               { out <- .unknown.header(VR, implicit, fid, endian,
                                        skipSequence) })
      }
    }
    if (is.null(SQ)) {
      hdr <- rbind(hdr, c(group, element, name, VR$code, out$length,
                          paste(out$value, collapse=" "), ""))
      if (debug) {
        cat("", seek.old, hdr[nrow(hdr), ], sep="\t", fill=TRUE)
      }
    } else {
      if (!skipSequence) {
        hdr <- rbind(hdr, c(group, element, name, VR$code, out$length,
                            paste(out$value, collapse=" "),
                            ifelse(is.null(SQ), "", paste(SQ, collapse=" "))))
        if (debug) {
          cat("", seek.old, hdr[nrow(hdr), ], sep="\t", fill=TRUE)
        }
      }
    }
    if (out$length > file.size) {
      warning(sprintf("DICOM tag (%s,%s) has length %d bytes which is greater than the file size (%d bytes).",
                      group, element, out$length, file.size))
    }
    if (name == "SequenceDelimitationItem" && !is.null(SQ)) {
      if ((lSQ <- length(SQ)) > 1) {
        SQ <- SQ[-lSQ] # remove the last (group,element) doublet
      } else {
        SQ <- NULL # set SQ to NULL
      }
    }
    if (!is.null(SQ)) {
      if (any(seek(fid) == EOS)) {
        SQ <- SQ[-which(seek(fid) == EOS)] # remove sequence(s)
        SQ <- eval(ifelse(length(SQ) < 1, expression(NULL), SQ)) # set SQ to NULL
        EOS <- EOS[-which(seek(fid) == EOS)] # remove endOfSequence(s)
        EOS <- eval(ifelse(length(EOS) < 1, expression(NULL), EOS)) # set EOS to NULL
      }
    }
  }
  close(fid)

  hdr <- as.data.frame(hdr, stringsAsFactors=FALSE)
  names(hdr) <- c("group", "element", "name", "code", "length",
                  "value", "sequence")
  hdr$length <- as.numeric(hdr$length)
  is.sequence <- nchar(hdr$sequence) > 0
  if (pixel.data && pixelData) {
    nr <- as.numeric(hdr$value[hdr$name == "Rows" & ! is.sequence])
    nc <- as.numeric(hdr$value[hdr$name == "Columns" & ! is.sequence])
    bytes <- as.numeric(hdr$value[hdr$name == "BitsAllocated" & ! is.sequence]) / 8
    length <- as.numeric(hdr$length[hdr$name == "PixelData" & ! is.sequence])
    total.bytes <- nr*nc*bytes
    if (total.bytes != length) {
      k <- length / total.bytes
      if (k == trunc(k)) {
        warning("3D DICOM file detected!")
        img <- array(out$img[1:(nc*nr*k)], c(nc,nr,k))
        img <- aperm(img, c(2,1,3))
        if (flipud) {
          img <- img[nr:1,,]
        }
      } else {
        stop("Number of bytes in PixelData does not match dimensions")
      }
    } else {
      img <- t(matrix(out$img[1:(nc*nr)], nc, nr))
      if (flipud) {
        img <- img[nr:1,]
      }
    }
  } else {
    img <- NULL
  }

  ## Warnings?
  options(warn=oldwarn)
  list(hdr=hdr, img=img)
}

readDICOM <- function(path, verbose=FALSE, counter=100, recursive=TRUE,
                      exclude=NULL, ...) {
  if (recursive) {
    filenames <- list.files(path, full.names=TRUE, recursive=TRUE)
  } else {
    filenames <- list.files(path, full.names=TRUE)
  }
  if (! is.null(exclude)) {
    filenames <- grep(exclude, filenames, ignore.case=TRUE, value=TRUE,
                      invert=TRUE)
  }
  nfiles <- length(filenames)
  nch <- nchar(as.character(nfiles))
  if (verbose) {
    cat("  ", nfiles, "files to be processed!", fill=TRUE)
  }
  headers <- images <- vector("list", nfiles)
  names(images) <- names(headers) <- filenames
  for (i in 1:nfiles) {
    if (verbose && (i %% counter == 0)) {
      cat("  ", sprintf(paste("% ", nch, "d", sep=""), i),
          "files processed...", fill=TRUE)
    }
    dcm <- readDICOMFile(filenames[i], ...)
    if (! is.null(dcm$img)) {
      images[[i]] <- dcm$img
    }
    if (! is.null(dcm$hdr)) {
      headers[[i]] <- dcm$hdr
    }
  }
  list(hdr=headers, img=images)
}

dicom2analyze <- function(dcm, datatype=4, reslice=TRUE, DIM=3,
                          descrip="SeriesDescription", ...) {
  require("oro.nifti")
  switch(as.character(DIM),
         "2" = { img <- create3D(dcm, ...) },
         "3" = { img <- create3D(dcm, ...) },
         "4" = { img <- create4D(dcm, ...) },
         stop("Dimension parameter \"DIM\" incorrectly specified."))
  if (DIM %in% 3:4 && reslice) {
    img <- swapDimension(img, dcm)
  }
  aim <- anlz(img, datatype=datatype)
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
  require("oro.nifti")
  switch(as.character(DIM),
         "2" = { img <- create3D(dcm, ...) },
         "3" = { img <- create3D(dcm, ...) },
         "4" = { img <- create4D(dcm, ...) },
         stop("Dimension parameter \"DIM\" incorrectly specified."))
  if (DIM %in% 3:4 && reslice) {
    img <- swapDimension(img, dcm)
  }
  nim <- nifti(img, datatype=datatype)
  if (is.null(attr(img,"pixdim"))) {
    ## (x,y) pixel dimensions
    pixelSpacing <- extractHeader(dcm$hdr, "PixelSpacing", FALSE)
    nim@"pixdim"[2:3] <- header2matrix(pixelSpacing, 2)[1,]
    ## z pixel dimensions
    nim@"pixdim"[4] <- ifelse(nim@"dim_"[1] > 2,
                              extractHeader(dcm$hdr, "SliceThickness")[1],
                              1)
  } else {
    nim@"pixdim"[2:4] <- attr(img,"pixdim")
  }
  ## description
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
  ## aux_file
  if (! is.null(aux.file)) {
    if (nchar(descrip.string) > 24) {
      warning("aux_file is greater than 24 characters and will be truncated")
    }
    nim@"aux_file" <- aux.file
  }
  ## units
  if (length(units) == 2) {
    nim@"xyzt_units" <- space.time2xyzt(units[1], units[2])
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

