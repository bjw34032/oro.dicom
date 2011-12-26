##
## Copyright (c) 2010-2011 Brandon Whitcher
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

dicomTable <- function(hdrs, stringsAsFactors=FALSE, collapse="-",
                       colSort=TRUE, verbose=FALSE, debug=FALSE) {
  myMerge <- function(df1, df2) {
    if (anyDuplicated(names(df1)) != 0) {
      warning("Duplicated group-element tags have been removed!")
      df1 <- df1[, ! duplicated(names(df1))]
    }
    if (! all(names(df2) %in% names(df1))) {
      newCols <- names(df2)[! names(df2) %in% names(df1)]
      ## newcols <- setdiff(names(df2), names(df1)) # removes duplicates!
      newDf <- as.data.frame(lapply(newCols, function(i, x) rep(NA, x),
                                    x = nrow(df1)))
      names(newDf) <- newCols
      df1 <- cbind(df1, newDf)
    }
    if (anyDuplicated(names(df2)) != 0) {
      warning("Duplicated group-element tags have been removed!")
      df2 <- df2[, ! duplicated(names(df2))]
    }
    if (! all(names(df1) %in% names(df2))) {
      newCols <- names(df1)[! names(df1) %in% names(df2)]
      ## newCols <- setdiff(names(df1), names(df2)) # removes duplicates!
      newDf <- as.data.frame(lapply(newCols, function(i, x) rep(NA, x),
                                    x = nrow(df2)))
      names(newDf) <- newCols
      df2 <- cbind(df2, newDf)
    }
    rbind(df1, df2)
  }
  ## Use first record to establish data.frame
  csv <- data.frame(matrix(hdrs[[1]]$value, 1, nrow(hdrs[[1]])),
                    stringsAsFactors=stringsAsFactors)
  names(csv) <-
    paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[1]]$sequence)),
          as.vector(apply(hdrs[[1]][,1:3], 1, paste, collapse=collapse)),
          sep="")
  ## Loop through all records and "merge" them
  if ((nhdrs <- length(hdrs)) > 1) {
    if (verbose) {
      cat(" ", nhdrs, "files to be processed by dicomTable()", fill=TRUE)
      tpb <- txtProgressBar(min=0, max=nhdrs, style=3)
    }
    for (l in 2:nhdrs) {
      if (debug) {
        cat("  l =", l, fill=TRUE)
      }
      if (verbose) {
        setTxtProgressBar(tpb, l)
      }
      temp <- data.frame(matrix(hdrs[[l]]$value, 1, nrow(hdrs[[l]])),
                         stringsAsFactors=stringsAsFactors)
      names(temp) <-
        paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[l]]$sequence)),
              as.vector(apply(hdrs[[l]][,1:3], 1, paste, collapse=collapse)),
              sep="")
      old.nrow <- nrow(csv)
      csv <- myMerge(csv, temp)
      if (nrow(csv) == old.nrow) {
        warning("Duplicate row was _not_ inserted in data.frame (csv)")
        csv <- rbind(csv, NA)
      }
    }
    if (verbose) {
      close(tpb)
    }
    row.names(csv) <- names(hdrs)
  }
  if (colSort) {
    return(csv[, order(names(csv))])
  } else {
    return(csv)
  }
}

extractHeader <- function(hdrs, string, numeric=TRUE, names=FALSE,
                          inSequence=TRUE) {
  if (is.data.frame(hdrs)) {
    hdrs <- list(hdrs)
  }
  out.list <- lapply(hdrs,
                     function(hdr, string, inSequence) {
                       if (inSequence) {
                         sequence <- FALSE
                       } else {
                         sequence <- nchar(hdr$sequence) > 0
                       }
                       index <- which(hdr$name %in% string & !sequence)
                       if (sum(index) > 0) {
                         hdr$value[index]
                       } else {
                         NA
                       }
                     }, string=string, inSequence=inSequence)
  out.names <- names(out.list)
  out.vec <- unlist(out.list)
  if (numeric) {
    out.vec <- as.numeric(out.vec)
  }
  if (names) {
    names(out.vec) <- out.names
  } else {
    out.vec <- as.vector(out.vec)
  }
  return(out.vec)
}

header2matrix <- function(hdr, ncol, sep=" ", byrow=TRUE) {
  matrix(as.numeric(unlist(strsplit(hdr, sep))), ncol=ncol, byrow=byrow)
}

matchHeader <- function(hdr, string) {
  ifelse(is.na(hdr), FALSE, regexpr(string, hdr, ignore.case=TRUE) > -1)
}

writeHeader <- function(dtable, filename, ...) {
  write.table(dtable, filename, quote=FALSE, sep="\t", ...)
}

nextHeader <- function(dcm, string, reference, str.warning,
                       htmlfile=NULL, heading=3, numeric=FALSE) {
  header <- extractHeader(dcm$hdr, string=string, numeric=numeric)
  for (i in 1:length(reference)) {
    if (any(matchHeader(header, string=reference[i]))) {
      if (! is.null(htmlfile)) {
        require("hwriter")
        hwrite(str.warning, htmlfile, heading=3)
      } else {
        warning(str.warning)
      }
      return(expression(next))
    }
  }
  invisible()
}
