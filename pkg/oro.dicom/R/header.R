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

dicomTable <- function(hdrs, stringsAsFactors=FALSE, collapse="-") {
  ## Use first record to establish data.frame
  csv <- data.frame(matrix(hdrs[[1]]$value, 1, nrow(hdrs[[1]])),
                    stringsAsFactors=stringsAsFactors)
  names(csv) <-
    paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[1]]$sequence)),
          as.vector(apply(hdrs[[1]][,1:3], 1, paste, collapse=collapse)),
          sep="")
  ## Loop through all records and "merge" them
  if (length(hdrs) > 1) {
    for (l in 2:length(hdrs)) {
      temp <- data.frame(matrix(hdrs[[l]]$value, 1, nrow(hdrs[[l]])),
                         stringsAsFactors=stringsAsFactors)
      names(temp) <-
        paste(sub("^-", "", gsub("[^0-9]+", "-", hdrs[[l]]$sequence)),
              as.vector(apply(hdrs[[l]][,1:3], 1, paste, collapse=collapse)),
              sep="")
      if (length(names(csv)) == length(names(temp)) &&
          all(names(csv) == names(temp))) {
        csv <- rbind(csv, temp)
      } else {
        csv <- merge(csv, temp, all=TRUE)
      }
    }
  }
  row.names(csv) <- names(hdrs)
  return(csv)
}

extractHeader <- function(hdrs, string, numeric=TRUE, names=FALSE) {
  if (is.data.frame(hdrs)) {
    hdrs <- list(hdrs)
  }
  out.list <- lapply(hdrs,
                     function(hdr) {
                       index <- which(hdr$name %in% string &
                                      hdr$sequence == "")
                       if(sum(index) > 0) {
                         hdr$value[index]
                       } else {
                         NA
                       }
                     })
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
        hwrite(str.warning, htmlfile, heading=3)
      } else {
        warning(str.warning)
      }
      return(expression(next))
    }
  }
  invisible()
}
