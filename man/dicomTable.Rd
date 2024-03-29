% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/header.R
\name{dicomTable}
\alias{dicomTable}
\title{Construct Data Frame from DICOM Headers}
\usage{
dicomTable(
  hdrs,
  stringsAsFactors = FALSE,
  collapse = "-",
  colSort = TRUE,
  verbose = FALSE,
  debug = FALSE
)
}
\arguments{
\item{hdrs}{List object of DICOM headers.}

\item{stringsAsFactors}{Logical variable to be passed to \code{data.frame}.}

\item{collapse}{Character string used to \code{paste} DICOM group, element
and value fields.}

\item{colSort}{Logical variable (default = \code{TRUE}) to sort column names
in the table.}

\item{verbose}{Flag to provide text-based progress bar (default =
\code{FALSE}).}

\item{debug}{Logical variable (default = \code{FALSE}) that regulates to
display of intermediate processing steps.}
}
\value{
Data frame where the rows correspond to images and the columns
correspond to the UNION of all DICOM fields across all files in the list.
}
\description{
A data frame is created given the valid DICOM fields provided by the user.
}
\references{
Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with
the DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
Software}, \bold{44} (6), 1--28.

Digital Imaging and Communications in Medicine (DICOM)\cr
\url{https://www.dicomstandard.org}
}
\author{
Brandon Whitcher \email{bwhitcher@gmail.com}
}
\keyword{misc}
