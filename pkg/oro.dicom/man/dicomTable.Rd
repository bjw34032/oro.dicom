% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{dicomTable}
\alias{dicomTable}
\title{Construct Data Frame from DICOM Headers}
\description{
  A data frame is created given the valid DICOM fields provided by the
  user.
}
\usage{
dicomTable(hdrs, stringsAsFactors = FALSE, collapse = "-", colSort = TRUE,
           debug = FALSE)
}
\arguments{
  \item{hdrs}{List object of DICOM headers.}
  \item{stringsAsFactors}{Logical variable to be passed to
    \code{data.frame}.}
  \item{collapse}{Character string used to \code{paste} DICOM group,
    element and value fields.}
  \item{colSort}{Logical variable (default = \code{TRUE}) to sort column
    names in the table.}
  \item{debug}{Logical variable (default = \code{FALSE}) that regulates
    to display of intermediate processing steps.}
}
%\details{}
\value{
  Data frame where the rows correspond to images and the columns
  correspond to the UNION of all DICOM fields across all files in the
  list.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
%\seealso{}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
