% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{readDICOM}
\alias{readDICOM}
\title{Read All DICOM Files in a Directory}
\description{
  All DICOM files are imported and a text file summarizing their content
  recorded.
}
\usage{
readDICOM(path, verbose = FALSE, counter = 100, recursive = TRUE,
          exclude = NULL, ...)
}
\arguments{
  \item{path}{Path name to the DICOM directory.}
  \item{verbose}{Flag to provide limited progress information.}
  \item{counter}{Frequency of updating the diaplay about the progress of
    file processing.}
  \item{recursive}{Search recursively down from the given path name.}
  \item{exclude}{Exclude file names containing this character string.}
  \item{...}{Arguments to be passed to \code{readDICOMFile}.}
}
\details{
  A \code{for} loop is used to process each DICOM file contained in the
  directory(ies).
}
\value{
  A list structure with two major components:
  \item{img}{All images associated with the DICOM directory(ies).}
  \item{hdr}{All header files associated with the DICOM directory(ies).}
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{readDICOMFile}}
}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
