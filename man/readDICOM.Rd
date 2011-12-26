% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{readDICOM}
\alias{readDICOM}
\alias{dicomSeparate}
\title{Read All DICOM Files in a Directory}
\description{
  All DICOM files are imported and a text file summarizing their content
  recorded.
}
\usage{
readDICOM(path, recursive = TRUE, exclude = NULL, verbose = FALSE,
          counter, ...)
}
\arguments{
  \item{path}{Path name to the DICOM directory.}
  \item{recursive}{Search recursively down from the given path name.}
  \item{exclude}{Exclude file names containing this character string.}
  \item{verbose}{Flag to provide text-based progress bar.}
  \item{counter}{Ignored.}
  \item{...}{Arguments to be passed to \code{readDICOMFile}.}
}
\details{
  A \code{for} loop is used to process each DICOM file contained in the
  directory(ies).  If only a single file is specified in the path,
  \code{readDICOM} will read that file only.
}
\value{
  A list structure with two major components:
  \item{img}{All images associated with the DICOM directory(ies).}
  \item{hdr}{All header files associated with the DICOM directory(ies).}
}
\references{
  Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with the
  DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
    Software}, \bold{44} (6), 1--28.
  \url{http://www.jstatsoft.org/v44/i06}

  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{readDICOMFile}}
}
\examples{
\dontrun{
dcmImage <- readDICOM(system.file("hk-40", package="oro.dicom"),
                      verbose=TRUE)
}
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
