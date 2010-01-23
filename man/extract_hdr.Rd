\name{extractHeader}
\alias{extractHeader}
\title{Extract Single Field from DICOM Headers}
\description{
  A particular DICOM field is extracted for a collection of DICOM
  headers.
}
\usage{
extractHeader(hdrs, string, numeric=TRUE, names=FALSE)
}
\arguments{
  \item{hdrs}{List object of DICOM headers.}
  \item{string}{DICOM field name.}
  \item{numeric}{Logical; values are converted to numbers when
    \code{TRUE}.}
  \item{names}{Logical; file names are kept with elements of the
    vector.}
}
\details{
  The DICOM field is extracted from each DICOM header and placed into a
  vector.
}
\value{
  Vector of values from the requested DICOM field.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{dicomInfo}}
}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
