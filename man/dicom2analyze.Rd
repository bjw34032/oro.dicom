\name{dicomToAnalyze}
\alias{dicomToAnalyze}
\title{Convert DICOM Header to Analyze}
\description{
  A subset of header information from DICOM is placed into Analyze 7.5
  format.
}
\usage{
dicomToAnalyze(hdr, img)
}
\arguments{
  \item{hdr}{DICOM header}
  \item{img}{4D data volume}
}
\value{
  A list containing Analyze header information.
}
\details{
  See the references.
}
\references{
  Analyze 7.5\cr
  \url{http://www.mayo.edu/bir/PDF/ANALYZE75.pdf}
  
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
