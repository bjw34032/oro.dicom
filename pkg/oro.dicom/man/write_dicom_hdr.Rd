\name{writeDicomHeader}
\alias{writeDicomHeader}
\title{Write DICOM Table to ASCII File}
\description{
  A wrapper to \code{write.table} specifically for DICOM tables.
}
\usage{
writeDicomHeader(dtable, filename, ...)
}
\arguments{
  \item{dtable}{The DICOM table.}
  \item{filename}{Name of the file to be created.}
  \item{...}{Additional parameters to be passed to \code{write.table}.}
}
\details{
  This function is a straightforward wrapper to \code{write.table}.
}
\value{
  None.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{write.table}}
}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
