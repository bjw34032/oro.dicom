\name{dicomInfo}
\alias{dicomInfo}
\alias{dicom.fields}
\alias{dicom.VRfields}
\title{Read Single DICOM File}
\description{
  All information, both header and image, is read into a list structure
  from a DICOM file.
}
\usage{
dicomInfo(fname, endian = "little", flipud = TRUE, skip128 = TRUE,
          DICM = TRUE, warn=-1)
}
\arguments{
  \item{fname}{is the file name of the DICOM image (with suffix).}
  \item{endian}{is the endian-ness of the file (default is \code{"little"}).}
  \item{flipud}{is a logical variable for vertical flipping of the image
    (default is \code{TRUE}).}
  \item{skip128}{is a logical variable to skip the first 128 bytes in
    the file (default = \code{TRUE}).}
  \item{DICM}{is a logical variable to read in the four bytes "DICM"
    (default = \code{TRUE}).}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{options} for more details.}
}
\details{
}
\value{
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}\cr
  \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
}
\seealso{
  \code{\link{dicomSeparate}}
}
\examples{
x <- dicomInfo(system.file("dcm/Abdo.dcm", package="DICOM")[1])
image(x$img, col=grey(0:64/64), axes=FALSE, xlab="", ylab="")

x <- dicomInfo(system.file("dcm/Spine1.dcm", package="DICOM")[1])
image(x$img, col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
