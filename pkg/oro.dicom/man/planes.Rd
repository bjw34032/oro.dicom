\name{Orthogonal Planes}
\alias{is.axial}
\alias{is.coronal}
\alias{is.sagittal}
\title{Orthogonal Planes}
\description{
  Functions to test the orientation for a single slice.
}
\usage{
is.axial(imageOrientationPatient, axial = c("L","R","A","P"))
is.coronal(imageOrientationPatient, coronal = c("L","R","H","F"))
is.sagittal(imageOrientationPatient, sagittal = c("A","P","H","F"))
}
\arguments{
  \item{imageOrientationPatient}{A vector of length six taken from the
    DICOM header field \dQuote{ImageOrientationPatient}.}
  \item{axial}{Characters that are valid in defining an \sQuote{axial}
    slice.}
  \item{coronal}{Characters that are valid in defining a
    \sQuote{coronal} slice.}
  \item{sagittal}{Characters that are valid in defining a
    \sQuote{sagittal} slice.}
}
%\details{}
\value{
  Logical value.
}
%\references{}
\seealso{
  \code{\link{getOrientation}}
}
\examples{
x <- dicomInfo(system.file("dcm/Abdo.dcm", package="oro.dicom"))
iop <-
  header2matrix(extractHeader(x$hdr, "ImageOrientationPatient", FALSE), 6)
is.axial(iop)
is.coronal(iop)
is.sagittal(iop)

x <- dicomInfo(system.file("dcm/Spine1.dcm", package="oro.dicom"))
iop <-
  header2matrix(extractHeader(x$hdr, "ImageOrientationPatient", FALSE), 6)
is.axial(iop)
is.coronal(iop)
is.sagittal(iop)
}
\author{
  Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>
}
\keyword{misc}
