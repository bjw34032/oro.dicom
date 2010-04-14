% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{Create}
\alias{create3D}
\alias{create4D}
\title{Create Arrays from DICOM Headers/Images}
\description{
  A DICOM list structure is used to produce a multi-dimensional array
  representing a single acquisition of medical imaging data.
}
\usage{
create3D(dcm, mode = "integer", transpose = TRUE, pixelData = TRUE,
         mosaic = FALSE, mosaicXY = NULL)
create4D(dcm, mode = "integer", transpose = TRUE, pixelData = TRUE,
         mosaic = FALSE, mosaicXY = NULL, W = NULL)
}
\arguments{
  \item{dcm}{is the DICOM list structure (if \code{pixelData} =
    \code{TRUE}) or the DICOM header information (if \code{pixelData} =
    \code{FALSE}).}
%  \item{X,Y,Z}{are the dimensions of the array.}
%  \item{sliceLocation}{is a vector of slice locations from the DICOM header.}
%  \item{patientPosition}{is a character string denoting the patient
%    position in the scanner.}
  \item{mode}{is a valid character string for \code{storage.mode}.}
  \item{transpose}{is necessary to switch the definition of rows and
    columns from DICOM to R.}
  \item{pixelData}{is a logical variable (default = \code{TRUE}) that is
    associated with the DICOM image data being pre-loaded.}
  \item{W}{is the fourth dimension of the array.}
  \item{mosaic}{is a logical variable (default = \code{FALSE}) to denote
    storage of the data in Siemens \sQuote{Mosaic} format.}
  \item{mosaicXY}{is a vector of length two that provides the (x,y)
    dimensions of the individual images.  Default behavior is to use the
    AcquisitonMatrix to determine the (x,y) values.}
%  \item{\ldots}{\ldots}
}
\value{
  Multi-dimensional array of medical imaging data.
}
%\details{}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{array}}, \code{\link{dicomInfo}},
  \code{\link{storage.mode}}
}
\examples{
## pixelData = TRUE
## The DICOM image data are read from dicomSeparate()
dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"))
dcmImage <- create3D(dcmList)
image(dcmImage[,,1], col=grey(0:64/64), axes=FALSE, xlab="", ylab="")
imagePositionPatient <- attributes(dcmImage)$ipp
dSL <- abs(diff(imagePositionPatient[,3]))
plot(dSL, ylim=range(range(dSL) * 1.5, 0, 10), xlab="Image", ylab="mm",
     main="Difference in Slice Location")

## pixelData = FALSE
## The DICOM image data are read from create3D()
## This may save on memory for large batches of DICOM data
dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"),
                         pixelData=FALSE)
dcmImage <- create3D(dcmList, pixelData=FALSE)
image(dcmImage[,,1], col=grey(0:64/64), axes=FALSE, xlab="", ylab="")

## mosaic = TRUE
dcm <- dicomInfo(system.file("dcm/MR-sonata-3D-as-Tile.dcm",
                             package="oro.dicom"))
dcmImage <- create3D(dcm, mode="integer", mosaic=TRUE)
image(dcmImage[,,dim(dcmImage)[3]/2], col=grey(0:64/64), axes=FALSE,
      xlab="", ylab="")
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
