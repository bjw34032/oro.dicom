% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{dicom2nifti}
\alias{dicom2nifti}
\title{Convert DICOM Header to NIfTI}
\description{
  A subset of header information from DICOM is placed into NIfTI-1
  format.
}
\usage{
dicom2nifti(img, hdr, units=c("mm","sec"), rescale=FALSE, 
            descrip="SeriesDescription", ...)
}
\arguments{
  \item{img}{Multi-dimensional data volume}
  \item{hdr}{DICOM header}
  \item{units}{Spatial and temporal units for \code{xyzt}}
  \item{rescale}{Should slope and intercept parameters be extracted from
    the DICOM headers and saved?}
  \item{descrip}{DICOM header field(s) to be included in the
    \code{descrip}}
  \item{...}{Arguments to be passed to \code{nifti}}
}
\value{
  An object of class \code{nifti}.
}
\details{
  See the references.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}

  NIfTI-1\cr
  \url{http://nifti.nimh.nih.gov/nifti-1}
}
\seealso{
  \code{\link{dicom2analyze}}, \code{\link[dcemriS4]{nifti}}
}
\examples{
\dontrun{
## pixelData = TRUE
dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"))
dcmImage <- create3D(dcmList, mode="integer")
## Accumulate voxel dimensions
pS <- header2matrix(extractHeader(dcmList$hdr, "PixelSpacing", FALSE), 2)[1,]
sT <- extractHeader(dcmList$hdr, "SliceThickness")[1]
pixdim <- c(0, pS, sT, rep(1,4))
require("oro.nifti")
dcmNifti <- nifti(dcmImage, datatype=4, pixdim=pixdim)
par(bg="black")
image(dcmNifti)
orthographic(dcmNifti)
}
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{file}
