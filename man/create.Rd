% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{Create}
\alias{create3D}
\title{Create Arrays from DICOM Headers/Images}
\description{
  A DICOM list structure is used to produce a multi-dimensional array
  representing a single acquisition of medical imaging data.
}
\usage{
create3D(dcm, mode = "double", transpose = TRUE, pixelData = TRUE,
         path = NULL)
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
  \item{path}{is the path to the directory holding the DICOM files (only
    necessary when \code{pixelData} = \code{FALSE}).}
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
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
