% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{dicom2analyze}
\alias{dicom2analyze}
\title{Convert DICOM Header to Analyze}
\description{
  A subset of header information from DICOM is placed into Analyze 7.5
  format.
}
\usage{
dicom2analyze(dcm, reslice = TRUE, descrip = "SeriesDescription",
              ...)
}
\arguments{
  \item{dcm}{DICOM object containing both header and image information.}
  \item{reslice}{Logical variable (default = \code{TRUE}) indicating if
    the data volume should be resliced.}
  \item{descrip}{DICOM header field(s) to be included in the
    \code{descrip}}
  \item{...}{Arguments to be passed to \code{anlz}}
}
\value{
  An object of class \code{anlz}.
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
  \code{\link{dicom2nifti}}, \code{\link[oro.nifti]{anlz}}
}
\examples{
\dontrun{
dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"))
require("oro.nifti")
dcmAnlz <- dicom2analyze(dcmList, datatype=4, mode="integer")
image(dcmAnlz)
orthographic(dcmAnlz)
}
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
