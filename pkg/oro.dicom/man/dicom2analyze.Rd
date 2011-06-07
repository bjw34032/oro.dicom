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
dicom2analyze(dcm, datatype = 8, reslice = TRUE, DIM = 3,
              descrip = "SeriesDescription", ...)
}
\arguments{
  \item{dcm}{DICOM object containing both header and image information.}
  \item{datatype}{is an integer that denotes the type of data contained
    in each voxel.  See \code{convert.datatype.anlz} or the ANALYZE
    documentation for more details.}
  \item{reslice}{Logical variable (default = \code{TRUE}) indicating if
    the data volume should be resliced.}
  \item{DIM}{The dimension of the array to be used (default = 3D).}
  \item{descrip}{DICOM header field(s) to be included in the
    \code{descrip}}
  \item{\dots}{Arguments to be passed to \code{anlz}}
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
  \code{\link[oro.nifti]{convert.datatype.anlz}}, \code{\link{dicom2nifti}},
  \code{\link[oro.nifti]{anlz}}
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
