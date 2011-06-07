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
dicom2nifti(dcm, datatype=4, units=c("mm","sec"), rescale=FALSE,
            reslice=TRUE, sform=TRUE, DIM=3,
            descrip="SeriesDescription", aux.file=NULL, ...)
}
\arguments{
  \item{dcm}{DICOM object containing both header and image information.}
  \item{datatype}{is an integer that denotes the type of data contained
    in each voxel.  See \code{convert.datatype} or the NIfTI
    documentation for more details.}
  \item{units}{Spatial and temporal units for \code{xyzt}}
  \item{rescale}{Should slope and intercept parameters be extracted from
    the DICOM headers and saved?}
  \item{reslice}{Logical variable (default = \code{TRUE}) indicating if
    the data volume should be resliced.}
  \item{sform}{Logical variable (default = \code{TRUE}) indicating if
    the 3D image orientation should be used.}
  \item{DIM}{The dimension of the array to be used (default = 3D).}
  \item{descrip}{DICOM header field(s) to be included in the
    \code{descrip} slot for the \code{nifti} class object.}
  \item{aux.file}{Character string to be included in the \code{aux_file}
    slot for the \code{nifti} class object.}
  \item{\dots}{Arguments to be passed to \code{nifti}}
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
  \code{\link[oro.nifti]{convert.datatype}}, \code{\link{dicom2analyze}},
  \code{\link[oro.nifti]{nifti}}
}
\examples{
\dontrun{
dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"))
require("oro.nifti")
dcmNifti <- dicom2nifti(dcmList, datatype=4, mode="integer")
image(dcmNifti)
orthographic(dcmNifti)
}
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
