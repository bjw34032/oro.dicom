% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{swapDimension}
\alias{swapDimension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reslice Data Volume Using DICOM Header Fields
}
\description{
  The input data volume (assumed to be three-dimensional) is re-sliced
  so that each slice is in the axial plane.  Orientation is preserved so
  that orthographic viewing is standardized.
}
\usage{
swapDimension(img, dcm, digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{Multidimensional array (assumed to be three-dimensional
    only).}
  \item{dcm}{DICOM header/image object associated with the
    multidimensional array.}
  \item{digits}{Number of significant digits used in testing
    \code{unique}-ness of values in DICOM header fields.}
}
%\details{}
\value{
  Multidimensional array with (potentially) permuted dimensions because
  of the reslicing operation.  An additional attribute
  \dQuote{\code{pixdim}} is provided in order to facilitate conversion
  from DICOM to NIFTI/ANALYZE.
}
%\references{}
\author{
Brandon Whitcher \email{bjw34032@users.sourceforge.net}
}
%\note{}
\seealso{
  \code{\link{dicom2nifti}}, \code{\link{getOrientation}}
}
%\examples{}
\keyword{misc}
