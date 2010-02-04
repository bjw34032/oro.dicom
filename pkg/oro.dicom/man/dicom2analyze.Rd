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
dicom2analyze(img, hdr, descrip="SeriesDescription", ...)
}
\arguments{
  \item{img}{Multi-dimensional data volume}
  \item{hdr}{DICOM header}
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
  \code{\link{dicom2nifti}}, \code{\link[dcemriS4]{anlz}}
}
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{file}
