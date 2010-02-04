% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{dicomSeparate}
\alias{dicomSeparate}
\title{Read All DICOM Files in a Directory}
\description{
  All DICOM files are imported and a text file summarizing their content
  recorded.
}
\usage{
dicomSeparate(path, debug=TRUE, ...)
}
\arguments{
  \item{path}{Path name to the DICOM directory.}
  \item{debug}{Flag to provide limited progress information.}
  \item{...}{Arguments to be passed to \code{dicomInfo}.}
}
\details{
  ...
}
\value{
  A list structure with two major components:
  \item{img}{All images associated with the DICOM directory.}
  \item{hdr}{All header files associated with the DICOM directory.}
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{dicomInfo}}
}
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{file}
