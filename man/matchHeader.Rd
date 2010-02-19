% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{matchHeader}
\alias{matchHeader}
\title{Match String to DICOM Header Field}
\description{
  ...
}
\usage{
matchHeader(hdr, string)
}
\arguments{
  \item{hdr}{is the result from extracting information from a DICOM
    header field; e.g., using \code{extractHeader}.}
  \item{string}{is a character string to be matched with the DICOM
    header.}
}
\details{
  ...
}
\value{
  A logical vector of length \code{length(hdr)}.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{extractHeader}}
}
\examples{
x <- dicomInfo(system.file("dcm/Abdo.dcm", package="oro.dicom"))
modality <- extractHeader(x$hdr, "Modality", numeric=FALSE)
matchHeader(modality, "mr") # case insensitive by default
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{file}
