% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{DICOM Headers}
\alias{floating.header}
\alias{other.header}
\alias{pixeldata.header}
\alias{sequence.header}
\alias{signed.header}
\alias{unknown.header}
\alias{unsigned.header}
\title{Internal Functions for DICOM Headers}
\description{
  These functions are only used in \code{dicomInfo} to handle the
  various types of header fields.
}
\usage{
floating.header(VR, implicit, fid, endian)
other.header(fid, implicit, endian)
pixeldata.header(hdr, implicit, fid, endian)
signed.header(VR, implicit, fid, endian)
sequence.header(group, element, fid, implicit, endian, skipSQ, SQ, EOS)
unknown.header(VR, implicit, fid, endian)
unsigned.header(VR, implicit, fid, endian)
}
\arguments{
  \item{VR}{...}
  \item{implicit}{...}
  \item{fid}{...}
  \item{endian}{...}
  \item{hdr}{...}
  \item{group}{...}
  \item{element}{...}
  \item{skipSQ}{...}
  \item{SQ}{...}
  \item{EOS}{...}
}
%\details{}
\value{
  List with two elements:
  \item{length}{...}
  \item{value}{...}
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
    \code{\link{dicomInfo}}, \code{\link{dicomSeparate}}
}
% \examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
%\keyword{file}
