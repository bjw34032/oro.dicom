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
  These functions are only used internally, via \code{dicomInfo}, to
  handle the various types of DICOM header fields.
}
\usage{
floating.header(VR, implicit, fid, endian)
other.header(fid, implicit, endian)
pixeldata.header(hdr, implicit, fid, endian)
signed.header(VR, implicit, fid, endian)
sequence.header(group, element, fid, implicit, endian, skipSQ, SQ, EOS)
unknown.header(VR, implicit, fid, endian, skipSQ)
unsigned.header(VR, implicit, fid, endian)
}
\arguments{
  \item{VR}{is the }
  \item{implicit}{is a logical value that denotes if the VR is implicit
    or not.}
  \item{fid}{is the \code{file} object.}
  \item{endian}{is a character string denoting the \dQuote{endianness}
    of the data, where possible values are \code{big}, \code{little} and
  \code{swap}.}
  \item{hdr}{is the DICOM header information provided in a \code{list}
    structure.}
  \item{group}{is the first hexadecimal number associated with the DICOM
    header element.}
  \item{element}{is the second hexadecimal number associated with the
    DICOM header element.}
  \item{skipSQ}{is a logical value (SQ = seqence).}
  \item{SQ}{is a character string (default = \code{NULL}).}
  \item{EOS}{is a logical value (EOS = end of sequence).}
}
%\details{}
\value{
  List with two elements:
  \item{length}{is the length of the \code{value} in bytes.}
  \item{value}{is the value of the DICOM header field.}
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
    \code{\link{readDICOM}}, \code{\link{readDICOMFile}}
}
% \examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
%\keyword{file}
