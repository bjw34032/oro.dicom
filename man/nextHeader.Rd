% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{nextHeader}
\alias{nextHeader}
\title{Check String Against DICOM Header Field to Produce Error Message
  or NEXT}
\description{
  ...
}
\usage{
nextHeader(dcm, string, reference, str.warning, htmlfile = NULL, heading = 3, numeric = FALSE)
}
\arguments{
  \item{dcm}{is the DICOM list structure.}
  \item{string}{is a character string to be matched with the DICOM
    header.}
  \item{reference}{is the scalar/vector of character strings to check
    against the DICOM header output.}
  \item{str.warning}{is a text string for the warning.}
  \item{htmlfile}{is the \pkg{hwriter} object for the HTML file (default
    = \code{NULL}.}
  \item{heading}{is the HTML tag <H?> (default = \code{3}).}
  \item{numeric}{is the argument to be passed to \code{matchHeader}.}
}
\details{
  ...
}
\value{
  An expression to be evaluated.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{extractHeader}}, \code{\link{matchHeader}}
}
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{file}
