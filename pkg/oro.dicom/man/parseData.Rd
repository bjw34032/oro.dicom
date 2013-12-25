% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{Parse DICOM Data}
\alias{parsePixelData}
\alias{parseSpectroscopyData}
\title{Parse DICOM Pixel or Spectroscopy Data}
\description{
  These subroutines process the information contained after the DICOM header and
  process this information into an image (2D or 3D) or complex-valued vector.
}
\usage{
parsePixelData(rawString, hdr, endian = "little", flipupdown = TRUE)
parseSpectroscopyData(rawString, hdr, endian = "little")
}
\arguments{
  \item{rawString}{is a vector of \code{raw} values taken directly from
    the DICOM file.}
  \item{hdr}{is the list object of DICOM header information.}
  \item{endian}{is the endian-ness of the file (default is \code{"little"}).}
  \item{flipupdown}{is a logical variable for vertical flipping of the image
    (default is \code{TRUE}).}
}
\details{
  A \code{while} loop is used to traverse the unknown number of DICOM
  header fields contained in a single file.  Information contained in
  \dQuote{sequences} may be included/excluded according to the logical
  variable \code{skipSequence} (default = \code{TRUE}).

  A resursive implementation of the code breaks the DICOM file into 
  segments and calls itself to parse each segment.  
}
\value{
  A list containing two elements:
  \describe{
    \item{hdr}{all DICOM header fields (with or without
      \dQuote{sequence} information).}
    \item{img}{the \sQuote{image} information.}
  }
}
\source{
  See references.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}\cr
  \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
}
\seealso{
  \code{\link{parseDICOMHeader}}, \code{\link{readDICOMFile}}.
}
%\examples{}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
\keyword{file}
