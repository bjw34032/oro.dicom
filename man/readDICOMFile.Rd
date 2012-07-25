% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{readDICOMFile}
\alias{readDICOMFile}
\alias{dicomInfo}
\title{Read Single DICOM File}
\description{
  All information, both header and image, is read into a list structure
  from a DICOM file.
}
\usage{
readDICOMFile(fname, endian = "little", flipud = TRUE, skip128 = TRUE,
              DICM = TRUE, skipSequence = TRUE, pixelData = TRUE,
              warn = -1, debug = FALSE)
}
\arguments{
  \item{fname}{is the file name of the DICOM image (with suffix).}
  \item{endian}{is the endian-ness of the file (default is \code{"little"}).}
  \item{flipud}{is a logical variable for vertical flipping of the image
    (default is \code{TRUE}).}
  \item{skip128}{is a logical variable to skip the first 128 bytes in
    the file (default = \code{TRUE}).}
  \item{DICM}{is a logical variable to read in the four bytes "DICM"
    (default = \code{TRUE}).}
  \item{skipSequence}{is a logical variable to skip all content
    contained in SequenceItem tags (default = \code{TRUE}).}
  \item{pixelData}{is a logical variable (default = \code{TRUE}) on
    whether or not the PixelData should be read from the DICOM files.}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{options} for more details.}
  \item{debug}{is a logical variable (default = \code{FALSE}) that
    regulates to display of intermediate processing steps.}
  %\item{rowBuffer}{is the maximum number of DICOM header fields allowed
  %  to be read in (adjust when necessary).}
}
\details{
  A \code{while} loop is used to traverse the unknown number of DICOM
  header fields contained in a single file.  Information contained in
  \dQuote{sequences} may be included/excluded according to the logical
  variable \code{skipSequence} (default = \code{TRUE}).
}
\value{
  A list containing two elements:
  \describe{
    \item{hdr}{all DICOM header fields (with or without
      \dQuote{sequence} information).}
    \item{img}{the \sQuote{image} information.}
  }
}
\references{
  Whitcher, B., V. J. Schmid and A. Thornton (2011).  Working with the
  DICOM and NIfTI Data Standards in R, \emph{Journal of Statistical
    Software}, \bold{44} (6), 1--28.
  \url{http://www.jstatsoft.org/v44/i06}

  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}\cr
  \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
}
\seealso{
  \code{\link{readDICOM}}
}
\examples{
x <- readDICOMFile(system.file("dcm/Abdo.dcm", package="oro.dicom"))
graphics::image(t(x$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
                main="Abdo.dcm")

x <- readDICOMFile(system.file("dcm/Spine1.dcm", package="oro.dicom"))
graphics::image(t(x$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
                main="Spine1.dcm")

\dontrun{ # not yet 2012.07.24;  still debugging
x.jpeg <- readDICOMFile(system.file('jepg/cervicalKyphosis.dcm',
                        package='oro.dicom'))
}
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{file}
