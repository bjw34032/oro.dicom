% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{header2matrix}
\alias{header2matrix}
\title{Converts DICOM Header Field to a Matrix}
\description{
  Converts a vector of DICOM header information, assuming there are
  multiple entries per element of the vector, into a matrix.
}
\usage{
header2matrix(hdr, ncol, sep = " ", byrow = TRUE)
}
\arguments{
  \item{hdr}{is the result from extracting information from a DICOM
    header field; e.g., using \code{extractHeader}.}
  \item{ncol}{is the number of columns.}
  \item{sep}{is the character string required to split entries in the
    header field.}
  \item{byrow}{is a logical variable (default = \code{TRUE}) telling the
    routine to populate the matrix by rows then columns.}
}
%\details{}
\value{
  Matrix with \code{length(hdr)} rows and \code{ncol} columns.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{extractHeader}}, \code{\link{matrix}}
}
\examples{
x <- dicomInfo(system.file("dcm/Abdo.dcm", package="oro.dicom"))
pixelSpacing <- extractHeader(x$hdr, "PixelSpacing", numeric=FALSE)
pSmat <- header2matrix(pixelSpacing, ncol=2)
IOP <- extractHeader(x$hdr, "ImageOrientationPatient", numeric=FALSE)
IOPmat <- header2matrix(IOP, ncol=6)
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
