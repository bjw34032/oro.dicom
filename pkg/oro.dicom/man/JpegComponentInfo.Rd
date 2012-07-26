% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{JpegComponentInfo}
\alias{JpegComponentInfo}
\title{
  Constructor for JPEG Component Information
}
\description{
  The lossless JPEG code in C from Cornell includes "typedef struct
  JpegComponentInfo" consisting of 4 variables of class "short".  In R,
  these are represented as raw vectors of length 2.
}
\usage{
JpegComponentInfo(componentID, componentIndex,
                  hSampFactor, vSampFactor)
}
\arguments{
  \item{componentID}{
    identifier for this component (0..255)
  }
  \item{componentIndex}{
    index in SOF or cPter->compInfo[] (see "jpeg.h" in the LJPG C
    code).
  }
  \item{hSampFactor,vSampFactor}{
    horizontal and verticalsampling factors
  }
}
\details{
  Constructor function for an object of class \code{JpegComponentInfo}.
  For details, see Huang and Smith (1994).
}
\value{
  A list with 4 components matching the arguments of this constructor
  function, each being a raw vector of length 2.  The input arguments
  must either be of class raw or capable of coersion via
  \code{\link{as.raw}}.
}
\references{
  Kongji Huang and Brian Smith (1994) "Lossless JPEG Codec", v. 1.0
  \url{ftp://ftp.cs.cornell.edu/pub/multimed/ljpg.tar.Z}
}
\seealso{
  \code{\link{as.raw}}
}
\examples{
jpegCI <- JpegComponentInfo(1, 2, 3, 4)

jpegCI. <- list(componentID=as.raw(1), componentIndex=as.raw(2),
                hSampFactor=as.raw(3), vSampFactor=as.raw(4))
class(jpegCI.) <- 'JpegComponentInfo'

\dontshow{stopifnot(}
all.equal(jpegCI, jpegCI.)
\dontshow{)}

}
\author{Spencer Graves}
\keyword{classes}
