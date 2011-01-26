% Part of the oro.dicom package for R
% Distributed under the BSD license: see oro.dicom/COPYING
%
% $Id: $

\name{dec2base}
\alias{dec2base}
\alias{dec2hex}
\title{Convert Decimal to Base N Number in String}
\description{
  This function converts the nonnegative integer to the specified base.
}
\usage{
dec2base(n, base, len=0)
dec2hex(n, len=0)
}
\arguments{
  \item{n}{Non-negative integer.}
  \item{base}{Number between 2 and 36.}
  \item{len}{Length of the character string.}
}
\value{
  The returned argument is a string.
}
\details{
  This function converts the nonnegative integer \code{n} to the
  specified base, where \code{n} must be a nonnegative integer smaller
  than \eqn{2^52}, \code{base} must be an integer between 2 and 36 and
  \code{len} suggests the length of the character string.
}
%\references{}
%\seealso{}
\examples{
(x <- dec2base(23, 2))
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
