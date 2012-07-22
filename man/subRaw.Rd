\name{subRaw}
\alias{subRaw}
\title{
  Pattern Replacement for Raw Vectors
}
\description{
  \code{subRaw} searchs for substring \code{pattern} within a raw vector
  \code{x} and replaces it if found with \code{replacement}.
}
\usage{
subRaw(pattern, replacement, x, offset=1L, ignore.case=FALSE,
       fixed=TRUE, all=FALSE)
}
\arguments{
  \item{pattern}{
    raw vector containing the pattern to be matched and replaced by
    \code{replacement}.

    This is passed to \code{\link{grepRaw}}, which supports the use of
    regular expressions in \code{pattern} if \code{fixed} = FALSE.
  }
  \item{replacement}{
    a replacement for matched pattern.
  }
  \item{x}{
    a raw vector where matches are sought, or an object which can be
    coerced by \code{charToRaw} to a raw vector.
  }
  \item{offset}{
    An integer passed to \code{grepRaw} specifying the offset from which
    the search should start.  Must be positive.  The beginning of line
    is defined to be at that offset so \code{"^"} will match there.  If
    all, offset only applies to the first match.
  }
  \item{ignore.case}{
    logical passed to \code{grepRaw}:  if 'FALSE', the pattern matching
    is _case sensitive_ and if 'TRUE', case is ignored during matching.
  }
  \item{fixed}{
    logical passed to \code{grepRaw}:  If 'TRUE', 'pattern' is a pattern
    to be matched as is.

    NOTE:  The default is \code{TRUE} here, even though the default for
    the \code{fixed} argument in \code{grepRaw} is FALSE.
  }
  \item{all}{
    logical.  If 'TRUE' all matches are returned, otherwise just the
    first one.
  }
}
\details{
  1.  found <- grepRaw(pattern, x=x, offset=offset,
  ignore.case=ignore.case, fixed=fixed)

  2.  if(length(found)>0) find the match and replace x with
  c(x[1:(start-1)], replace, c[(end+1):length(x)]).  Also, if(all)
  recurse.

  3.  return(x)
}
\value{
  a raw vector matching \code{x} except with any subsequence of bytes
  matching \code{pattern} replaced by \code{replacement}.
}
%\references{}
\seealso{
  \code{\link{array}}, \code{\link{readDICOM}},
  \code{\link{storage.mode}}
}
\author{Spencer Graves}
\examples{
##
## example
##
abcdef <- charToRaw('abcdef')
ab. <- subRaw(charToRaw('ab'), charToRaw('d'), abcdef)

# replace first two bytes
\dontshow{stopifnot(}
all.equal(ab., charToRaw('dcdef') )
\dontshow{)}

# replace lst two bytes
\dontshow{stopifnot(}
all.equal(subRaw(charToRaw('ef'), charToRaw('g'), abcdef),
          charToRaw('abcdg') )
\dontshow{)}

# replace nothing
\dontshow{stopifnot(}
all.equal(subRaw(charToRaw('noop'), charToRaw('c'), abcdef),
          abcdef )
\dontshow{)}

# replace all
ababab <- charToRaw('ababab')
ab3 <- subRaw(charToRaw('ab'), charToRaw('g'), ababab,
               all=TRUE)

\dontshow{stopifnot(}
all.equal(ab3, charToRaw('ggg') )
\dontshow{)}
}
%\keyword{}
