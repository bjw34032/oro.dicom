dec2base <- function(n, base, len=0) {
  symbols <- c(as.character(0:9), LETTERS)
  max.len <- trunc(log(max(n, 1)) / log(base)) + 1
  max.len <- max(max.len, len)
  ## determine digits for each number
  power <- rep(1, length(n)) * base^((max.len-1):0)
  n <- n * rep(1, max.len)
  digits <- floor((n %% (base*power)) / power)
  ## convert digits to symbols
  paste(symbols[digits+1], collapse="")
}

dec2hex <- function(n, len=0)
  dec2base(n, 16, len)

