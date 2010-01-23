str2time <- function(tt) {
  tt <- as.numeric(tt)
  hh <- as.integer(trunc(tt / 10000))
  tt <- tt %% 10000
  mm <- as.integer(trunc(tt / 100))
  ss <- tt %% 100
  list(txt = sprintf("%02i:%02i:%08.5f", hh, mm, ss),
       time = 3600*hh + 60*mm + ss)
}
  
str2date <- function(dd)
  format(as.Date(dd, "%Y%m%d"), "%d %b %Y")
