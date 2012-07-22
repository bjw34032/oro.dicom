subRaw <- function(pattern, replacement, x, offset=1L,
                   ignore.case=FALSE, fixed=TRUE, all=FALSE){
##
## find
##
  where <- grepRaw(pattern, x=x, offset=offset,
                   ignore.case=ignore.case, fixed=fixed)
  if((length(where)>0) && is.numeric(where)){
##
## what
##
    what <- grepRaw(pattern, x=x, offset=offset,
          ignore.case=ignore.case, value=TRUE, fixed=fixed)
    nwhat <- length(what)
##
## replace
##
    nx <- length(x)
    x <- c(x[seq(1, length=where-1)], replacement,
           x[seq(where+nwhat, length=nx-where-nwhat+1)])
    if(all) x <- subRaw(pattern, replacement, x,
                        ignore.case=ignore.case, fixed=fixed, all=all)
  }
##
## return
##
  x
}

