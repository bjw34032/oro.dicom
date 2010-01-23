dicom.VRfields <- function() {
  VR <- matrix(c("AE", "Application Entity", 16, 0,
                 "AS", "Age String", 4, 1, 
                 "AT", "Attribute Tag", 4, 1, 
                 "CS", "Code String", 16, 0, 
                 "DA", "Date", 8, 1, 
                 "DS", "Decimal String", 16, 0, 
                 "DT", "Date Time", 26, 0, 
                 "FL", "Floating Point Single", 4, 1, 
                 "FD", "Floating Point Double", 8, 1, 
                 "IS", "Integer String", 12, 0, 
                 "LO", "Long Strong", 64, 0, 
                 "LT", "Long Text", 10240, 0, 
                 "OB", "Other Byte String", 0, 0, 
                 "OW", "Other Word String", 0, 0, 
                 "PN", "Person Name", 64, 0, 
                 "SH", "Short String", 16, 0, 
                 "SL", "Signed Long", 4, 1, 
                 "SQ", "Sequence of Items", 0, 0, 
                 "SS", "Signed Short", 2, 1, 
                 "ST", "Short Text", 1024, 0, 
                 "TM", "Time", 16, 0, 
                 "UI", "Unique Identifier UID", 64, 0, 
                 "UL", "Unsigned Long", 4, 1, 
                 "UN", "Unknown", 0, 0, 
                 "US", "Unsigned Short", 2, 1, 
                 "UT", "Unlimited Text", 0, 0),
               ncol=4, byrow=TRUE)
  dimnames(VR) <- list(NULL, c("code", "name", "bytes", "fixed"))
  VR <- data.frame(VR, stringsAsFactors=FALSE)
  VR$bytes <- as.numeric(VR$bytes)
  VR$fixed <- as.numeric(VR$fixed)
  return(VR)
}

dicom.fields <- function() {
  data("dicom.fields.data")
  dicom.fields.data
}
