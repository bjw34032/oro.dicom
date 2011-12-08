.onAttach <- function(lib, pkg) {
  txt <- paste("\n", pkg,": Rigorous - DICOM Input / Output (version = ",
               as.character(sessionInfo()$otherPkgs$oro.dicom["Version"]),
               ")\n", sep="", fill=TRUE)
  packageStartupMessage(txt)
}

##.onLoad <- function(lib, pkg) {
##  data("dicom.dic")
##}

