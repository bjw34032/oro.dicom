.onAttach <- function (lib, pkg) {
  cat("\n", pkg,": Rigorous - DICOM Input / Output (version = ",
      as.character(sessionInfo()$otherPkgs$oro.dicom["Version"]), ")\n",
      sep="", fill=TRUE)
}

#.onLoad <- function (lib, pkg) {
#  data("dicom.dic")
#}

