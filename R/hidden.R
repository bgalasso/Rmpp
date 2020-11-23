.onAttach <- function(libname, pkgname) {
  packageStartupMessage("## ========================================= ##")
  RFver <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                    fields = "Version")
  packageStartupMessage(paste
                       ("##", pkgname, RFver, "                                 ##"))
  packageStartupMessage("## ----------------------------------------- ##")
  packageStartupMessage("##  Copyright (C) 2020                       ##")
  packageStartupMessage("##  B. Galasso, Y. Zemel and M. de Carvalho  ##")
  packageStartupMessage("## ========================================= ##")
  packageStartupMessage("")
}
