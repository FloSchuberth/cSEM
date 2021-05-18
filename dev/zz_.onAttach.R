.onAttach <- function(libname, pkgname) {
  
  version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                      fields = "Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage("Caution: ", pkgname, " is under active developement! Breaking changes may occur in the future.")
  packageStartupMessage("Please report improvements and bugs to: manuel-rademaker@outlook.de")
}