.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Welcome to GlycoDiveR")

  theme_set(theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                               panel.grid = element_blank()))
}
