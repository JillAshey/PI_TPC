# Load base installer helpers
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Named list: package = source
# "cran" = install from CRAN, "bioc" = Bioconductor, "gh" = GitHub repo
packages <- list(
  dplyr       = "cran",
  purrr       = "cran",
  minpack.lm  = "cran",
  tidyverse   = "cran",
  nls.multstart = "cran",
  broom       = "cran",
  LoLinR      = "gh:colin-olito/LoLinR",  # GitHub source
  readr       = "cran",
  lubridate   = "cran",
  fuzzyjoin   = "cran",
  stringr     = "cran",
  future      = "cran",
  furrr       = "cran"
)

install_and_load <- function(pkg, source){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " from ", source)
    if (source == "cran") {
      install.packages(pkg)
    } else if (source == "bioc") {
      BiocManager::install(pkg)
    } else if (startsWith(source, "gh:")) {
      repo <- sub("gh:", "", source)
      remotes::install_github(repo)
    }
  }
  library(pkg, character.only = TRUE)
}

# Loop over all
invisible(mapply(install_and_load, names(packages), packages))
