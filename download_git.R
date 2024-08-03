# Script for downloading packages from GitHub (mostly ChatGPT).
# This is required because renv can't handle GitHub installs...

packages <- list("briatte/ggnet" = "da9a7cf", "thomasp85/scico" = "v1.5.0")

# Function to install a specific version of a package from GitHub
install_github_version <- function(repo, ref) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org/")
  }
  devtools::install_github(repo, ref = ref, upgrade = "never")
}

# Install the specified versions of the packages from GitHub
for (pkg in names(packages)) {
  ref <- packages[[pkg]]
  pkg_name <- basename(pkg)
  
  message(sprintf("Installing %s version %s from GitHub repository %s", pkg_name, ref, pkg))
  install_github_version(pkg, ref)
}

message("All specified packages are installed with the correct versions from GitHub.")
