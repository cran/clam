## ----eval=FALSE---------------------------------------------------------------
# install.packages('clam')

## -----------------------------------------------------------------------------
library(clam)

## ----eval=FALSE---------------------------------------------------------------
# clam()

## ----echo=FALSE, fig.width=5, fig.height=5------------------------------------
base_temp_dir <- tempdir()
clam_dir <- file.path(base_temp_dir, "clam_runs")
dir.create(clam_dir, recursive = TRUE, showWarnings = FALSE)
clam_dir <- clam:::assign_coredir(clam_dir, core="Example", ask=FALSE)
clam(, ask=FALSE, coredir=clam_dir)

## ----echo=FALSE, eval=FALSE---------------------------------------------------
# clam(type=4, outliers=6, hiatus=470, dmax=800, slump=c(220, 250))

## ----fig.width=5, fig.height=5------------------------------------------------
clam(ask=FALSE, type=4, outliers=6, hiatus=470, dmax=800, slump=c(220, 250), coredir=clam_dir)

