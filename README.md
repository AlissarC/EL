This repository contains the data and R code used to perform analyses for "Soil nitrogen supply exerts largest influence on leaf nitrogen in environments with the greatest leaf nitrogen demand" by Cheaib et al. (2024). The manuscript is currently in revision at _Ecology Letters_. Folder contents are as follows:

`data` = An initial file `data/traits.csv` compiles all climate, leaf traits, and biomass data. Relevant calculations for the manuscript are calculated in `scripts/leaf_processing.R`, which gives a .csv file (`data/traits_nofence_chi_sub.csv`) that is then used for all analyses. A metadata file for `data/traits.csv` is also included in this folder 

`scripts` = contains R scripts used for leaf trait calculations (`scripts/leaf_processing.R`), data analysis (`scripts/leafN_analysis.R`, `scripts/delta_N_AGB_analysis.R`), and supplemental plots (`scripts/diagram.R`)

You can install plotbiomes from GitHub by visiting https://github.com/valentinitnelav/plotbiomes/blob/develop/README.md, and as follow:
using devtools package:
if (!require(devtools)) install.packages("devtools")
devtools::install_github("valentinitnelav/plotbiomes")
using remotes package:
if (!require(remotes)) install.packages("remotes")
remotes::install_github("valentinitnelav/plotbiomes")


R package information: R version 4.3.2, R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for
  Statistical Computing, Vienna, Austria. <https://www.R-project.org/>)