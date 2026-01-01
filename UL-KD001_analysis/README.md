# UL-KD001 ANALYSIS SUBDIRECTORY DESCRIPTION

This subdirectory houses all relevant files for analysing the preprocessed 
UL-KD001 data. 

## Files Overview

- rang_graph_analysis.rds: R-specific file containing the graph of the environment 
(packages and versions) created by the 'rang' package. Currently, no renv.lock
file can be created given the used packages, so this .rds file is a good 
compromise. 

```r

# The file can be read as:
graph <- readRDS('rang_graph.rds')
  
# Additionally, the specific packages can be installed as (not tested):
pkgs <- paste0(graph$pkgs$package, '@', graph$pkgs$version)
pak::pkg_install(pkgs)

# In the future, it may be possible to run the following for more convenience:
rang::export_renv(graph)
rang::dockerize(graph)

# Currently, the above fails. 

# The WRS package doesn't come from CRAN but can be downloaded from: https://github.com/nicebread/WRS. 

```

- analysis_helpers.R: R script with helper functions for the present analyses but
also in general. It is loaded by UL-KD001_analysis.R and UL-KD001_report.R. 

- UL-KD001_analysis.R: main R script for analysing the preprocessed data data. 
It contains all analyses that were performed. 

- UL-KD001_report.R: R script that was used to generate tables for the supplementary
materials




