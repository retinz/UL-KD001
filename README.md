# UL-KD001 STUDY

This project contains files used for preprocessing and analysing data of 
a nonrandomised controlled trial, comparing a ketogenic diet with a clean diet
in terms of mood, sleep, ADHD symptoms, and cognitive performance. 

This study was carried out by the LBI research group at Leiden University. 

This directory is the root of the project. This is also marked by the '.here'
file, which indicates to R (package 'here') that this is the root (main) directory. 

## Files Overview

- .gitignore: R file specifying which files should not be synced by git and GitHub

- UL-KD001_analysis: houses all relevant files for analysing UL-KD001 data; 
assumes preprocessed data

- UL-KD001_preprocessing: contains relevant files for preprocessing pseudo-anonymised
UL-KD001 raw data

- rang_graph.rds: R-specific file containing all the packages and their versions 
(graph of the environment) used by the scripts in this repository. The graph was 
created by the 'rang' package. Currently, no renv.lock file can be created given 
the packages used, so this .rds file is a good compromise. The section below outlines
how the graph file can be loaded and the packages installed.

```r

# Install pak if it isn't already installed
if (!requireNamespace('pak', quietly = TRUE)) {
  install.packages('pak')
}

# Read graph file
graph <- readRDS('rang_graph.rds')

# Custom recursive function to dig out packages and versions
extract_pkgs <- function(node) {
  # If we hit a dataframe containing package info, paste them together
  if (is.data.frame(node) && 'x' %in% names(node) && 'x_version' %in% names(node)) {
    return(paste0(node$x, '@', node$x_version))
  } else if (is.list(node)) {
    # If it is a list, keep digging deeper
    return(unlist(lapply(node, extract_pkgs)))
  }
}

# Extract everything and keep only the unique package@version strings
pkgs <- unique(extract_pkgs(graph))

# Install all packages
pak::pkg_install(pkgs) # This particular step wasn't tested

# The WRS package doesn't come from CRAN but can be downloaded from: https://github.com/nicebread/WRS. 

```

The graph file was created in the following way. 

```r

library(rang)

# 1. Define the files to scan
files_to_scan <- c(
  'UL-KD001_preprocessing/questionnaires_prep_pseudo.R', 
  'UL-KD001_preprocessing/taskswitch_prep_pseudo.R',
  'UL-KD001_analysis/UL-KD001_analysis_core.R',
  'UL-KD001_analysis/UL-KD001_analysis_extended.R',
  'UL-KD001_analysis/analysis_helpers.R',
  'UL-KD001_analysis/UL-KD001_report.R')

# 2. Scan the scripts to catch all used packages
scanned_pkgs <- unique(renv::dependencies(files_to_scan, quiet = TRUE)$Package)

# 3. Get the exact local versions for the scanned packages
scanned_vers <- vapply(scanned_pkgs, function(pkg) as.character(packageVersion(pkg)), character(1))
all_pkgs <- paste0(scanned_pkgs, '@', scanned_vers)

# 4. Drop any packages you don't want (including base R packages)
base_pkgs <- rownames(installed.packages(priority = 'base'))
drop <- c('WRS', base_pkgs)

# Create a regex pattern to match exactly 'PackageName@'
pattern <- paste0('^(', paste(drop, collapse = '|'), ')@')
keep <- !grepl(pattern, all_pkgs)
final_pkgs <- all_pkgs[keep]

# 5. Resolve the final graph
graph <- rang::resolve(
    pkgs = final_pkgs,
    snapshot_date = '2026-04-17',
    query_sysreqs = FALSE
)

# 6. Fix R version if needed
graph$r_version <- as.character(getRversion())

# 7. Generate graph RDS
saveRDS(graph, file='rang_graph.rds')

```

