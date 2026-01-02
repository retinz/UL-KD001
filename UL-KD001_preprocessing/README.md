# UL-KD001 PREPROCESSING SUBDIRECTORY DESCRIPTION

This subdirectory houses all relevant files for preprocessing the UL-KD001 data. 

## Files Overview

- renv.lock: file containing the graph of the environment 
(packages and versions) created by the 'rang' package.

```r

# To restore the environment, run the following code. 
# The renv.lock file is assumed to be in the current working directory. 
install.packages('renv')
renv::restore()

```

- questionnaires_prep_pseudo.R: R script for preprocessing the questionnaire data. 
It is the newer version which uses the already pseudo-ananonymised data. 

- taskswitch_prep_pseudo.R: R script for preprocessing the cognitive task data. 
It is the newer version which uses the already pseudo-ananonymised data. 


