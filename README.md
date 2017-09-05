# FindFocals - An R package to detect focal gains and deletions in SNP array data


## Description

*FindFocals* is a simple package to detect small focal gains and deletions in 
SNP array data. It uses a naive approach to find clusters of SNPs above or 
below the expected range.

## Notes: 

  * The package is still in early development and may or may no work as expected.
  * This package has been developed for a specific need in our lab. It might not
    be developed further or maintained in any way.

## Install

The easiest way to install the package is using devtools

    library(devtools)
    install_github("bernatgel/FindFocals")

## How to use it

Load the data with `loadSNPData` and then process the results with `findFocals`.
The resulting object will contain the raw data and the focal gains and 
deletions. Two additional functions are available to plot the raw data 
(`plotSNPData`) and the complete results (`plotFocals`).





