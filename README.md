# VIS2026
Code related to a paper submitted to VIS 2026

This repository contains a zip file containing MATLAB code related to Section 3
of the submitted paper.

R code associated with ARA plots, analyzed in section 4, is available at CRAN:
https://cran.r-project.org/web/packages/aramappings/
The code currently implements the closed-form solution in Section 4.2, and the 
disambiguation approach in Section 4.3.1. The variants in Section 4.3.2 are not
implemented in the current version of the package. Instead, this repository
contains MATLAB file ara_infinity_norm_uniqueness.m that implements the approach
using the CVX toolbox.
