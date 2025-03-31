This is the code for reproducing the results, figures and tables appearing in:

"Scalable Fitting Methods for Multivariate Gaussian Additive Models with Covariate-dependent Covariance Matrices"

by Gioia et al. (2025).

To reproduce the results, one should set the working directory to that of the Main.R file and then run the code in the following files:

-   Main.R: this file runs all the simulations and fits all the models described in the paper. Running the whole file takes over a week on a 12-cores workstation. Memory usage is proportional to the number of cores on which the code is parallelised. On 12 cores, 200GBytes of RAM are needed.

-   Plot.R: this file uses the output of Main.R to produce all the plots and tables appearing in the paper and in the supplementary material.

Note that Main.R will install several R packages. Mostly of these are packages available on CRAN, with the exception of:

-   a modified version of the mgcv package. This is best installed on a local folder as done by Main.R.

-   the SCM package for fitting multivariate Gaussian GAMs. Main.R will install it from Github.

-   the electbook package which provides the GEFCom data. Main.R will install it from Github.
