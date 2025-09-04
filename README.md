This code in R accompanies the manuscript: Jansma, Landi and Hui, Bayesian inference captures metabolite-bacteria interactions in a microbial community (DOI....).
The script "ODE_Bayesian_inference_Metabolite_bacteria_Interactions.R" is used to simulate bacterial communities and to estimate interaction coefficients using Rstan.
The script "ODE_Bayesian_inference_Metabolite_bacteria_Interactions_Diagnostics_Analysis.R" is used for all analyses and diagnostics and includes figure construction.
Comments and explanations are included in the scripts.

The data generated can be found here: DOI:10.5281/zenodo.17053090

The R libraries needed to execute both script fully are:
library(deSolve)
library(rstan)
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

--------------------------------------------------------------
The library Rstan is dependend on the C++ compiler Rtools. Parameter estimation in the manuscript is performed using Rtools 44, which can be installed via: https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html).

**If you have a version of R starting with R4.3, you need to install Rtools 43 (https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html).
If you have a version of R starting with R4.2, you need to install Rtools 42 (https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html).**

After the installation of Rtools, Rstan can be installed. First any old versions of Rstan need to be removed. For this, go to R and run:
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")

After that, still in R, run:
install.packages("rstan", repos = https://cloud.r-project.org/, dependencies = TRUE)

To check the installation, run:
example(stan_model, package = "rstan", run.dontrun = TRUE)

**It should not give any errors.**

Rstan can be accessed via the R package cmdstan. To install the package, still in R, run:
devtools::install_github("stan-dev/cmdstanr")

Then run:
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

This should give the message:
**The C++ toolchain required for CmdStan is setup properly!**

After that run:
install_cmdstan()

Sometimes other dependencies need to be installed before being able to use Rstan, with the most likely packages being:
install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))

---------------------------------------------------------------

All simulated communities as well as all data regarding parameter estimation can be found on this github page
All simulated communities can be recreated using set.seed(**42**)
