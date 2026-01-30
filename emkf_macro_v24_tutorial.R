#############################################################################
#                                                                           #
# U.S. National Center for Health Statistics                                #
# Rmkf: R prototype of the SAS eMKF v2.4 macro                              #
# (using R 'nimble' version >= 1.3.0; see https://r-nimble.org)             #
#                                                                           #  
# Suggested citation:                                                       #
# Talih M, Patel P, Rossen LM. An R-NIMBLE implementation of the            #
# enhanced modified Kalman filter (eMKF) tool for small domain              #
# estimation (version 2.4 2026-01-30). National Center for                  # 
# Health Statistics. 2026. https://github.com/CDCgov/eMKF.                  #
#                                                                           #
#                                                                           #  
# Example usage and test cases                                              #
#                                                                           #
# Last modified: 30-Jan-2026                                                #
# Makram Talih, Ph.D.                                                       #
#                                                                           #
#############################################################################

###################################
# Set search path(s) for packages #
###################################

## Example .libPaths() statement to use in nonstandard environments, 
## e.g., the CSP (uncomment and edit as needed)

# .libPaths(c("R:/win-library/4.4", paste0("R:/Personal/", Sys.getenv("USERNAME"),"/4.4")))

##################################
# Load Rmkf package dependencies #
##################################

## Replace 'lib = NULL' with desired library path location, 
## e.g., 'lib = .libPaths()[2]', for nonstandard environments.

installRmkfPackages <- function(lib = NULL) {   
  
  ## Default library path if none specified 
  if (missing(lib) || is.null(lib)) { 
    lib <- .libPaths()[1] 
  }
  
  ## 'stargazer' package for R-side text and html table formatting and output
  if (!requireNamespace("stargazer", quietly=TRUE)) {
    install.packages("stargazer", lib=lib, quiet = TRUE)
  }
  
  ## 'msm' package for R-side implementation of the truncated normal
  if (!requireNamespace("msm", quietly=TRUE)) {
    install.packages("msm", lib=lib, quiet = TRUE)
  }
  
  ## 'coda' package for R-side MCMC post-processing functionality
  if (!requireNamespace("coda", quietly=TRUE)) {
    install.packages("coda", lib=lib, quiet = TRUE)
  }
  
  ## 'nimble' package (assumes C++ compiler and related tools are installed, 
  ## e.g., via Rtools; see https://r-nimble.org/download.html).
  ## Note that the Rtools and R versions need to match 
  ## (e.g., Rtools version 4.4 for R versions 4.4.x).
  
  if (!requireNamespace("nimble", quietly=TRUE)) {
    install.packages("nimble", lib=lib, quiet = TRUE)
  }
  if (packageVersion("nimble") < "1.3.0") {
    stop("Error: Please install 'nimble' version 1.3.0 or later.")
  }
  
} #installRmkfPackages

installRmkfPackages()

####################################
# User-specified working directory #
####################################

## If not specified, output tables will be written to the default working directory getwd()
WORKDIR = "D:\\eMKF"
setwd(WORKDIR)

## Change directory path for location of data and macro files as needed
DATADIR = "D:\\eMKF\\MKFdata" 
MACRDIR = "D:\\eMKF\\MKFmacro"

#############################################################
# User-specified temporary directory for DLLs and C++ files #
#############################################################

## Only required in nonstandard environments if *.dll files cannot be written to default tempdir().
## White space in directory name (e.g., "R:\\My Temp Directory") is not recommended.
TEMPDIR = "D:\\eMKFtemp" 
if (!dir.exists(TEMPDIR)) {
  dir.create(TEMPDIR)
}
Sys.setenv(TMPDIR=tools::file_path_as_absolute(TEMPDIR))

######################################################################################
# Example from NHIS -- used in the eMKF Evaluation Report (DOI: 10.15620/cdc/157497) #
######################################################################################

NHISdata <- read.csv(paste(DATADIR, "nhis2rev2.csv", sep="\\"), as.is=TRUE, fill=FALSE)

## Code in 'emkf_macro_v24.R' must be sourced prior to every Rmkf() call.
## Rmkf() removes custom samplers at the end of each call, so those must be recompiled.

source(paste(MACRDIR, "emkf_macro_v24.R", sep="\\"))

## Rmkf call mimics the SAS eMKF call, with most argument names and functionality carried over.
## See eMKF Guidance Report (DOI: 10.15620/cdc/157496) for a walk-through of SAS eMKF arguments.
## Below are some of the high-level Rmkf parameters that the user is likely to want to customize.
## For the full set of Rmkf call parameters, the user should refer to the README file.

Rmkf(out = "param",                # Default value is "param". This is the prefix to use for  
                                   # naming the output datasets produced by Rmkf().
     data = "NHISdata",            # Input dataset name should be in quotes as it is passed 
                                   # to Rmkf() by reference, as in SAS, instead of by value.
     group = "Population",         # Column name for grouping variable
     by = "",                      # Column name for stratification variable (none, here, 
                                   # so the empty string "" is used)
     time = "Year",                # Column name for timepoint variable
     outcome = "HYPERTEN",         # Column name for outcome variable
     se = "HYPERTEN_SE",           # Column name for SE of outcome variable
     neff =  "HYPERTEN_NEFF",      # Column name for (effective) sample size, to use when 
                                   # randomVars = TRUE. Ignored when randomVars = FALSE.
     bayesModel = "common_linear", # Default value is 'bma_cubic' if breakPoint is unspecified 
                                   # and 'bma_linear' is specified. If not "", bayesModel must 
                                   # be one of 'bma_cubic', 'indep_cubic', 'common_cubic', 
                                   # 'bma_quad', 'indep_quad', 'common_quad', 'bma_linear', 
                                   # 'indep_linear', 'common_linear', or 'dropped' 
     ARmodel = "indep_ar",         # Default value is 'common_ar'. If not "", ARmodel must be 
                                   # one of 'common_ar', 'common_arh', or 'indep_ar'.
     randomVars = TRUE,            # Both fixed (FALSE) and random (TRUE; default) sampling 
                                   # variances can be accommodated.
     comparedTo = "",              # Default value is "". When not "" (more resource-intensive), 
                                   # disparities are calculated, and comparedTo must be one of 'ALL', 
                                   # 'MIN', 'MAX', or a valid group name from the 'group' column.
     chains = 4,                   # Number of chains for MCMC (default is 4)
     thin = 5,                     # Thinning rate (steps skipped between samples) (default is 1)
     nbi = 10000,                  # Number of burn-in steps (default is 10000)
     nmc = 50000,                  # Number of MCMC samples after burn-in (default is 50000)
     modelPrint = TRUE,            # Set to TRUE to print compilation progress and basic information
                                   # regarding the MCMC.
     mcmclog = TRUE,               # Default value is FALSE. Set to TRUE (more resource-intensive) 
                                   # to retain/examine the full posterior samples in addition to 
                                   # the summaries produced by Rmkf. This is needed for producing 
                                   # diagnostic plots.
     mcmcHPDintervals = TRUE,      # Default value is FALSE. Set to TRUE to calculate the highest 
                                   # posterior density (HPD) 95% credible intervals for the scalar 
                                   # model parameters. Posterior means and SDs are always calculated.
     outloc = WORKDIR,             # Working R directory with read/write access to store output 
                                   # HTML tables and the CSV file.
     dllloc = TEMPDIR,             # Temporary R directory with read/write access to store DLLs 
                                   # and C++ files that are generated by 'nimble'.
     clearDLLs = TRUE              # Set to TRUE (default) to unload DLLs and clear C++ files from
                                   # specified temporary directory *after* user exits the R session.
)

##
## Example diagnostics plots that can be generated using 'coda' when mcmclog = TRUE,
## assuming Rmkf input parameter 'out' was set to "param", which is the default value.
##

## Autocorrelation plots for stratum 1 -- all chains
pdf(paste(WORKDIR, "param_eMKFautocorr.pdf", sep="\\"), 
    width=11, height=8.5, paper="letter")
par(mfrow=c(3,3))
coda::autocorr.plot(param_postSample$rep1, auto.layout=FALSE, ask=FALSE)
dev.off()

## Density plots for stratum 1 -- chain 3 only
pdf(paste(WORKDIR, "param_eMKFdensplot.pdf", sep="\\"), 
    width=11, height=8.5, paper="letter")
par(mfrow=c(3,3))
coda::densplot(param_postSample$rep1$chain3, show.obs = FALSE)
dev.off()

##
## Other 'coda' graphical diagnostics can be generated similarly.
##


######################################################################################
# Example from NHANES -- used in the eMKF Guidance Report (DOI: 10.15620/cdc/157496) #
######################################################################################

NHANESobesity <- read.csv(paste(DATADIR, "nhanes9920obesity.csv", sep="\\"), as.is=TRUE, fill=FALSE)

## Source code and Rmkf call with default parameter values whenever those are left unspecified
source(paste(MACRDIR, "emkf_macro_v24.R", sep="\\"))
Rmkf(out = "obes",
     data = "NHANESobesity",
     group = "Population", 
     by = "Age", 
     time = "Year", 
     outcome = "Obesity",
     se = "SE_obesity",
     neff = "NEFF_obesity", 
     bayesModel = "bma_quad",
     comparedTo = "MIN",
     outloc = WORKDIR,
     dllloc = TEMPDIR
)


#######################################################################
# Simulated data with level shift in trend -- adapted from the NHANES #
#  data used in the eMKF Evaluation Report (DOI: 10.15620/cdc/157497) #
#######################################################################

## Read-in simulated data from common or independent trends across groups 
## (uncomment desired selection)

NHANESobesity <- read.csv(paste(DATADIR, "nhanes9923_simtrends_redesign.csv", sep="\\"), as.is=TRUE, fill=FALSE)
#NHANESobesity <- read.csv(paste(DATADIR, "nhanes9923_simtrends_grpsp_redesign.csv", sep="\\"), as.is=TRUE, fill=FALSE)

## Uncomment desired "true" underlying trend from which data were simulated
#NHANESobesity1 <- NHANESobesity[NHANESobesity$trend == "cubic",]
NHANESobesity1 <- NHANESobesity[NHANESobesity$trend == "quadratic",]
#NHANESobesity1 <- NHANESobesity[NHANESobesity$trend == "linear",]

## Source code and Rmkf call with default parameter values whenever those are left unspecified
source(paste(MACRDIR, "emkf_macro_v24.R", sep="\\"))
Rmkf(out = "break",
     data = "NHANESobesity1",
     group = "raceethn",
     by = "agegroup2",
     time = "year",
     breakPoint = 2012,            # The breakPoint argument specifies the timepoint at which 
                                   # segment 2 *starts*, defaulting to NA when unspecified.
     breakType = "level_break",    # Unless set to 'full_break' by the user, breakType defaults to 
                                   # 'level_break' when breakPoint is specified.
     outcome = "obesity",
     se = "se.obesity",
     neff = "neff.obesity",
     bayesModel = "bma_cubic",     # The default value for bayesModel is 'bma_linear' when
                                   # breakPoint is specified, but this can be overridden by the user.
                                   # Rmkf will check to make sure the requested model and polynomial
                                   # trend can be fit given the number of timepoints in each segment.
     thin = 2,
     nbi = 20000,
     nmc = 100000,
     mcmcHPDintervals = TRUE,
     
     outloc = WORKDIR,
     dllloc = TEMPDIR
)

## MT. 30-Jan-2026.
