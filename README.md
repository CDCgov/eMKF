# The Enhanced Modified Kalman Filter (eMKF) tool for small domain estimation [version 2.4 2026-01-30]

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the CDC mission. GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Overview

This project contains the SAS and R code to implement the **enhanced Modified Kalman Filter (eMKF)**. The enhanced MKF procedure enables production of model-based estimates for small populations where direct estimates may lack precision, improving the availability of data for assessing and monitoring health disparities. The enhanced MKF procedure and macro build on the earlier Modified Kalman Filter procedure (Setodji et al, 2011; Lockwood et al, 2011) to accommodate nonlinear time trends, irregularly spaced time points, and random sampling variances for the underlying population subgroup means, rates, or proportions. In version 2.4, the eMKF macro also allows for a trend break to be specified at a given timepoint, to capture survey redesigns or other interruptions or changes in the underlying data. Additionally, the eMKF macro version 2.4 ensures the autocorrelation coefficient for the AR(1) random effects remains nonnegative, which is required when time points are fractional (e.g., 2018.6). Bayesian estimation in the SAS eMKF macro is implemented adaptably and transparently using PROC MCMC and related SAS 9.4 procedures. The R eMKF macro ('Rmkf') uses R package 'nimble' to implement Bayesian estimation. Model averaging in the SAS and R eMKF macros uses a Bayesian mixture prior approach and renders predictions more robust to polynomial trend misspecification. Various other features in the SAS eMKF macro also improve its functionality, flexibility, and usability relative to the earlier SAS MKF macro; an outline of these improvements is included below, and comprehensive technical guidance for using the SAS eMKF macro is provided in Talih et al (2024). An evaluation of the eMKF approach in the context of small subpopulation data from the National Center for Health Statistics is available in Rossen et al (2024). A comparison between the R and SAS implementations of the eMKF macro version 2.4 is also provided below.

The SAS code to implement the eMKF macro version 2.4 is available from [emkf_macro_v24.sas](SAS-macro/emkf_macro_v24.sas), along with several use cases (in the [Testing-and-implementation](Testing-and-implementation) folder) and sample data sets (in the [Sample-data-files](Sample-data-files) folder). All data included here are public-use and/or simulated data. The R prototype of the eMKF macro for version 2.4 is available from [emkf_macro_v24.R](R-macro/emkf_macro_v24.R). The R implementation ('Rmkf') relies on the 'nimble' package; see below for more information. A tutorial on how to use the R eMKF macro can be found in [emkf_macro_v24_tutorial](R-macro/emkf_macro_v24_tutorial.R).

### Requires

* SAS 9.4 for the SAS eMKF macro
* R 4.4.0 or higher for the R eMKF macro

### Running the SAS eMKF macro

Sample SAS programs illustrating how to run the eMKF macro are included in [Testing-and-implementation](Testing-and-implementation) folder. 

* Please refer to the SAS file [Compare-eMKF-to-MKF.sas](Testing-and-implementation/Compare-eMKF-to-MKF.sas) for sample code to implement the eMKF and replicate the functionality of the earlier MKF macro.
* Please refer to the SAS file [Example-eMKF-NHANES-Data.sas](Testing-and-implementation/Example-eMKF-NHANES-Data.sas) for sample code to implement the eMKF to estimate trends in adult obesity prevalence by race/ethnicity using data from the National Health and Nutrition Examination Survey, 1999–March 2020.
* Please refer to the SAS file [Example-eMKF-Simulated-Mortality-Data.sas](Testing-and-implementation/Example-eMKF-Simulated-Mortality-Data.sas) for sample code to implement the eMKF to estimate trends in mortality rates due to external causes of death by state, age, and race/ethnicity using simulated mortality data from the National Vital Statistics System, 1999–2020.
  * Users are advised to consider log-transforming vital rates (after possibly adding a small offset, say 0.5, to avoid taking the logarithm of zero) prior to running the eMKF macro, as the underlying model assumes normal distributions for the fixed and random effects. 
* Please refer to the SAS file [Example-eMKF-NSFG-Data.sas](Testing-and-implementation/Example-eMKF-NSFG-Data/Example-eMKF-NSFG-Data.sas) and R file [Example-eMKF-NSFG-Data.R](Testing-and-implementation/Example-eMKF-NSFG-Data/Example-eMKF-NSFG-Data.R) for sample code to implement the eMKF to estimate trends in risk of pregnancy loss by metropolitan status and maternal age using data from the National Survey of Family Growth, 2006–2019.  
  * Forrest SE, Rossen LM, Ahrens KA. Trends in Risk of Pregnancy Loss Among US Women by Metropolitan Status, 2000–2018. Paediatric and Perinatal Epidemiology. 2025. <https://doi.org/10.1111/ppe.70066>.

### Running the R eMKF macro

A tutorial with several examples on how to run the R eMKF macro ('Rmkf') is available from [emkf_macro_v24_tutorial.R](R-macro/emkf_macro_v24_tutorial.R).

### eMKF macro requirements

Input data to both the SAS and R eMKF macros are required to be in stacked (long) format, with each row representing a time- and group-specific estimate. Additional columns required include timepoint and population group identifiers (e.g., racial/ethnic group), standard errors (SEs), and (effective) sample sizes (NEFFs) for survey data. The variable for timepoint (typically year or survey cycle) must be numeric. A stratification variable can also be included (e.g., age group). 

#### Other requirements of the SAS and R eMKF macros:
* There should generally be at least *k+4* time points to estimate a degree *k* polynomial trend. eMKF version 2.4 will check the model specification against the number of available data points to avoid over-fitting.
  * For linear trends, 5 timepoints are needed.
  * For quadratic trends, 6 time points are needed.
  * For cubic trends, 7 time points are needed.
* For trend-break modeling in eMKF version 2.4, only one breakpoint can be specified; it must be one of the values of the time variable/column in the input dataset, and the macro will interpret it as the starting point of the second segment. 
  * Polynomial trend combinations for segments 1 and 2 will be determined by the type of trend break specified (level-shift only vs. full trend break), the specified trend model for segment 1, and the number of available timepoints; see below.
* There should be no missing data for any subgroup. For example, if there were 0 individuals sampled for a given group for one timepoint, analysts will need to aggregate the data over larger time periods or groups to ensure that there are no group-by-time cells without any data/sample.
* Estimates of 0 are allowed, but 0 SEs are imputed by the macro using the average SE across the non-zero SEs for that group/stratum. However, if estimates and corresponding SEs for a given group/stratum are 0 across all timepoints included, analysts should consider aggregating the data into larger groups to ensure that there are some non-zero estimates for all groups/strata. As a preliminary solution, the macro will impute such 0 SEs using the average SE across strata for each given timepoint.
* As of version 2.4 of the SAS eMKF macro, the number of groups cannot be larger than 200 and the number of groups times the number of time points cannot exceed 5,000 per stratum; see below. There are no similar restrictions in the R eMKF macro.

Please refer to the documentation provided in the eMKF Guidance Report (doi:10.15620/cdc/157496) for more details; see <https://www.cdc.gov/nchs/data/series/sr_02/sr02-209.pdf>. 

### Suggested citations

Talih M, Patel P, Rossen LM. An R-NIMBLE implementation of the enhanced modified Kalman filter (eMKF) tool for small domain estimation [version 2.4 2026-01-30]. National Center for Health Statistics. 2026. Available from: <https://github.com/CDCgov/eMKF>.

Talih M, Patel P, Rossen LM. The enhanced modified Kalman filter (eMKF) tool for small domain estimation [version 2.4 2026-01-30]. National Center for Health Statistics. 2026. Available from: <https://github.com/CDCgov/eMKF>.

Talih M, Rossen LM, Patel P, Earp M, Parker JD. The enhanced modified Kalman filter (eMKF) tool for small domain estimation [version 1.4 2024-08-10]. National Center for Health Statistics. 2024. Available from: <https://github.com/CDCgov/eMKF>.

## Related documents

* eMKF Guidance Report (doi:10.15620/cdc/157496); see <https://www.cdc.gov/nchs/data/series/sr_02/sr02-209.pdf>.
   * Appendix II describes the parameter settings, defaults and functionality for the SAS eMKF macro version 1.4.
* eMKF Evaluation Report (doi:10.15620/cdc/157497); see <https://www.cdc.gov/nchs/data/series/sr_02/sr02-208.pdf>.

## Methodological differences between the SAS eMKF macro and the earlier SAS MKF macro

### Time points
* eMKF allows for time points to be unequally spaced as well as fractional. As of version 2.4, to ensure a nonnegative autocorrelation coefficient for the AR(1) random effects, a truncated normal prior distribution is used. Previously, negative draws were discarded as they resulted in invalid AR(1) covariance matrices when time points were fractional, e.g., 2018.6.
* In the maximum likelihood estimation (MLE) setting, the SAS eMKF macro updates the recursion formulas used to determine the mean squared error (MSE)-optimal estimators (also known as best linear unbiased predictors or BLUPs) of the true states to allow for an arbitrary lag between successive time points instead of lag=1 (Lockwood et al, 2011).

### Polynomial trends
* The SAS eMKF macro allows for quadratic and cubic time trends to be fitted in both the Bayesian and MLE-based estimation settings.
* By default, eMKF does not allow fitting a degree k polynomial trend (k=0,1,2,3) unless there are k+4 available time points.
* eMKF returns an error if there is only one available time point per group.
* eMKF pre-transforms trend coefficients using an orthogonal polynomial design matrix for comparability of coefficients in linear, quadratic, and cubic trend models. In the SAS eMKF macro, the coefficients are reverse-transformed before the program exits so the user only sees the "raw" or untransformed coefficients.

### Sampling variances
* eMKF allows for random sampling variances in the Bayesian setting; see Polettini (2017) for an overview.
* eMKF uses (effective) sample size (NEFF) as degrees of freedom in chi-squared distribution of group- and time-specific sampling variances, therefore NEFF must be supplied when survey data are used. 

### Disparities calculations
* In the Bayesian setting, eMKF estimates all pairwise differences and ratios between groups at the latest time point.
* Additional disparities measures (highest and lowest rates, maximal rate difference and ratio, summary rate difference and ratio) are also calculated.
* The user can further calculate any other measures desired from the posterior draws, which the user can request to be saved to the workspace.
* These features differ from the earlier MKF where only pairwise differences were calculated and the full Markov chain Monte Carlo (MCMC) samples from the joint posterior distribution were not available.

### Bayesian estimation setting
* The SAS eMKF macro implements "independent" and "common" trend options, which were left out of the earlier MKF macro, in addition to the "full" hierarchical Bayesian model (Setodji et al, 2011).
* eMKF implements Bayesian model averaging using a mixture prior approach, up to linear (3 possible models), quadratic (5 possible models), or cubic (7 possible models).
* The SAS eMKF macro replaces the call to the external .exe file in the earlier MKF macro (which consisted of pre-compiled C code) with a call to SAS PROC MCMC.
* The SAS eMKF macro implements Gibbs sampling in SAS PROC MCMC by calling user defined samplers (UDSs) that are custom-built and precompiled using SAS PROC FCMP.
* eMKF replaces the z-score-based convergence diagnostic used in the earlier MKF with a robust version of the Gelman-Rubin diagnostic (Vehtari et al, 2021).
* eMKF applies the Gelman-Rubin diagnostic to all model parameters, not just for the true state predictions (etas) as in the earlier MKF.
* eMKF defaults to more stringent threshold of 1.01 instead of 1.10 for the Gelman-Rubin diagnostic, as per recommendation in Vehtari et al (2021). 
* eMKF defaults to 4 chains instead of 3, as per recommendation in Vehtari et al (2021). Each chain is further split in 2 to compute the Gelman-Rubin diagnostic.
* eMKF uses closed-form expressions for the determinant, inverse, and Cholesky decomposition of the autoregressive (AR) variance-covariance matrix whenever possible to speed up calculations.
* The SAS eMKF macro allows the user to select the built-in slice sampler in SAS PROC MCMC to use instead of the traditional random walk Metropolis-Hastings sampler for sampling AR parameters and standard deviation hyperparameters.

### MLE-based estimation setting
* In the SAS eMKF macro, any subset of the seven allowable models (indep_cubic, indep_quad, indep_linear; common_cubic, common_quad, common_linear; and dropped) can be averaged; see [Appendix II Table in Guidance report](https://www.cdc.gov/nchs/data/series/sr_02/sr02-209.pdf]). 
* However, the code checks the specified models and adds a common "descendent" if is not already included, so as to have a reference model for Bayes factors. Specifically:
    * If both the indep_quad and common_cubic models are specified, then the common_quad model will be added if needed.
    * If both the indep_linear and common_cubic models are specified, then the common_linear model will be added if needed.
    * If both the indep_linear and common_quad models are specified, then the common_linear model will be added if needed.  
* The SAS eMKF macro increases the maximum iteration (maxiter) option for SAS PROC NLMIXED to 400 instead of 200 (SAS default) when dealing with two outcomes to improve convergence when k = 2,3.
* The SAS eMKF macro initializes parameters to pass to SAS PROC NLMIXED using the appropriate 'by' group stratum/replication (SAS PROC REG). This differs from the earlier MKF macro where only the first stratum/replication was used to initialize the regression coefficients across strata/replications.
* The SAS eMKF macro initializes parameters to pass to SAS PROC NLMIXED using the appropriate degree k polynomial regression (SAS PROC REG), including for k=0. This differs from the earlier MKF macro where for k=0 (dropped), the intercept values were initialized at those from the linear regression y=a+b×time instead of y=a.
* For the dropped (k=0) case, the SAS eMKF macro only keeps the column vector of 1s in the X matrix (and subsequent matrix calculations for the MSE). This differs from the earlier MKF macro where both the 1s and times (ts) column were kept in the X matrix.
* The number of groups was limited to 15 in the earlier SAS MKF macro due to the use of SAS PROC IML's function BLOCK to create block diagonal matrices in the calculation of MSEs; this was corrected in the SAS eMKF macro to allow an arbitrary number of groups.

### Macro usability
* The SAS eMKF macro includes extensive comments and streamlines the code for readability.
* The SAS eMKF macro allows the user additional flexibility in customizing model output and diagnostics, and streamlines the SAS workspace.
* The SAS eMKF macro checks for errors in macro parameter specification, including length of character strings for prefix of output datasets and the maximum number of groups and/or data points that the code implementation can handle.
* The "Std. Error" label in output tables was replaced with "RMSE" (root mean squared error) in eMKF to avoid confusion.
* Zero SEs and effective sample sizes (when applicable) are imputed using the average across timepoints for the given group and stratum. Any remaining zero SEs and effective sample sizes (when applicable) are imputed timepoint by timepoint using the average across strata for the given group. 

## New in eMKF macro version 2.4
* eMKF version 2.4 allows for a trend break to be specified at a given break point T in both the Bayesian (SAS and R implementations) and MLE settings (SAS implementation only).
   * Break point T must be one of the values of the time variable/column.
   * Segment 1 will consist of timepoints 1, 2, ..., T-1
   * Segment 2 will consist of timepoints T, T+1, ..., n.
   * Trend break will be either a level-shift (default) or a full trend break.
* Polynomial trend combinations for segments 1 and 2 will be determined by the type of trend break specified (level-shift only vs. full trend break), the specified trend model for segment 1, and the number of available timepoints.
  * Models for segments 1 and 2 will be either both bma-, both full-, both indep-, both common-, or both dropped (intercepts-only) models.
  * The polynomial degree for segment 2’s trend will not exceed that of segment 1.
  * For each segment, the polynomial degree will only be as large as can be supported by the available number of timepoints in that segment, and will also take into account the random effects AR model option specified by the user.
* eMKF version 2.4 allows specifying a common AR autocorrelation parameter but independent AR variance parameters by group.
  * This new option (ARmodel = common_arh) offers an in-between model for the random effects, adding to the common AR parameters option (common_ar) and the fully independent AR parameters option (indep_ar).
* SAS eMKF macro version 2.4 allows for the random effects variance-covariance matrix in the MLE setting (SAS implementation only) to be structured according to the common_ar, common_arh, or indep_ar structures.
  * In the SAS eMKF macro version 1.4, only the common_ar structure was available in the MLE setting for use with SAS PROC NLMIXED. 
  * Additionally, SAS eMKF macro version 2.4 users can no longer specify values for the AR autocorrelation and variance parameters to pass to SAS PROC NLMIXED; those will be estimated.
* eMKF version 2.4 uses a truncated normal prior distribution to ensure the autocorrelation coefficient(s) for the AR(1) random effects remain between 0 and 1 instead of -1 and 1. This avoids the need to discard negative draws that resulted in invalid AR(1) covariance matrices when time points were fractional, e.g., 2018.6, because raising a negative number to a fractional power is a mathematically ill-defined operation.
* To streamline calculations, the SAS eMKF macro version 2.4 allows the user to omit model datasets with fit details/statistics from SAS PROC NLMIXED, including the covariance matrix of parameters, which requires reversing the orthogonal transformation prior to exiting.
* SAS eMKF macro version 2.4 also gives users more control over the choice of optimization algorithm used in SAS PROC NLMIXED in the MLE setting.
* SAS eMKF macro version 2.4 streamlines code by using the macro IN operator instead of repeatedly compiling strings of OR conditions.
* SAS eMKF macro version 2.4 streamlines output by allowing the user to suppress the voluminous runtime notes from being written to the log.
* SAS eMKF macro version 2.4 also adds user-facing status updates on the progress of the macro to the log.
* Whereas SAS eMKF version 1.4 still allowed users to provide data in "format 1", eMKF version 2.4 only allows for "format 2".
* To streamline code that could result in macro variables exceeding the SAS character limit, SAS eMKF macro version 2.4 revises cutoffs from 204 to 200 for the number of groups and 5508 to 5000 for the number of data points.

## Main differences between the R and SAS implementations of eMKF macro version 2.4                                                                
* Only one outcome variable can be specified in the R eMKF macro.
* Only Bayesian estimation is available in the R eMKF macro; the MLE option is not available.
* Only one Bayesian model can be specified in the R eMKF macro; Bayesian model averaging (BMA) must be requested via the keywords bma_***.
* The fully-Bayesian model (from the earlier SAS MKF macro), available in the SAS eMKF macro, is not implemented in the R eMKF macro.
* The slice sampler option is not implemented in the R eMKF macro.
* The R eMKF macro uses a reflective random walk sampler for bounded variables, which is not available in the SAS eMKF macro.
* The R eMKF macro uses a target acceptance rate of 0.44 for the random walk sampler (this is hard-coded into nimble::sampler_RW using optimalAR <- 0.44), whereas the SAS eMKF macro allows the user to specify lower target acceptance rates.
* The R eMKF macro implements standalone samplers for the intercepts in BMA cases, whereas the SAS eMKF macro uses joint samplers for the intercepts and regression coefficients.
* Unlike the SAS eMKF macro, the R eMKF macro does not reverse-transform the regression coefficients after orthogonal polynomials are requested (default); users should interpret the model coefficients accordingly.
* Unlike the SAS eMKF macro, the R eMKF macro does not require a maximum length for dataset or column names.
* Unlike the SAS eMKF macro, the R eMKF macro does not require a maximum number of groups or data points.
* Unlike the SAS eMKF macro, the R eMKF macro gives its users the option to calculate the HPD credible intervals for all scalar variables, in addition to the posterior means and SDs of those variables.
* Unlike the SAS eMKF macro, the R eMKF macro does not implement graphical diagnostics directly; however, users can call 'coda' diagnostic plots on the saved posterior sample; see the examples provided in [emkf_macro_v24_tutorial.R](R-macro/emkf_macro_v24_tutorial.R).

### References

de Valpine, P., D. Turek, C.J. Paciorek, C. Anderson-Bergman, D. Temple Lang, and R. Bodik. 2017. Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26: 403-413. <https://doi.org/10.1080/10618600.2016.1172487>.

de Valpine P, Paciorek C, Turek D, Michaud N, Anderson-Bergman C, Obermeyer F, Wehrhahn Cortes C, Rodrìguez A, Temple Lang D, Paganin S (2025). _NIMBLE: MCMC, Particle Filtering, and Programmable Hierarchical Modeling_.  <https://doi.org/10.5281/zenodo.1211190>, R package version 1.4.0, <https://cran.r-project.org/package=nimble>.

de Valpine P, Paciorek C, Turek D, Michaud N, Anderson-Bergman C, Obermeyer F, Wehrhahn Cortes C, Rodrìguez A, Temple Lang D, Paganin S (2025). _NIMBLE User Manual_. <https://doi.org/10.5281/zenodo.1211190>, R package manual version 1.4.0, <https://r-nimble.org>.

Lockwood JR, McCaffrey DF, Setodji CM, Elliott MN. Smoothing across time in repeated cross-sectional data. Stat Med 30(5):584–94. 2011. <https://dx.doi.org/10.1002/sim.3897>.

Polettini S. A generalised semiparametric Bayesian Fay–Herriot model for small area estimation shrinking both means and variances. Bayesian Anal 12(3):729–52. 2016. <https://dx.doi.org/10.1214/16-BA1019>.
    
Rossen LM, Talih M, Patel P, Earp M, Parker JD. Evaluation of an enhanced modified Kalman filter approach for estimating health outcomes in small subpopulations. National Center for Health Statistics. Vital Health Stat 2(208). 2024. doi:10.15620/cdc/157496. Available from: <https://www.cdc.gov/nchs/data/series/sr_02/sr02-208.pdf>.

Setodji CM, Lockwood JR, McCaffrey DF, Elliott MN, Adams JL. The Modified Kalman Filter macro: User’s guide. RAND Technical Report No. TR-997-DHHS. 2011. Available from: <https://www.rand.org/pubs/technical_reports/TR997.html>.

Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical guidance for using the modified Kalman filter in small-domain estimation at the National Center for Health Statistics. National Center for Health Statistics. Vital Health Stat 2(209). 2024. doi:10.15620/cdc/157496. Available from: <https://www.cdc.gov/nchs/data/series/sr_02/sr02-209.pdf>.

Vehtari A, Gelman A, Simpson D, Carpenter B, Bürkner PC. Rank-normalization, folding, and localization: An improved Ȓ for assessing convergence of MCMC (with discussion). Bayesian Anal 16(2):667–718. 2021. <https://dx.doi.org/10.1214/20-BA1221>.

  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
