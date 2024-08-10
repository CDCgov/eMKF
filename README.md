# The Enhanced Modified Kalman Filter (eMKF) tool for small domain estimation

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Overview

This project contains the SAS code to implement the **enhanced Modified Kalman Filter (eMKF)**. The enhanced MKF procedure enables production of model-based estimates for small populations where direct estimates may lack precision, improving the availability of data for assessing and monitoring health disparities. The enhanced MKF procedure and macro build on the earlier Modified Kalman Filter procedure (Setodji, Lockwood, McCaffrey, Elliott & Adams, 2011) to accommodate nonlinear time trends, irregularly spaced time points, and random sampling variances for the underlying population subgroup means, rates, or proportions. Bayesian estimation in the eMKF macro is implemented adaptably and transparently using PROC MCMC and related SAS 9.4 procedures. Model averaging in the eMKF macro uses a Bayesian mixture prior approach, and renders predictions more robust to polynomial trend misspecification. Various other features in the eMKF macro also improve its functionality, flexibility, and usability relative to the earlier macro; an outline of these improvements is included below. 

This project contains the SAS macro to implement the eMKF, [emkf_macro.sas](emkf_macro.sas), along with several examples of SAS code to implement the eMKF macro in the [Testing-and-implementation](Testing-and-implementation) folder, with some sample data sets included in [Sample-data-files](Sample-data-files). All data included here are public-use and/or simulated data. 

### Requires

* SAS 9.4

### Running the eMKF macro

Sample SAS programs illustrating how to run the eMKF macro are included in [Testing-and-implementation](Testing-and-implementation). 

* Please refer to the SAS file [Compare-eMKF-to-MKF.sas](Testing-and-implementation/Compare-eMKF-to-MKF.sas), for sample code to implement the eMKF and replicate the functionality of the older MKF macro.
* Please refer to the SAS file [Example-eMKF-NHANES-Data.sas](Testing-and-implementation/Example-eMKF-NHANES-Data.sas), for sample code to implement the eMKF to estimate trends in adult obesity prevalence by race/ethnicity using data from the National Health and Nutrition Examination Survey, 1999-March 2020.
* Please refer to the SAS file [Example-eMKF-Simulated-Mortality-Data.sas](Testing-and-implementation/Example-eMKF-Simulated-Mortality-Data.sas), for sample code to implement the eMKF to estimate trends in mortality rates due to external causes of death by state, age, and race/ethnicity using simulated mortality data from the National Vital Statistics System, 1999-2020.

#### eMKF macro requirements

Input data to the enhanced MKF macro are required to be in stacked (long) format, with each row representing a time- and group-specific estimate. Additional columns required include timepoint and population group identifiers (e.g., racial/ethnic group), standard errors (SE), and (effective) sample sizes for survey data. The variable for timepoint (typically year or survey cycle) must be numeric. A stratification variable can also be included (e.g., age group). 

Other requirements of the eMKF macro:
* There should be at least *k+4* time points to estimate a degree *k* polynomial trend.
  * For linear trends, 5 timepoints are needed.
  * For quadratic trends, 6 time points are needed.
  * For cubic trends, 7 time points are needed.
* There should be no missing data for any subgroup. For example, if there were 0 individuals sampled for a given group for one timepoint, analysts will need to aggregate the data over larger time periods or groups to ensure that there are no group-by-time cells without any data/sample.
* Estimates of 0 are allowed, but 0 SEs are imputed by the macro using the average SE across the non-zero SEs for that group/stratum. However, if estimates and corresponding SEs for a given group/stratum are 0 across all timepoints included, analysts should consider aggregating the data into larger groups to ensure that there are some non-zero estimates for all groups/strata. As a preliminary solution, the macro will impute such 0 SEs using the average SE across strata for each given timepoint.

Please refer to the documentation provided in <mark>[Link to Series 2 report]</mark> for more details. 

## Related documents

* <mark>[Link to Series 2 report]</mark>
* <mark>[Link to Series 2 report]</mark>
* <mark>[Link to Appendix Table II-1 in Guidance report] (Default eMKF macro parameter settings)</mark>


### Outline of methodological differences between eMKF (version 1.4 2024-08-10) and the original MKF

* Time points
  * eMKF allows for time points to be unequally spaced.
  * In the MLE-based setting, eMKF updates the recursion formulas used to determine the MSE-optimal estimators (BLUPs) of the true states to allow for an arbitrary lag between successive time points instead of lag=1 (see Lockwood et al 2011; DOI: 10.1002/sim.3897).

* Polynomial trends
  * eMKF allows for quadratic and cubic time trends to be fitted in both the Bayesian and MLE-based estimation settings.
  * By default, eMKF does not allow fitting a degree k polynomial trend (k=0,1,2,3) unless there are k+4 available time points.
  * eMKF returns an error if there is only one available time point per group.
  * eMKF pre-transforms trend coefficients using an orthogonal polynomial design matrix for comparability of coefficients in linear, quad, and cubic trend models. The coefficients are reverse-transformed before the program exits so the user only sees the "raw" coefficients.

* Sampling variances
  * eMKF allows for random sampling variances in the Bayesian setting (see Polettini 2017; DOI 10.1214/16-BA1019, for an overview).
  * eMKF uses (effective) sample size (neff) as degrees of freedom in chi-squared distribution of group- and time-specific sampling variance, therefore neff must be supplied. 

* Disparities calculations
  * In the Bayesian setting, eMKF estimates all pairwise differences and ratios between groups at the latest time point.
  * Additional disparities measures (highest and lowest rates, maximal rate difference and ratio, summary rate difference and ratio) are also calculated.
  * The user can further calculate any other measures he/she desires from the posterior draws, which the user can request to be saved to the workspace.
  * These features differ from MKF where only pairwise differences were calculated and the full MCMC samples from the joint posterior distribution were not available.

* Bayesian estimation setting
  * eMKF implements "independent" and "common" trend options, which were left out of the RAND version of the macro, in addition to the "full" hierarchical Bayesian model (see RAND User's Guide "TR997_compiled").
  * eMKF implements Bayesian model averaging using a mixture prior approach, up to linear (3 possible models), quadratic (5 possible models), or cubic (7 possible models).
  * eMKF replaces the call to the external .exe file (which consisted of pre-compiled C code) with a call to PROC MCMC.
  * eMKF implements Gibbs sampling in PROC MCMC by calling user defined samplers (UDSs) that are custom-built and precompiled using PROC FCMP.
  * eMKF replaces z-score-based convergence diagnostic used in MKF with robust version of the Gelman-Rubin diagnostic (see Vehtari et al 2021; DOI 10.1214/20-BA1221).
  * eMKF applies the Gelman-Rubin diagnostic to all model parameters, not just for the true state predictions (etas) as in the original MKF.
  * eMKF defaults to more stringent threshold of 1.01 instead of 1.10 for the Gelman-Rubin diagnostic, as per recommendation in Vehtari et al (2021). 
  * eMKF defaults to 4 chains instead of 3, as per recommendation in Vehtari et al (2021). Each chain is further split in 2 to compute the Gelman-Rubin diagnostic.
  * eMKF uses closed-form expressions for the determinant, inverse, and Cholesky decomposition of the AR variance-covariance matrix whenever possible to speed up calculations.
  * eMKF allows the user to select the built-in slice sampler in PROC MCMC to use instead of the traditional random walk MH sampler for sampling AR parameters + SD hyperparameters.

* MLE-based estimation setting
  * In eMKF, any subset of the seven allowable models (indep_cubic, _quad, _linear; common_cubic, _quad, _linear; and dropped) can be averaged. 
  * However, the code checks the specified models and adds a common "descendent" if is not already included, so as to have a reference model for Bayes factors.
  * eMKF increases the maxiter option for PROC NLMIXED to 400 instead of 200 (default) when dealing with two outcomes to improve convergence when k = 2,3.
  * eMKF initializes parameters to pass to PROC NLMIXED using the appropriate 'by' group stratum/replication (PROC REG). This differs from MKF where only the first stratum/replication was used to initialize the regression coefficients across strata/replications.
  * eMKF initializes parameters to pass to PROC NLMIXED using the appropriate degree k polynomial regression (PROC REG), including for k=0. This differs from MKF where for k=0 (dropped), the intercept values were initialized at those from the linear regression y=a+b*t instead of y=a.
  * For the dropped (k=0) case, eMKF only keeps the column vector of 1s in the X matrix (and subsequent matrix calculations for the MSE). This differs from MKF where both the 1s and ts column were kept in the X matrix.
  * The number of groups was limited to 15 in the original MKF due to the use IML function block to create block diagonal matrices in the calculation of MSEs; this was corrected in eMKF to allow an arbitrary number of groups.

* Macro usability
  * eMKF includes extensive comments and streamlines the code for readability.
  * eMKF allows the user additional flexibility in customizing model output and diagnostics, and streamlines the SAS workspace.
  * eMKF checks for errors in macro parameter specification, including length of character strings for prefix of output datasets and the maximum number of groups and/or data points that the code implementation can handle.
  * Std. Error label in output table was replaced with RMSE in eMKF to avoid confusion.
  * Zero SEs and effective sample sizes (when applicable) are imputed using the average across timepoints for the given group and stratum. Any remaining zero SEs and effective sample sizes (when applicable) are imputed timepoint by timepoint using the average across strata for the given group. 


### References

Lockwood JR, McCaffrey DF, Setodji CM, Elliott MN. Smoothing across time in repeated cross-sectional data. Stat Med 30(5):584–94. 2011. DOI: https://dx.doi.org/10.1002/sim.3897.

Polettini S. A generalised semiparametric Bayesian Fay–Herriot model for small area estimation shrinking both means and variances. Bayesian Anal 12(3):729–52. 2016. DOI: https://dx.doi.org/10.1214/16-BA1019.
    
Rossen LM, Talih M, Patel P, Earp M, Parker JD. Evaluation of an enhanced modified Kalman filter approach for estimating health outcomes in small subpopulations. National Center for Health Statistics. Vital Health Stat 2(208). 2024. DOI: https://dx.doi.org/10.15620/cdc/157496.

Setodji CM, Lockwood JR, McCaffrey DF, Elliott MN, Adams JL. The Modified Kalman Filter macro: User’s guide. RAND Technical Report No. TR-997-DHHS. 2011. Available from: https://www.rand.org/pubs/technical_reports/TR997.html.

Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical guidance for using the modified Kalman filter in small-domain estimation at the National Center for Health Statistics. National Center for Health Statistics. Vital Health Stat 2(209). 2024. DOI: https://dx.doi.org/10.15620/cdc/157496.

Vehtari A, Gelman A, Simpson D, Carpenter B, Bürkner PC. Rank-normalization, folding, and localization: An improved Ȓ for assessing convergence of MCMC (with discussion). Bayesian Anal 16(2):667–718. 2021. DOI: https://dx.doi.org/10.1214/20-BA1221.

  
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
