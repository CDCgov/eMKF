# SAS eMKF Macro Version 2.4

## Overview

See the eMKF Guidance Report (DOI: [10.15620/cdc/157496](https://stacks.cdc.gov/view/cdc/157496)) for background and how to use the SAS eMKF macro.

The NCHS eMKF macro expands RAND's MKF macro to:

- **[NEW in v2.4]** Allow for trend break (level shift or full break) at a specified break point
- **[NEW in v2.4]** Better accommodate fractional timepoints, restrict AR(1) auto-correlation parameter for random effects to the interval [0,1) instead of (-1,1)
- **[NEW in v2.4]** Implement truncated normal random number generator with PROC FCMP and use truncated normal priors for all AR(1) parameters and hyper-parameters
- Allow for unequally-spaced time points
- Allow for quadratic and cubic trends in both MLE-based and Bayesian estimation settings
- Allow for model averaging over (orthogonal) polynomial trend model sequences up to cubic in both MLE-based and Bayesian settings
- Allow for common as well as group-specific AR parameters rho and tausq for the random effects in the Bayesian setting and **[NEW in v2.4]** MLE-based setting
- **[NEW in v2.4]** Added specification for common autocorrelation parameter rho but independent variance parameters tausq by group for the random effects
- Allow for random sampling variances in the Bayesian setting
- Implement Gibbs sampling in PROC MCMC using user defined samplers (UDSs) that are precompiled with PROC FCMP, replacing .exe file (C code)
- Expand calculations of between-group disparities at the latest time point in the Bayesian setting

**Original macro defined:** 03-16-2007, Revised 05-08-09 D. McCaffrey. Revised 03-15-2010 [RAND Statistics Group]  
**MKF macro + eMKF expansion:** 2023 Q1-Q2 by M. Talih + expanded in 2025 Q1-Q2 to allow for trend break specification (eMKF v2.4)

## Macro Definition

```sas
%macro mkf(
    /* Core Data Parameters */
    data        = ,                  /* Name of the dataset. Must be in long/stacked format */
    outcome     = ,                  /* The outcome of interest */
    se          = ,                  /* The standard error of the outcome */
    neff        = ,                  /* (Effective) sample sizes for random sampling variance case */
    outcome2    = ,                  /* Second outcome (if applicable) */
    se2         = ,                  /* Standard error of second outcome */
    neff2       = ,                  /* Effective sample size for second outcome */
    by          = ,                  /* Allows models to run for multiple strata; can be used for simulations */
    group       = ,                  /* The different groups (e.g., race/ethnicity group) variable */
    time        = ,                  /* The time variable */
    
    /* Trend Break Parameters [NEW in v2.4] */
    breakPoint  = ,                  /* Specifies value of timepoint at which trend break occurs (default is none) */
    breakType   = ,                  /* level_break (default) or full_break */
    
    /* Modeling Options */
    slopes      = ,                  /* Type of slope: indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, dropped, or empty */
    Bayesmodel  = ,                  /* Bayesian model type: defaults to bma_cubic when breakPoint is empty, bma_linear otherwise */
    BayesmodelAvg = YES,            /* Whether to perform Bayesian model averaging */
    randomVars  = YES,              /* Treat sampling variances as random (YES) or fixed (NO) */
    ARmodel     = common_ar,        /* AR(1) model: common_ar (default), common_arh, or indep_ar */
    
    /* Output Options */
    xtrakeep    = ,                  /* Any variable to keep in the data while running models */
    out         = param,             /* Prefix for output datasets */
    nlmixedDetails = NO,            /* Whether to include model datasets with fit details/statistics */
    comparedto  = ,                  /* Options for estimating disparities/differences in Bayesian setting */
    comparedata = ,                  /* Dataset for disparity calculations */
    modelprint  = NO,               /* YES to print proc nlmixed/mcmc results */
    finalprint  = YES,              /* YES to print MKF estimates for last time point */
    
    /* Output Formatting */
    pdigit      = 4,                /* Number of decimal digits for printed outputs */
    
    /* MLE-Based Model Parameters */
    nlmixedDF   = 10000,            /* Non-linear model degrees of freedom */
    nlmixedTech = ,                 /* Optimization algorithm (defaults to NEWRAP for one outcome, QUANEW for two) */
    
    /* Bayesian MCMC Parameters */
    chains      = 4,                /* Number of chains (default per Vehtari et al 2021) */
    GRthreshold = 1.01,             /* Gelman-Rubin convergence diagnostic threshold */
    seed        = 1235,             /* Random number generating seed */
    maxtune     = 50,               /* Maximum number of proposal tuning loops */
    ntu         = 1000,             /* Number of tuning iterations per MCMC proposal tuning phase */
    nbi         = 10000,            /* Number of burn-in MCMC iterations */
    nmc         = 50000,            /* Number of post-burn-in MCMC iterations */
    thin        = 1,                /* Thinning rate */
    accepttol   = ,                 /* Tolerance for target acceptance probabilities */
    targetaccept= ,                 /* Target acceptance rate for random walk Metropolis */
    propcov     = ,                 /* Method for constructing initial covariance matrix */
    init        = reinit,           /* Option for generating initial values (default: REINIT) */
    slicesampler= NO,               /* YES to use slice sampler instead of MH algorithm */
    checkSampleSize = YES,          /* YES to check sample size before proceeding */
    orpoly      = YES,              /* YES for orthogonal polynomial transformation */
    
    /* Prior Parameters - Segment 1 */
    malpha      = ,  palpha   = ,   /* Prior mean and precision for intercepts */
    mbeta1      = 0, pbeta1   = ,   /* Prior mean and precision for linear coefficients */
    mbeta2      = 0, pbeta2   = ,   /* Prior mean and precision for quadratic coefficients */
    mbeta3      = 0, pbeta3   = ,   /* Prior mean and precision for cubic coefficients */
    beta1l      = 0, beta1u   = ,   /* Bounds for Unif(a,b) prior for SD of linear coefficients */
    beta2l      = 0, beta2u   = ,   /* Bounds for Unif(a,b) prior for SD of quadratic coefficients */
    beta3l      = 0, beta3u   = ,   /* Bounds for Unif(a,b) prior for SD of cubic coefficients */
    
    /* Prior Parameters - Segment 2 [NEW in v2.4] */
    s2malpha    = ,  s2palpha = ,   /* Prior mean and precision for intercepts (segment 2) */
    s2mbeta1    = 0, s2pbeta1 = ,   /* Prior mean and precision for linear coefficients (segment 2) */
    s2mbeta2    = 0, s2pbeta2 = ,   /* Prior mean and precision for quadratic coefficients (segment 2) */
    s2mbeta3    = 0, s2pbeta3 = ,   /* Prior mean and precision for cubic coefficients (segment 2) */
    s2beta1l    = 0, s2beta1u = ,   /* Bounds for SD of linear coefficients (segment 2) */
    s2beta2l    = 0, s2beta2u = ,   /* Bounds for SD of quadratic coefficients (segment 2) */
    s2beta3l    = 0, s2beta3u = ,   /* Bounds for SD of cubic coefficients (segment 2) */
    
    /* AR(1) and Variance Prior Parameters */
    mrho        = 0, prho     = 1,  /* Prior mean and precision for transformed rho */
    taul        = 0.0001, tauu= ,   /* Bounds for Unif(a,b) prior for tau */
    vshape      = ,  vscale   = ,   /* Shape and scale for inverse gamma prior of variance */
    
    /* Dirichlet Prior Parameters */
    wdirichlet  = NO, wshape  = 2,  /* Whether to use Dirichlet prior for model weights */
    
    /* Advanced Options */
    mcmcplot    = NO,               /* YES to include within-chain trace/diagnostics plots */
    mcmclog     = NO,               /* YES to retain full posterior samples */
    saslognotes = NO,               /* YES to include SAS notes in log file */
    cmploc      = work.funcs        /* Location of CMP library */
) / minoperator;
```

## Parameter Categories

### 1. Printing Options

- **modelprint**: YES to print PROC NLMIXED and/or PROC MCMC results, NO otherwise (Default: NO)
- **finalprint**: YES to print MKF estimates for last time point, NO otherwise (Default: YES)

### 2. Data and Variables

**Note:** eMKF v2.4 streamlining - "format 2" (long/stacked format) is the only supported option

- **data**: Name of the dataset (must be in long/stacked format)
- **outcome**: The outcome of interest
- **se**: The standard error of the outcome
- **neff**: (Effective) sample sizes for random sampling variance case in Bayesian model (ignored in MLE-based models)
- **time**: The time variable
- **by**: Allows models to run for multiple strata simultaneously; can be used for simulations
- **group**: The different groups variable (e.g., race/ethnicity group)
- **breakPoint**: **[NEW in v2.4]** Specifies timepoint value at which trend break occurs (default is none)
- **breakType**: **[NEW in v2.4]** Indicates whether trend break impacts intercepts only (`level_break` or empty [default]), or all regression coefficients (`full_break`)

### 3. Modeling Option: The Linear Mixed Model

**slopes**: Type of slope needed. eMKF options:

- `indep_cubic`: Values of parameters b1, b2, and b3 computed for each group
- `indep_quad`: b3=0. Values of b1 and b2 computed for each group
- `indep_linear`: b3=0 and b2=0. Value of b1 computed for each group
- `common_cubic`: Values of b1, b2, and b3 assumed same across groups
- `common_quad`: b3=0. Values of b1 and b2 assumed same across groups
- `common_linear`: b3=0 and b2=0. Value of b1 assumed same across groups
- `dropped`: Model without time trend (b3=0, b2=0, b1=0)
- (empty): DEFAULT when using one outcome - results in Bayesian approach only

**Time point requirements:**
- Cubic trend requires at least 5 time points
- Quadratic trend requires at least 4 time points
- Linear trend requires at least 3 time points
- At least 2 time points required

**Model Averaging:**
- If multiple options selected, model averaging is computed
- When using two outcomes, DEFAULT is model average over all models up to cubic
- Combinations checked to ensure common reference model for Bayes factor calculations

### 4. Modeling Option: The Bayesian Model

**Bayesmodel**: Bayesian model to fit. eMKF options:

- `bma_cubic`: (v2.4 DEFAULT with no trend break) Bayesian model averaging using mixture of polynomial trends up to unconstrained cubic
- `bma_quad`: Bayesian model averaging using mixture of polynomial trends up to unconstrained quadratic
- `bma_linear`: (v2.4 DEFAULT with trend break) Bayesian model averaging using mixture of polynomial trends up to unconstrained linear
- `full_cubic`: Fully Bayesian cubic trends with multilevel priors
- `full_quad`: Fully Bayesian quadratic trends (b3=0) with multilevel priors
- `full_linear`: Fully Bayesian linear trends (b3=0, b2=0) with multilevel prior
- `indep_cubic`: Bayesian cubic trends with uninformative (flat) priors
- `indep_quad`: Bayesian quadratic trends (b3=0) with uninformative priors
- `indep_linear`: Bayesian linear trends (b3=0, b2=0) with uninformative prior
- `common_cubic`: Bayesian model with common cubic trend and uninformative priors
- `common_quad`: Bayesian model with common quadratic trend (b3=0) and uninformative priors
- `common_linear`: Bayesian model with common linear trend (b3=0, b2=0) and uninformative prior
- `dropped`: Bayesian model with no time trend (b3=0, b2=0, b1=0)

**BayesmodelAvg**: If YES (DEFAULT) and multiple models specified, model averaging up to the largest unconstrained polynomial trend is conducted

**ARmodel**: AR(1) model specification
- `common_ar` (default): Groups share AR(1) parameters rho and tausq
- `common_arh` **[NEW in v2.4]**: Each group has own variance parameter tausq, but shares common autocorrelation coefficient rho
- `indep_ar`: Each group has own set of stationary AR(1) parameters rho and tausq from common prior

**randomVars**: Default is YES. Sampling variances treated as scaled chi-squared random variables with inverse gamma prior

### 5. Output Options

- **out**: Output prefix name (Default: param)
  - Outputs saved as prefix + suffix:
    - `_pred`: Kalman prediction including original values and parameters
    - `_bayes`: Kalman prediction from Bayesian modeling
- **nlmixedDetails**: (v2.4 streamlining) Default NO - omits model datasets with fit details, retains only predictions and coefficients
- **xtrakeep**: Variables to keep in data while running models (e.g., weights, labels)
- **comparedata**, **comparedto**: Options for estimating disparities/differences in Bayesian setting

### 6. Advanced Parameters

**For experienced users only**

#### MLE-Based Model Parameters
- **nlmixedDF**: Non-linear model degrees of freedom (Default: 10,000)
- **nlmixedTech**: (v2.4 streamlining) Optimization algorithm (defaults to NEWRAP for one outcome, QUANEW for two)

#### Bayesian MCMC Parameters
- **chains**: Number of chains (Default: 4, per Vehtari et al 2021; DOI 10.1214/20-BA1221)
- **GRthreshold**: Gelman-Rubin convergence diagnostic threshold (Default: 1.01, per Vehtari et al)
- **mcmcplot**: YES to include within-chain trace/diagnostics plots (Default: NO)
- **mcmclog**: YES to retain full posterior samples in work directory (Default: NO - dataset will be large)
- **seed**: Random seed for reproducibility
- **maxtune**: Maximum proposal tuning loops (empty uses PROC MCMC default of 24; 0 skips tuning)
- **ntu**: Tuning iterations per MCMC proposal phase (empty uses PROC MCMC default of 500)
- **nbi**: Burn-in MCMC iterations (empty uses PROC MCMC default of 1000; 0 skips burn-in)
- **nmc**: Post-burn-in MCMC iterations (empty uses PROC MCMC default of 1000)
- **thin**: Thinning rate (empty uses PROC MCMC default of 1)
- **accepttol**: Tolerance for target acceptance probabilities (empty uses PROC MCMC defaults)
- **targetaccept**: Target acceptance rate for random walk Metropolis (empty uses PROC MCMC defaults)
- **propcov**: Method for constructing initial covariance matrix (empty uses PROC MCMC default of IND)
- **init**: Initial value generation option (eMKF default: REINIT - resets parameters after tuning)
- **slicesampler**: YES to use slice sampler instead of MH algorithm (Default: NO due to computational load)
- **checkSampleSize**: YES (default) to check adequate sample size before proceeding
- **orpoly**: YES (default) for orthogonal polynomial transformation using SAS IML orpol function

#### Prior Parameters

**NEW in eMKF v2.4:** User can optionally specify separate values for segment 2 in trend break scenarios

**Segment 1:**
- **malpha**, **palpha**: Prior mean and precision for intercepts
- **mbeta1**, **pbeta1**: Prior mean and precision for linear coefficient(s)
- **mbeta2**, **pbeta2**: Prior mean and precision for quadratic coefficient(s)
- **mbeta3**, **pbeta3**: Prior mean and precision for cubic coefficient(s)
- **beta1l**, **beta1u**: Bounds for Unif(a,b) prior for SD of linear coefficients (hyperprior for full_ models)
- **beta2l**, **beta2u**: Bounds for Unif(a,b) prior for SD of quadratic coefficients (hyperprior for full_cubic/full_quad)
- **beta3l**, **beta3u**: Bounds for Unif(a,b) prior for SD of cubic coefficients (hyperprior for full_cubic)

**Segment 2 (prefix with s2):**
- **s2malpha**, **s2palpha**: Prior mean and precision for intercepts
- **s2mbeta1**, **s2pbeta1**: Prior mean and precision for linear coefficient(s)
- **s2mbeta2**, **s2pbeta2**: Prior mean and precision for quadratic coefficient(s)
- **s2mbeta3**, **s2pbeta3**: Prior mean and precision for cubic coefficient(s)
- **s2beta1l**, **s2beta1u**: Bounds for SD of linear coefficients
- **s2beta2l**, **s2beta2u**: Bounds for SD of quadratic coefficients
- **s2beta3l**, **s2beta3u**: Bounds for SD of cubic coefficients

**AR(1) and Variance:**
- **mrho**, **prho**: Prior mean and precision for transformed rho (psi = -ln[(1-rho)/(1+rho)])
- **taul**, **tauu**: Bounds for Unif(a,b) prior for tau (SD of innovation variance tausq)
- **vshape**, **vscale**: Shape and scale parameters for inverse gamma prior when randomVars = YES

**Dirichlet Prior:**
- **wdirichlet**: (v2.4 streamlining) Whether to use Dirichlet prior for model weights (Default: NO - not fully tested)
- **wshape**: Common shape parameter for Dirichlet prior on model indicators (Default: 2)

#### System Options
- **saslognotes**: **[NEW in v2.4]** YES or empty to include SAS notes in log; NO to turn off (warnings/errors kept; Default: NO)
- **cmploc**: Location of CMP library (sasuser.funcs or work.funcs if sasuser is write-protected)

## Important Notes

- If parameters are empty, data will be used to inform starting values
- Parameters should only be modified by experienced users unless specified
- v2.4 streamlining: _rho_ and _tausq_ can no longer be user-supplied - they are estimated instead
- Regression coefficients will be reverse-transformed if orpoly=YES
- Prior parameters are assumed to be for orthogonal polynomial regression coefficients if orpoly=YES
