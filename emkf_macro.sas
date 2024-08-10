/*
 * Version 1.4 10-Aug-2024
 *
 * eMKF: Expansion of RAND's MKF macro to:
 *
 *  - allow for unequally-spaced time points.
 *	- allow for quadratic and cubic trends in both MLE-based and Bayesian estimation settings.
 *  - allow for model averaging over (orthogonal) polynomial trend model sequences up to cubic in both MLE-based and Bayesian settings.
 *  - allow for common as well as group-specific AR parameters rho and tausq for the random effects in the Bayesian setting.
 * 	- allow for random sampling variances in the Bayesian setting.
 *  - implement Gibbs sampling in PROC MCMC using user defined samplers (UDSs) that are precompiled with PROC FCMP, replacing .exe file (C code).
 *  - expand calculations of between-group disparities at the latest time point in the Bayesian setting.
 *   
 * See README.md file for additional details on methodological differences between this eMKF version and RAND's MKF.
 *
 * Comments and code inserted for the eMKF version in this file are prefixed with *eMKF or otherwise indicated. 
 * All other comments and code are from the developers of the original MKF macro.
 *
 * Makram Talih, Ph.D. (NCHS contractor)
 * 
 * Macros defined in this file:
 *
 * - MKF (main) .............................................................. line    41
 * - bayesfit (Bayesian estimation workhorse) ................................ line  3645
 * - bayesBMA (Bayesian model averaging workhorse, via mixture prior) ........ line  4925
 * - htrp (ML-based estimation workhorse) .................................... line  6205
 * - htrp2d (ML-based estimation workhorse for 2 outcomes) ................... line  7697
 * - reformat (set up of dataset for analysis) ............................... line  9664
 * - _counts_, zeros, thevacompr, etc. (various utility macros) .............. line 10082
 * - gibbs_uds_compile_EP (Gibbs sampler for true state predictions) ......... line 10245
 * - gibbs_uds_compile_RP (Gibbs sampler for random sampling variances) ...... line 10327
 * - gibbs_uds_compile_MP (Gibbs samplers for mean hyper-parameters) ......... line 10364
 * - gibbs_uds_compile_CP (Gibbs samplers for regression coefficients) ....... line 10540
 * - gibbs_uds_compile_FP (Gibbs samplers for model flags) ................... line 12689
 *
 */

data _null_;
run;

/* MKF macro -- eMKF expansion in 2023 Q1-Q2 by M. Talih
Macro defined 03-16-2007, Revised 05-08-09 D. McCaffrey. Revised 03-15-2010
It allows for the estimation of parameters from the non-linear mixed effect model.
In addition it allows for the estimation of model averaging as well as bayesian estimates

1- Printing (modified for eMKF)

modelprint  : YES will print the proc nlmixed and/or proc mcmc results and NO will not. Default is NO.
finalprint  : YES will print MKF estimates for last time point at the end of the estimations and NO will not. Default is YES.

2- data and variables (format 2 strongly advised for eMKF, as eMKF was not fully tested with format 1)

data 	: name of the dataset.
		It can be a data where each year outcome is listed differently, i.e., Y1-Y8 (format 1)
        or it can be an outcome where one year is stacked up on top of the other (format 2) and in that case
        the time variable will need to be specified
outcome : the outcome of interest. It can be in two format: 
          a- multiple outcomes can be specified Y1-Y8 which means the data is in the format (1)
          b- a single outcome with data format (2) and the time variable will need to be specified
se      : the standard errors of the outcomes. It also needs to be in the a or b format to match the outcome format
neff	: (eMKF) (Effective) sample sizes of the outcomes for the random sampling variance case in the Bayesian model
		  If empty, but random sampling variances are requested, an error is returned.
		  Ignored in the MLE-based models.
time    : the time variable for data format (2). If missing, will assume format (1)
by      : allows models to run for multiple strata at the same time and can be used for simulations
group   : the different groups (e.g., race/ethnicity group) variable

3- Modeling option: The linear mixed model

slopes : the type of slope one needs: eMKF options are the following:

         indep_cubic	: The values of the parameters b1, b2, and b3 are computed for each group
         indep_quad		: b3=0. The values of the parameters b1 and b2 are computed for each group
         indep_linear   : b3=0 and b2=0. The value of the slope b1 is computed for each group
         common_cubic	: The values of each of the parameters b1, b2, and b3 are assumed to be the same across groups
         common_quad	: b3=0. The values of each of the parameters b1 and b2 are assumed to be the same across groups
         common_linear  : b3=0 and b2=0. The value of the slope b1 is assumed to be the same across groups
         dropped    	: A model without time trend is computed (b3=0, b2=0, b1=0).
         (empty)		: DEFAULT when using one outcome -- results in Bayesian approach only.

         eMKF: the number of time points will be checked and an error returned if:
	  	  * Cubic trend is requested when there are less than 5 time points
	  	  * Quadratic trend is requested when there are less than 4 time points
	  	  * Linear trend is requested when there are less than 3 time points
	  	  * There are less than 2 time points

         If multiple options are selected, model averaging is computed using those options
           * eMKF: When using two outcomes, DEFAULT is model average over all models up to cubic.
           * eMKF: Combinations will be checked to ensure a common reference (shared descendent) model for the Bayes factor calculations.
		  	** If both indep_quad and common_cubic are selected, then common_quad will be added if needed.
		  	** If both indep_linear and common_cubic are selected, then common_linear will be added if needed.
		  	** If both indep_linear and common_quad are selected, then common_linear will be added if needed.

4- Modeling option: The Bayesian model

Bayesmodel : This is the Bayesian model one wants to fit. eMKF options are the following:

		bma_cubic     : (DEFAULT) Bayesian model averaging implemented using a mixture of polynomial trend models up to unconstrained cubic.
		bma_quad	  : Bayesian model averaging implemented using a mixture of polynomial trend models up to unconstrained quadratic.
		bma_linear    : Bayesian model averaging implemented using a mixture of polynomial trend models up to unconstrained linear.
        full_cubic 	  : Fully Bayesian cubic trends. The b1s, b2s, and b3s are independent draws from multilevel priors.
        full_quad 	  : Fully Bayesian quadratic trends. Here, b3s are 0, and b1s and b2s are independent draws from multilevel priors.
        full_linear   : Fully Bayesian linear trends. Here, b3s and b2s are 0, and b1s are indep. draws from a multilevel prior.
        indep_cubic	  : Bayesian cubic trends. The b1s, b2s, and b3s are independent draws from uninformative (flat) priors.
        indep_quad	  : Bayesian quadratic trends. Here, b3s are 0 and b1s and b2s are independent draws from uninformative (flat) priors.
        indep_linear  : Bayesian linear trends. Here, b3s and b2s are 0, and b1s are independent draws from an uninformative (flat) prior.
        common_cubic  : Bayesian model with common cubic trend. The common b1, b2, and b3 are each drawn from an uninformative (flat) prior.
        common_quad	  : Bayesian model with common quadratic trend. b3=0 and the common b1 and b2 are each drawn from an uninformative (flat) prior.
        common_linear : Bayesian model with common linear trend. b3=0, b2=0, and the common b1 is drawn from an uninformative (flat) prior.
        dropped    	  : Bayesian model with no time trend is computed (b3=0, b2=0, b1=0).

BayesmodelAvg :	(eMKF) If multiple options are selected in Bayesmodel and BayesmodelAvg = NO or empty, 
				all selected models are included in the output datasets. 
			 	However, only the last model in the specified sequence is used for printing the MKF estimates for last time point.
 			 	If BayesmodelAvg = YES (DEFAULT) and more than one model is specified, then model averaging up to largest 
				unconstrained polynomial trend that was listed is conducted.

ARmodel : (eMKF) Default is common_ar. 
		  If indep_ar, each group will have its own set of AR parameters rho and tausq, with common prior distribution. 

randomVars : (eMKF) Default is YES. 
			 Sampling variances will be treated as scaled chi-squared random variables with an inverse gamma prior. 

5- Output prefix (modified for eMKF)

out     : The name for the output prefix. Default is param. 
		  Outputs (prefix + suffix) are saved. Here are some of the relevant suffixes:
           _pred    : Kalman prediction of the outcome of interest includes original values as well as parameters
           _bayes   : Kalman prediction of the outcome from bayesian modeling
         (e.g., for OUT=param then PARAM_PRED will be the Kalman prediction data of the outcome of interest.)

xtrakeep: Any variable one wants to keep in the data while running models: weights, ... (eMKF: could be used to retain labels for multiyear data)

comparedata, comparedto  : (eMKF) options from MKF allowing for estimating disparities/differences in Bayesian setting

6- Parameters that should only be changed by experienced users (modified for eMKF)

 pdigit	: Number of decimal digits for the printed outputs. Default is set to 4.

 MLE-based model parameters (eMKF: passed on to proc nlmixed)

 _rho_   : the value of the true rho that generated the data, if known. If not given, it will be estimated.
 _tausq_ : the value of the true tau-square that generated the data, if known. If not given, it will be estimated.
 DF      : Non-linear model degrees of freedom. The default is set pretty high at 1,0000

 Bayes model parameters (eMKF: passed on to proc mcmc)

 chains   			: number of chains to run for the Bayesian estimation. Default is 4, as per Vehtari et al (2021; DOI 10.1214/20-BA1221)
 GRthreshold		: threshold to use for Gelman-Rubin convergence diagnostic, usually no larger than 1.1. Default is 1.01, as per Vehtari et al.
 mcmcplot			: if YES, within-chain trace/diagnostics plots from proc mcmc will be included (default is NO)
 mcmclog   			: if YES, the full posterior samples will be retained in the work directory (default is NO, as dataset will be large)
 seed     			: random number generating seed that will allow the user to reproduce the same results in the Bayesian model
 maxtune			: maximum number of proposal tuning loops (if empty, proc mcmc default of 24 is used; if 0, tuning will be skipped)
 ntu				: number of tuning iterations to use in each MCMC proposal tuning phase (if empty, proc mcmc default of 500 is used)
 nbi	 			: number of burn-in in MCMC iterations (if empty, proc mcmc default of 1000 is used; if 0, burn-in will be skipped)
 nmc	 			: number of post-burn-in MCMC iterations (if empty, proc mcmc default of 1000 is used)
 thin				: controls thinning rate (if empty, proc mcmc default of 1 is used)
 accepttol			: tolerance for target acceptance probabilities (targetaccept +|- accepttol). If empty, proc mcmc defaults are used.
 targetaccept		: target acceptance rate for random walk Metropolis. If empty, proc mcmc defaults are used (see proc mcmc documentation). 
 propcov			: method used in constructing initial covariance matrix for the MH algorithm (see proc mcmc documentation).
					  If empty, proc mcmc default of IND will be used.
 init				: Option for generating initial values for the parameters. If empty, proc mcmc defaults to MODE (see proc mcmc documentation).
					  eMKF default is REINIT, which resets model parameters to the user-supplied initial values after tuning. 
 slicesampler		: YES to use the slice sampler instead of MH algorithm for parameters that are not included in the Gibbs sampling steps. 
					  Default is NO due to heavier computational load.
 checkSampleSize	: YES (default) to check sample size is large enough before proceeding. Set to NO to replicate earlier MKF from RAND.
 orpoly  			: YES (default) for pre-transforming the design matrix using SAS IML orpol function. NO for "raw" polynomials.
          			  If YES, regression coefficients will be reverse-transformed prior to macro end. 
					  However, prior parameters below are assumed to be for the coefficients of the orthogonal polynomial regression if orpoly=YES.
 malpha , palpha 	: prior mean and precision for alphas
 mbeta1				: prior mean for mean linear coefficient across groups -- used for cubic, quadratic, or linear trends (full_, indep_, common_)
 pbeta1 			: prior precision for mean linear coefficient across groups -- used for SD hyperprior(s) in all 3 full_ models
 mbeta2 			: prior mean for mean quadratic coefficient across groups -- used for cubic or quadratic trends (full_, indep_, common_)
 pbeta2			 	: prior precision for mean quadratic coefficient across groups -- used for SD hyperprior(s) in full_cubic or full_quad
 mbeta3 			: prior mean for mean cubic coefficient across groups -- used for cubic trends (full_, indep_, and common_)
 pbeta3			 	: prior precision for mean cubic coefficient across groups -- used for SD hyperprior(s) in full_cubic
 beta1l, beta1u		: bounds for U(a,b) prior for SD of linear coefficients across groups -- used for hyperprior(s) in all 3 full_ models
 beta2l, beta2u		: bounds for U(a,b) prior for SD of quadratic coefficients across groups -- used for hyperprior(s) in full_cubic or full_quad
 beta3l, beta3u		: bounds for U(a,b) prior for SD of cubic coefficients across groups -- used for hyperprior(s) in full_cubic
 mrho, prho			: prior mean and precision for transformed rho -- ie., psi = ln[(1-rho)/(1+rho)]
 taul, tauu			: bounds for U(a,b) prior for tau (SD of innovation variance tausq)
 vshape				: Shape parameter for inverse gamma prior distribution of the variance when randomVars = YES
 vscale				: Scale parameter for inverse gamma prior distribution of the variance when randomVars = YES
 wshape				: Common shape parameter to use for Dirichlet prior on model indicators in mixture priors (default = 2)
 cmploc				: Desired location of CMP library (sasuser.funcs or work.funcs if sasuser is write-protected)
*/

%macro mkf(  
             data			= , 
             outcome		= , 
             se				= , 
			 neff 			= ,
			 outcome2		= , 
             se2			= , 
			 neff2 			= ,
		     by				= ,  
             group			= , 
             time			= ,
             slopes			= ,
			 Bayesmodel		= bma_cubic,
			 BayesmodelAvg 	= YES,
			 randomVars		= YES,
			 ARmodel		= common_ar,
             xtrakeep		= ,
             out			= param, 
			 comparedto		= ,
			 comparedata	= ,
             modelprint 	= NO,
             finalprint 	= YES,
  		/* eMKF: The parameters from here down should not be modified unless the user really knows what s/he is doing*/
			 pdigit 		= 4,
			 _rho_			= , 
             _tausq_		= , 
             DF 			= 10000, 			 
			 chains  		= 4 ,
			 GRthreshold 	= 1.01,
			 seed			= 1235,
			 maxtune		= 50,
			 ntu			= 1000,
			 nbi			= 10000,
			 nmc			= 50000,
			 thin			= 1,			
			 accepttol		= ,	 
			 targetaccept	= ,
			 propcov  		= ,
			 init 			= reinit,
			 slicesampler 	= NO,
			 checkSampleSize= YES,
			 orpoly 		= YES,
  		/* eMKF: If parameters below are empty, the data will be used to inform the starting values */
			 malpha 		= , 
			 palpha 		= , 
			 mbeta1  		= 0,
			 pbeta1  		= , 
			 mbeta2  		= 0,
			 pbeta2  		= , 
			 mbeta3  		= 0,
			 pbeta3  		= , 
			 beta1l  		= 0,
			 beta1u  		= , 
			 beta2l  		= 0,
			 beta2u  		= ,
			 beta3l  		= 0,
			 beta3u  		= , 
			 mrho   		= 0,
			 prho   		= 1,
			 taul   		= 0.0001,
			 tauu   		= ,
			 vshape 		= ,
			 vscale 		= ,
			 wshape 		= 2,
			 mcmcplot 		= NO,
			 mcmclog 		= NO,
			 cmploc 		= work.funcs
            );

%local _oo1_ _oo2_ ui uii uj uk un um uvar newuvar ug uloc
	   flag1 flag2 flag3 flag4 flag5 flag6 flag7 flag1f flag2f flag3f flag1a flag2a flag3a 
	   crep run1 run2 run3 Bayesian _chainseed 
       xtrakeep22 toprint toprint2 _ssby _ssn _thekeeps _thekeepsb _thekeep1 _thekeep1b _thekeep2 _thekeep3 _thet 
	   emkfkeep emkfrename _BMAmodel _slopes _rtimess _igrp_ _t_ _comp2;

/* eMKF: Added error check for the length of &out prefix */
%if &out = %str() or %length(&out) > 16 %then %do;
	%put ERROR: Please specify a non-empty prefix (out=&out) that is no longer than 16 characters in length;
	proc iml;
		print "  Error Note:";
		print "  Please specify a non-empty prefix (out=&out) that is no longer than 16 characters in length";
	quit;
	%return;
%end;

/* eMKF: Added error check for the length of &outcome variable name */
%if &outcome = %str() or %length(&outcome) > 18 %then %do;
	%put ERROR: Please specify a non-empty variable name for outcome (&outcome) that is no longer than 18 characters in length;
	proc iml;
		print "  Error Note:";
		print "  Please specify a non-empty variable name for outcome (&outcome) that is no longer than 18 characters in length";
	quit;
	%return;
%end;

/* eMKF: Added error check for the length of &outcome2 variable name (if applicable) */
%if &outcome2 ^= %str() and %length(&outcome2) > 18 %then %do;
	%put ERROR: Please specify a variable name for outcome2 (&outcome2) that is no longer than 18 characters in length;
	proc iml;
		print "  Error Note:";
		print "  Please specify a variable name for outcome2 (&outcome2) that is no longer than 18 characters in length";
	quit;
	%return;
%end;

/* eMKF: Added error check if both slopes and bayesmodel are empty */
%if &Bayesmodel = %str() and &slopes = %str() %then %do;
	%put ERROR: There appears to be no estimation method specified: Please check!;
	proc iml;
		print "  Error Note:";
		print "  There appears to be no estimation method specified. Please check. ";
	quit;
	%return;
%end;

/*eMKF: Check for any mispellings/misspecification of the model(s) */
%if &slopes ^= %str() %then %do;
	%let uvar=; %let ui=0;
    %do ui=1 %to %_counts_(&slopes);
		%let uvar=;
		%let uvar=%scan(&slopes, &ui);
		%if %upcase(&uvar) ^=INDEP_CUBIC and %upcase(&uvar) ^=INDEP_QUAD and %upcase(&uvar) ^=INDEP_LINEAR and 
    		%upcase(&uvar) ^=COMMON_CUBIC and %upcase(&uvar) ^=COMMON_QUAD and %upcase(&uvar) ^=COMMON_LINEAR and %upcase(&uvar) ^=DROPPED %then %do;
				%put ERROR: &uvar is not a supported ML model specification; 
                %put ERROR- Model(s) must be one (or more) of indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, or dropped;
				proc iml;
					print "  Error Note:";
					print "  Specified ML model(s) not supported.";
                    print "  Model(s) must be one (or more) of indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, or dropped. ";
				quit;
				%return;
		%end;
	%end;
%end;

/*eMKF: Bayesian modeling flag and check for any mispellings/misspecification of the model(s) */
%let Bayesian=;
%if &Bayesmodel ^=%str() %then %do;
	%let Bayesian= Yes;
	%let uvar=; %let ui=0;
    %do ui=1 %to %_counts_(&Bayesmodel);
		%let uvar=;
		%let uvar=%scan(&Bayesmodel, &ui);
		%if %upcase(&uvar) ^=BMA_CUBIC and %upcase(&uvar) ^=BMA_QUAD and %upcase(&uvar) ^=BMA_LINEAR and 
			%upcase(&uvar) ^=FULL_CUBIC and %upcase(&uvar) ^=FULL_QUAD and %upcase(&uvar) ^=FULL_LINEAR and 
			%upcase(&uvar) ^=INDEP_CUBIC and %upcase(&uvar) ^=INDEP_QUAD and %upcase(&uvar) ^=INDEP_LINEAR and 
    		%upcase(&uvar) ^=COMMON_CUBIC and %upcase(&uvar) ^=COMMON_QUAD and %upcase(&uvar) ^=COMMON_LINEAR and %upcase(&uvar) ^=DROPPED %then %do;
			%put ERROR: &uvar is not a supported Bayesian model specification;
            %put ERROR- Model(s) must be one (or more) of bma_cubic, bma_quad, bma_linear, full_cubic, full_quad, full_linear, indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, or dropped;
			proc iml;
				print "  Error Note:";
				print "  Specified Bayesian model(s) not supported."; 
                print "  Model(s) must be one (or more) of bma_cubic, bma_quad, bma_linear, full_cubic, full_quad, full_linear, indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, or dropped. ";
			quit;
			%return;
		%end;
	%end;
%end;
%else %let Bayesian = No; /* eMKF: Added to allow for Bayes approach to be explicitly disabled */

/*eMKF: Set up internal &by variable name _reps, e.g., for use in calls to HTRP and HTRP2D */
%let _ssby=;
%if &by ^=%str() %then %let _ssby=_reps;

/*eMKF: Set up dataset name for differences/disparities if requested */
%if &comparedata = %str() and &comparedto ^= %str() %then %let comparedata = &out._diff;

/*eMKF: Set up internal names for the outcome(s) variable(s) */
%let _oo1_=; %let _oo2_=;
%if &outcome  ^=%str() %then %let _oo1_=%scan(&outcome, 1);
%if &outcome2 ^=%str() %then %let _oo2_=%scan(&outcome2, 1);

/*eMKF: Initialize macro variables for use in symbolic calculations/operations */
%let _thet=; %let _thekeeps=; %let _thekeepsb=; %let _thekeep1=; %let _thekeep1b=;
%let _thekeep2= &by &group &time; 
%let _thekeep3=;
%if &time = %str() %then %let _thekeep2 = &_thekeep2 _time; /*eMKF: this was introduced in MKF to deal with data in format 1*/
%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;
	%if %scan(&outcome,2) ^=%str() %then %let _thekeep2 = &_thekeep2 _y _se;
	%if %scan(&outcome,2)  =%str() %then %let _thekeep2 = &_thekeep2 &outcome &se;
	%if &by ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &by);
	%let _thekeep3= &_thekeep3 &out.(keep= &group);
	%if &time ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &time); /*eMKF: format 2 */
	%if &time  =%str() %then %let _thekeep3= &_thekeep3 &out.(keep= _time); /*eMKF: format 1 */
	%if %scan(&outcome,2) ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= _y) &out.(keep= _se);
	%if %scan(&outcome,2)  =%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &outcome) &out.(keep= &se);
%end;
%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	%if %scan(&outcome,2) ^=%str() %then %let _thekeep2 = &_thekeep2 _y _se _y2 _se2;  
	%if %scan(&outcome,2)  =%str() %then %let _thekeep2 = &_thekeep2 &outcome &se &outcome2 &se2;
	%if &by ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &by);
	%let _thekeep3= &_thekeep3 &out.(keep= &group);
	%if &time ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &time);
	%if &time  =%str() %then %let _thekeep3= &_thekeep3 &out.(keep= _time);
	%if %scan(&outcome,2) ^=%str() %then %let _thekeep3= &_thekeep3 &out.(keep= _y) &out.(keep= _se) &out.(keep= _y2) &out.(keep= _se2);
	%if %scan(&outcome,2)  =%str() %then %let _thekeep3= &_thekeep3 &out.(keep= &outcome) &out.(keep= &se) &out.(keep= &outcome2) &out.(keep= &se2);
%end;

/* eMKF: Bayesian estimation not yet implemented for two outcomes--use MLE-based estimation instead */
%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	%let Bayesian= No; 
	%let Bayesmodel= ;
	/* eMKF: If slopes parameter is unspecified, default to MA up to cubic */
	%if &slopes=%str() %then %let slopes= indep_cubic common_cubic indep_quad common_quad indep_linear common_linear dropped;
%end; 

/******************************************************/
/* eMKF: Reformat the data only once at the top level */
/******************************************************/
data _bayesdata_ ;
run;

/*eMKF: macro variable &_ssn will be used in call to bayesfit to indicate generic variable name for neff */
%let _ssn = ;

/*eMKF: reformat data if needed and calculate variables needed for processing */
%if %upcase(&Bayesian) = YES and %upcase(&randomVars) = YES %then %do;
	%if &neff = %str() %then %do; /*eMKF: if no effective sample sizes were provided, return an error */
		%put ERROR: (Effective) sample sizes neff must be specified to fit random sampling variances;
		proc iml;
			print "  Error Note:";
			print "  (Effective) sample sizes neff must be specified to fit random sampling variances. ";
		quit;
		%return;
	%end;
	%reformat(data=&data, outcome=&outcome, se=&se, neff=&neff, outcome2=&outcome2, se2=&se2, neff2=&neff2, 
			  group=&group, time=&time, by=&by, randomVars=YES, outformat= _bayesdata_ );
	%let _ssn = _n; 
%end;
%else %do;
	/*eMKF: Ignore effective samples sizes if random variances not requested in Bayesian setting, or for MLE-based estimation */
	%reformat(data=&data, outcome=&outcome, se=&se, outcome2=&outcome2, se2=&se2,
			  group=&group, time=&time, by=&by, randomVars=NO, outformat= _bayesdata_ );
	%let randomVars = NO;
%end;

proc sort data= _bayesdata_;
  by &by _group_ _time ;
run;

/*eMKF: Note that _nlmixdata_ (formerly _bayesdata00_) has renamed variables */
/*eMKF: Added columns with new variables instead of re-naming, as original MKF code was re-creating generic variables _y, _time, etc. */

/*eMKF: Create copies of _group_ and _rep variables for use in proc iml matrix calculations */
%if %upcase(&Bayesian) ^= YES %then %do;
    data _nlmixdata_; 
      set _bayesdata_; 
      _groupnum = _group_; 
      %if &by ^=%str() %then _reps	= _rep;; 
    run;
%end;

/* eMKF: Datasets to check number of time points and groups are consistent within and across strata */
data _junk_ _freq1_ _freq2_ _freq3_;
run;

/* eMKF: Original RAND macro assumed < 10 groups and timepoints at this point in the code. Here, we do away with this assumption. */
data _junk_;
  set _bayesdata_;
  if _y=. then delete;
  if _se=. then delete;
  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _y2=.  then delete;;;
  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _se2=. then delete;;;
run;

proc freq data=_junk_ noprint;
  tables _group_ /list out=_freq1_(rename=(count=ntime));
  %if &by ^=%str() %then by &by;;
run;
proc freq data=_junk_ noprint;
  tables _time /list out=_freq2_(rename=(count=ngroup));
  %if &by ^=%str() %then by &by;;
run;

/*eMKF: added if-clause to avoid creating dataset _finalprint_ if it remains unused */
%if %upcase(&finalprint) = YES %then %do; 

	options formchar="|----|+|---+=|-/\<>*"; /*eMKF: system options for correct SAS monospace font in printed tables*/

	data _finalprint_;
	run;
	data _finalprint_;
	  set _freq2_;
	  _rp_=1;
	run;
	data _finalprint_;
	  set _finalprint_;
	  keep=0;
	  %if &by ^=%str() %then by &by _time;;
	  %if &by  =%str() %then by _rp_ _time;;
	  %if &by ^=%str() %then if last.&by and last._time then keep=1;;;
	  %if &by  =%str() %then if last._rp_ and last._time then keep=1;;;
	  drop _rp_;
	run;
	data _finalprint_;
	  set _finalprint_;
	  if keep=1;
	  keep &by _time;
	run;

%end;

proc sort data=_freq1_ nodupkey;
  by &by ntime;
run;
proc sort data=_freq2_ nodupkey;
  by &by ngroup;
run;
data _freq1_;
  merge _freq1_ _freq2_;
  %if &by ^=%str() %then by &by;;
  drop _group_ _time percent;
run;

data _freq1_;
  set _freq1_;
  %if &by ^=%str() %then by &by;;
  %if &by ^=%str() %then if first.&by then id =0;;;
  stop=0;
  id +1;
run;

data _freq2_;
run;
data _freq2_;
  set _freq1_;
  if id=1;
run;
data _freq1_;
  set _freq1_;
  stop=1;
  if id=2;
run;

proc sort data= _freq1_ nodupkey;
  by &by id;
run;
data _freq2_;
  merge _freq2_ _freq1_(drop=id);
  %if &by ^=%str() %then by &by;;
run;

%let run1=0; %let run2=0;

data _null_;
  set _freq2_;
  one=1;
  if stop = 1 then call symput("run1" , stop);
  if stop ne 1 then call symput("run2" , one);
run;

%let run2= %eval(&run2 +0);
%let run1= %eval(&run1 +0);

/* eMKF: Precaution to ensure groups and timepoints are consistent across &by strata */
%if &by ^= %str() and &run1 = 0 %then %do;

  %let run3 = 0;
  proc freq data=_freq2_ noprint;
    tables ngroup /list out=_freq3_;
  run;
  data _freq3_;
	set _freq3_;
	stop = 0;
	_ng_ + 1;
    if _ng_ > 1 then stop = 1;
	call symput("run3", stop);
  run;
  %let run3= %eval(&run3 +0);

  %if &run3 = 0 %then %do;
    proc freq data=_freq2_ noprint;
      tables ntime /list out=_freq3_;
    run;
    data _freq3_;
	  set _freq3_;
	  stop = 0;
	  _nt_ + 1;
	  if _nt_ > 1 then stop = 1;
	  call symput("run3", stop);
    run;
  %end;

  %let run1= %eval(&run3 +0);
%end;

proc datasets nolist;
  delete _junk_ _freq1_ _freq2_ _freq3_;
run ;

%if &run1=0 and &run2=0 %then %do;
	%put ERROR: There appears to be no data to work with: Please check!;
	proc iml;
	  print "  Error Note:";
	  print "  There appears to be no data to work with. Please check. ";
	quit;
	%return; /*eMKF: added return functionality for easier debugging*/
%end;

%if &run1=1 and &run2=0 %then %do;
	%put ERROR: Please check your data: The number of valid timepoints is inconsistent across groups;
	proc iml;
	  print "  Error Note:";
	  print "  An error occurred with your data. ";
	  print "  Check the data to make sure there are no missing values for means/SEs or only zero SEs,  ";
	  print "  and that all the groups have data for the same number of timepoints (&time).  ";
	  print "  You may need to combine some groups and/or timepoints to avoid this problem.  ";
	quit;
	%return; /*eMKF: added return functionality for easier debugging*/
%end;

%if &run1=1 and &run2=1 %then %do;
	%put ERROR: Please check your data: The number of valid timepoints is inconsistent across groups and strata;
	%put ERROR- This problem occurred in at least one subgroup;
	proc iml;
	  print "  Error Note:";
	  print "  An error occurred with your data. ";
	  print "  Check the data to make sure there are no missing values for means/SEs or only zero SEs,  ";
	  print "  and that all the groups and strata have data for the same number of timepoints (&time).  ";
	  print "  You may need to combine some groups and/or timepoints to avoid this problem.  ";
	quit;
	%return; /*eMKF: added return functionality for easier debugging*/
%end;

%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	%let ug=0;
	data _null_;
	  set _bayesdata_;
	  one=1;
	  if _se2 ne . then call symput('ug',one);
	run;
	%let ug = %eval(&ug + 0);
	%if &ug =0 %then %let run1=1;
	%if &ug =0 %then %do;
		%put ERROR: Please check your data: You are using two outcomes &outcome and &outcome2;
		%put ERROR- Those outcomes or their standard errors should be of the same format: Check that and rerun the model!;
		proc iml;
		  print "  Error Note:";
		  print "  An error occurred with your data. ";
		  print "  Check the data to make sure that both outcomes and SEs are in the same format.  ";
		quit;
		%return; /*eMKF: added return functionality for easier debugging*/
	%end;
%end;

%if &run1=0 and &run2=1 %then %do;

	/* Setting up the number of time points and the number of groups */

	%let ug=0; /* eMKF: Number of groups */
	data _freqg_;
	run;
	proc freq data=&data noprint;
	  tables &group /list out=_freqg_;
	run;
	data _freqg_;
	  set _freqg_;
	  _group_ +1;
	  call symput('ug',_group_);
	  keep _group_ &group;
	run;
	%let ug= %eval(0 + &ug);

	%let un=0;  /*eMKF: Macro variable for the number of time points */
	data _freqn_;
	run;
	proc freq data=_bayesdata_ noprint;
	  tables _rtime /list out=_freqn_;
	run;
	data _freqn_;
	  set _freqn_;
	  _time+1;
	  call symput('un',_time);
	  keep _rtime _time;
	run;
	%let un= %eval(0 + &un);

	%let _rtimess = ; /*eMKF: Macro variable for the real times to use in calculations */
	data _freqn_;
	  set _freqn_;
	  retain _rts;
	  if _n_= 1 then _rts = cat(_rtime);
	  else _rts = catx(" ", _rts, _rtime);
	  call symput('_rtimess', _rts);
	  drop _rts;
	run;

	proc datasets nolist;
  		delete _freqn_ _freqg_ ;
	run ;

	%if &slopes ^=%str() %then %do;

		/* eMKF: Guard against too many groups/timepoints: code for current implementation cannot handle more than 204 groups or more than 5508 data points */
		%if (&ug > 204 or %eval(&ug*&un) > 5508) %then %do;
			%put ERROR: Code for current implementation of this macro cannot handle more than 204 groups or more than 5508 data points (group by time) per stratum;
			%put ERROR- Please reduce the number of groups and/or timepoints;
			proc iml;
			  print "  Error Note:";
			  print "  Code for current implementation of this macro cannot handle more than 204 groups or more than 5508 data points (group by time) per stratum. ";
			  print "  Please reduce the number of groups and/or timepoints. ";
			quit;
			%return;
		%end;

		/* eMKF: Make sure there are enough time points to fit the desired trends (up to cubic) */

		%if (&un < 2) %then %do; /*eMKF: minimal check for 2+ points */
			  %put ERROR: There are not enough timepoints to work with! At least 2 timepoints are required (4+ recommended) to use this macro;
			  proc iml;
			    print "  Error Note:";
			    print "  There are not enough timepoints to work with! At least 2 timepoints are required (4+ recommended) to use this macro. ";
			  quit;
			  %return;
		%end;

		%if %upcase(&checkSampleSize) = YES %then %do;

			%if (&un < 4) %then %do; /*eMKF: added check for 4+ points */
				  %put ERROR: There are not enough timepoints to work with! At least 4 timepoints are recommended to use this macro;
				  proc iml;
				    print "  Error Note:";
				    print "  There are not enough timepoints to work with! At least 4 timepoints are recommended to use this macro. ";
				  quit;
				  %return;
			%end;

			%let ui=0;
	        %do ui=1 %to %_counts_(&slopes);
				%let uvar=;
				%let uvar=%scan(&slopes, &ui);

				%if (%upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR) and (&un < 5) %then %do;
					  %put ERROR: At least 5 timepoints are recommended to fit linear trends: You may try dropping the time trends instead;
					  proc iml;
					    print "  Error Note:";
					    print "  At least 5 timepoints are recommended to fit linear trends. You may try dropping the time trends instead. ";
					  quit;
					  %return;
				%end;
				%if (%upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD) and (&un < 6) %then %do;
					  %put ERROR: At least 6 timepoints are recommended to fit quadratic trends: You may try linear trends instead;
					  proc iml;
					    print "  Error Note:";
					    print "  At least 6 timepoints are recommended to fit quadratic trends. You may try linear trends instead. ";
					  quit;
	                  %return;
				%end;
				%if (%upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC) and (&un < 7) %then %do;
					%put ERROR: At least 7 timepoints are recommended to fit cubic trends: You may try quadratic or linear trends instead;
				  	proc iml;
				   	  print "  Error Note:";
				   	  print "  At least 7 timepoints are recommended to fit cubic trends. You may try quadratic or linear trends instead. ";
				  	quit;
		  			%return;
				%end;
			%end;
		%end;

		%let uvar =; %let toprint2=;

		/* Compute the estimation for all different slopes estimation the user desires*/
		%let uj=0;
		%let uj= %eval(&uj + %_counts_(&slopes) );

		%if &uj=1 %then %do;

			%let uvar = &slopes;

			/* eMKF: Strings used to account for various combinations of available parameter estimates (up to cubic only) */
			%let emkfkeep = ; %let emkfrename = ;
			%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do; /* one outcome */
				%let emkfkeep = a; %let emkfrename = a=a_&uvar;
				%if %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC %then %do;
					%let emkfkeep = &emkfkeep b1 b2 b3;
					%let emkfrename = &emkfrename b1=b1_&uvar b2=b2_&uvar b3=b3_&uvar;
				%end;
				%if %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD %then %do;
					%let emkfkeep =  &emkfkeep b1 b2;
					%let emkfrename = &emkfrename b1=b1_&uvar b2=b2_&uvar;
				%end;
				%if %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR %then %do;
					%let emkfkeep = &emkfkeep b1;
					%let emkfrename = &emkfrename b1=b1_&uvar;
				%end;
			%end;
			%else %do; /* two outcomes */
				%let emkfkeep = o1a o2a; %let emkfrename = o1a=o1a_&uvar o2a=o2a_&uvar;
				%if %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC %then %do;
					%let emkfkeep = &emkfkeep o1b1 o2b1 o1b2 o2b2 o1b3 o2b3;
					%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar o1b2=o1b2_&uvar o2b2=o2b2_&uvar o1b3=o1b3_&uvar o2b3=o2b3_&uvar;
				%end;
				%if %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD %then %do;
					%let emkfkeep =  &emkfkeep o1b1 o2b1 o1b2 o2b2;
					%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar o1b2=o1b2_&uvar o2b2=o2b2_&uvar;
				%end;
				%if %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR %then %do;
					%let emkfkeep = &emkfkeep o1b1 o2b1;
					%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar;
				%end;
			%end;

			%let xtrakeep22=;
			%let xtrakeep22= &xtrakeep &outcome &se &group &time impute inputorder &by;
			%if &by ^=%str() %then %let xtrakeep22= &xtrakeep22 imputeb;
			%let toprint2=&slopes;

			/* eMKF: Simplified suffix assignment */
			%let newuvar = ;

			%if %upcase(&uvar) = INDEP_CUBIC   %then %let newuvar = _GC;
			%if %upcase(&uvar) = INDEP_QUAD    %then %let newuvar = _GQ;
			%if %upcase(&uvar) = INDEP_LINEAR  %then %let newuvar = _GL;
			%if %upcase(&uvar) = COMMON_CUBIC  %then %let newuvar = _1C;
			%if %upcase(&uvar) = COMMON_QUAD   %then %let newuvar = _1Q;
			%if %upcase(&uvar) = COMMON_LINEAR %then %let newuvar = _1L;
			%if %upcase(&uvar) = DROPPED	   %then %let newuvar = _0;

			%let _thekeep1 = &_thekeep1 pred_&uvar =pred&newuvar predse_&uvar=se&newuvar;
			%let _thekeep2 = &_thekeep2 pred_&uvar predse_&uvar;
			%let _thekeep3 = &_thekeep3 &out.(keep= pred_&uvar) &out.(keep= predse_&uvar);
			%let _thekeeps = &_thekeeps  pred&newuvar =&_oo1_._pred&newuvar  se&newuvar =&_oo1_._se&newuvar;

			/* eMKF: corrected issue in original MKF with re-labeling pred and se using outcome label when there are two outcomes */
			%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
				%let _thekeep1 = &_thekeep1 pred2_&uvar =pred2&newuvar pred2se_&uvar=se2&newuvar;
				%let _thekeep2 = &_thekeep2 pred2_&uvar pred2se_&uvar;
				%let _thekeep3 = &_thekeep3 &out.(keep= pred2_&uvar) &out.(keep= pred2se_&uvar);
				%let _thekeeps = &_thekeeps  pred2&newuvar =&_oo2_._pred&newuvar se2&newuvar= &_oo2_._se&newuvar;
			%end;

			%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;

				%htrp( data=_nlmixdata_, outcome=_y, se=_se, group=_groupnum, time=_time, by=&_ssby, 
					   xtrakeep=&xtrakeep22, orpoly=&orpoly,
					  _rho_=&_rho_ , _tausq_=&_tausq_ , bvalue= &uvar , DF=&DF, out=&out , print=&modelprint );

				data &out._pred;
				  set &out._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar )
				                 drop=_2loglike _rho_ _tausq_ 
									  &emkfkeep /* eMKF: Additional quadratic and cubic terms to drop */
 									  );; 
				run;

			%end;
			%else %do;

				%htrp2d( data=_nlmixdata_, outcome=_y, se=_se, outcome2=_y2, se2=_se2, group=_groupnum, time=_time, by=&_ssby,
						 xtrakeep=&xtrakeep22, orpoly=&orpoly,
						 _rho_=&_rho_ , _tausq_=&_tausq_ , bvalue= &uvar , DF=&DF, out=&out , print=&modelprint );

				data &out._pred;
				  set &out._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar 
				                         prediction2=pred2_&uvar predMSE2=pred2MSE_&uvar predSE2=pred2SE_&uvar )
				                 drop=_2loglike _rho_ _tausq_ 
									  &emkfkeep  /* eMKF: Additional quadratic and cubic terms to drop */
									  moddelta _tausq2_ _rho2_ err3 _se3 gamma gamma2 gammaeff predeff predeff2 
									  gamma_blup gamma_se pred_blup pred_blup2 
				                );; 
				run;

			%end;

			/* eMKF: clean-up. Deleting model-specific results but those could have useful information for the advanced user  */
			proc datasets nolist;
			  delete &out._covY &out._data &out._H &out._predvar ;
			run ;
			quit;

		%end; /*End of condition if uj = 1 */

		/* If there are more than one type of estimation methods requested then conduct a model averaging estimation */
		/* eMKF: 1=indep_cubic, 2=indep_quad, 3=indep_linear, 4=common_cubic, 5=common_quad, 6=common_linear, 7=dropped */

		%let flag1 =; %let flag2 =; %let flag3 =; %let  flag4=; %let flag5=; %let flag6=; %let flag7=;

 		%if &uj > 1 %then %do;

			/** eMKF: Check that a common reference/descendant model is included, and add to list if not:
			 ** - If both indep_quad   & common_cubic are selected, common_quad   will be added to the list if not already specified.
		  	 ** - If both indep_linear & common_cubic are selected, common_linear will be added to the list if not already specified.
		  	 ** - If both indep_linear & common_quad  are selected, common_linear will be added to the list if not already specified.
			 **/

			%let flag1=0; %let flag2=0; %let flag3=0; %let flag4=0; %let flag5=0; %let flag6=0; %let flag7=0;
			%let uvar=; %let ui = 0;
			%do ui=1 %to %_counts_(&slopes); /* eMKF: first pass to set flags */
				%let uvar=;
				%let uvar=%scan(&slopes, &ui);
				%if %upcase(&uvar) = INDEP_CUBIC   %then %let flag1 = 1;
				%if %upcase(&uvar) = INDEP_QUAD    %then %let flag2 = 1;
				%if %upcase(&uvar) = INDEP_LINEAR  %then %let flag3 = 1; 
				%if %upcase(&uvar) = COMMON_CUBIC  %then %let flag4 = 1;
				%if %upcase(&uvar) = COMMON_QUAD   %then %let flag5 = 1;
				%if %upcase(&uvar) = COMMON_LINEAR %then %let flag6 = 1; 
				%if %upcase(&uvar) = DROPPED 	   %then %let flag7 = 1; 
            %end;
			%let uvar=;

			%let _slopes = &slopes;

			%if &flag5 = 1 and &flag3 = 1 and &flag6 ^= 1 and &flag7 ^= 1 %then %do;
                %put ;	
				%put Both the common_quad and indep_linear models were specified with no shared descendant. Adding the common_linear model...;
			  	%let _slopes = &_slopes common_linear;
				%let flag6 = 1;
				%put &slopes;
				%put &_slopes;
			%end;
            %if &flag4 = 1 and &flag3 = 1 and &flag6 ^= 1 and &flag7 ^= 1 %then %do;
			    %put ;	
				%put Both the common_cubic and indep_linear models were specified with no shared descendant. Adding the common_linear model...;
			  	%let _slopes = &_slopes common_linear;
				%let flag6 = 1;
				%put &slopes;
				%put &_slopes;
			%end;
            %if &flag4 = 1 and &flag2 = 1 and &flag5 ^= 1 and &flag6 ^= 1 and &flag7 ^= 1 %then %do;				
			    %put ;	
				%put Both the common_cubic and indep_quad models were specified with no shared descendant. Adding the common_quad model...;
			  	%let _slopes = &_slopes common_quad;
				%let flag5 = 1;
				%put &slopes;
				%put &_slopes;
			%end;

			/*eMFK: End check for nested models */

			%let uvar=ModelAvg;
			%let _thekeep1 = &_thekeep1 pred_&uvar =pred_MA predse_&uvar=se_MA;
			%let _thekeep2 = &_thekeep2 pred_&uvar predse_&uvar;
			%let _thekeep3 = &_thekeep3 &out.(keep= pred_&uvar) &out.(keep= predse_&uvar);
			%let _thekeeps = &_thekeeps pred_MA=&_oo1_._pred_MA se_MA=&_oo1_._se_MA;

			%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
				%let _thekeep1 = &_thekeep1 pred2_&uvar =pred2_MA pred2se_&uvar=se2_MA;
				%let _thekeep2 = &_thekeep2 pred2_&uvar pred2se_&uvar;
				%let _thekeep3 = &_thekeep3 &out.(keep= pred2_&uvar) &out.(keep= pred2se_&uvar);
				%let _thekeeps = &_thekeeps pred2_MA=&_oo2_._pred_MA se2_MA=&_oo2_._se_MA;
			%end;
			
			data &out._pred &out._ests _loglike_;
			run;

			%let uvar=; %let ui=0;
			%do ui=1 %to %_counts_(&_slopes); 	/*eMKF: slopes replaced with _slopes */

				%let uvar=;
				%let uvar=%scan(&_slopes, &ui);

				/* eMKF: Strings used to account for various combinations of available parameter estimates (up to cubic only) */
				%let emkfkeep = ; %let emkfrename = ;
				%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do; /* one outcome */
					%let emkfkeep = a; %let emkfrename = a=a_&uvar;
					%if %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC %then %do;
						%let emkfkeep = &emkfkeep b1 b2 b3;
						%let emkfrename = &emkfrename b1=b1_&uvar b2=b2_&uvar b3=b3_&uvar;
					%end;
					%if %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD %then %do;
						%let emkfkeep =  &emkfkeep b1 b2;
						%let emkfrename = &emkfrename b1=b1_&uvar b2=b2_&uvar;
					%end;
					%if %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR %then %do;
						%let emkfkeep = &emkfkeep b1;
						%let emkfrename = &emkfrename b1=b1_&uvar;
					%end;
				%end;
				%else %do; /* two outcomes */
					%let emkfkeep = o1a o2a; %let emkfrename = o1a=o1a_&uvar o2a=o2a_&uvar;
					%if %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC %then %do;
						%let emkfkeep = &emkfkeep o1b1 o2b1 o1b2 o2b2 o1b3 o2b3;
						%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar o1b2=o1b2_&uvar o2b2=o2b2_&uvar o1b3=o1b3_&uvar o2b3=o2b3_&uvar;
					%end;
					%if %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD %then %do;
						%let emkfkeep =  &emkfkeep o1b1 o2b1 o1b2 o2b2;
						%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar o1b2=o1b2_&uvar o2b2=o2b2_&uvar;
					%end;
					%if %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR %then %do;
						%let emkfkeep = &emkfkeep o1b1 o2b1;
						%let emkfrename = &emkfrename o1b1=o1b1_&uvar o2b1=o2b1_&uvar;
					%end;
				%end;

				/* eMKF: Simplified suffix assignment */
				%let newuvar = ;

				%if %upcase(&uvar) = INDEP_CUBIC   %then %let newuvar = _GC;
				%if %upcase(&uvar) = INDEP_QUAD    %then %let newuvar = _GQ;
				%if %upcase(&uvar) = INDEP_LINEAR  %then %let newuvar = _GL;
				%if %upcase(&uvar) = COMMON_CUBIC  %then %let newuvar = _1C;
				%if %upcase(&uvar) = COMMON_QUAD   %then %let newuvar = _1Q;
				%if %upcase(&uvar) = COMMON_LINEAR %then %let newuvar = _1L;
				%if %upcase(&uvar) = DROPPED	   %then %let newuvar = _0;

				%let _thekeep1 = &_thekeep1 pred_&uvar =pred&newuvar  predse_&uvar=rmse&newuvar;
				%let _thekeep2 = &_thekeep2 pred_&uvar predse_&uvar;
				%let _thekeep3 = &_thekeep3 &out.(keep= pred_&uvar) &out.(keep= predse_&uvar);
				%let _thekeeps = &_thekeeps  pred&newuvar=&_oo1_._pred&newuvar  rmse&newuvar=&_oo1_._se&newuvar;

				%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
					%let _thekeep1 = &_thekeep1 pred2_&uvar =pred2&newuvar pred2se_&uvar=rmse2&newuvar;
					%let _thekeep2 = &_thekeep2 pred2_&uvar pred2se_&uvar;
					%let _thekeep3 = &_thekeep3 &out.(keep= pred2_&uvar) &out.(keep= pred2se_&uvar);
					%let _thekeeps = &_thekeeps  pred2&newuvar=&_oo2_._pred&newuvar  rmse2&newuvar=&_oo2_._se&newuvar;
				%end;

				%if %upcase(&uvar)=INDEP_CUBIC   %then %let flag1=YES;
				%if %upcase(&uvar)=INDEP_QUAD    %then %let flag2=YES;
				%if %upcase(&uvar)=INDEP_LINEAR  %then %let flag3=YES;
				%if %upcase(&uvar)=COMMON_CUBIC  %then %let flag4=YES;
				%if %upcase(&uvar)=COMMON_QUAD   %then %let flag5=YES;
				%if %upcase(&uvar)=COMMON_LINEAR %then %let flag6=YES;
				%if %upcase(&uvar)=DROPPED       %then %let flag7=YES;

				%let xtrakeep22=;
				%let xtrakeep22= &xtrakeep &outcome &se &group &time impute inputorder &by;
				%if &by ^=%str() %then %let xtrakeep22= &xtrakeep22 imputeb;

				%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %let xtrakeep22= &xtrakeep22 &outcome2 &se2;;

				%let toprint2=ModelAvg;

				%put Start &uvar;

				%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;
				  	%htrp( data=_nlmixdata_, outcome=_y, se=_se, group=_groupnum, time=_time, by=&_ssby, 
						   xtrakeep=&xtrakeep22, orpoly=&orpoly,
						   _rho_=&_rho_ , _tausq_=&_tausq_ , bvalue= &uvar , DF=&DF, out=_KF_&uvar , print=&modelprint );
				%end;
				%else %do;
					%htrp2d( data=_nlmixdata_, outcome=_y, se=_se, outcome2=_y2, se2=_se2, group=_groupnum, time=_time, by=&_ssby,  
							 xtrakeep=&xtrakeep22, orpoly=&orpoly,
							 _rho_=&_rho_ , _tausq_=&_tausq_ , bvalue= &uvar , DF=&DF, out=_KF_&uvar , print=&modelprint );

					data _KF_&uvar._pred;
				 	  set _KF_&uvar._pred(drop=moddelta err3 _se3 gamma gamma2 gammaeff predeff predeff2 gamma_blup gamma_se pred_blup pred_blup2);
					run;
			   	%end;

		       	%put End &uvar;

				/*eMKF: From this point on, _groupnum and _reps have been dropped: will work with generic _group_ and _rep directly */

			   	proc sort data=_KF_&uvar._pred;
				  by _rep _group_ _time;
				run;

				data _junk_;
				run;

				proc sort data=_KF_&uvar._pred out=_junk_ nodupkey;
				  by _rep _group_ ;
				run;

				proc sort data=_KF_&uvar._ests;
				  by _rep _group_ ;
				run;

				data _KF_&uvar._ests;
				  merge _KF_&uvar._ests _junk_(keep= _rep _group_);
				  by _rep _group_;
				run;

		 		proc sort data=_KF_&uvar._ests;
		  		  by _rep _group_;
		 		run;

		    	%if &ui=1 %then %do;
					%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;

			  			data &out._pred;
			    		  set _KF_&uvar._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar )
		                  		              drop=_2loglike _rho_ _tausq_ 
													&emkfkeep /* eMKF: additional quadratic and cubic terms to drop */
											  );; 
			   			run;

			  			data &out._ests;
			    		  set _KF_&uvar._ests(rename=(_2loglike=_2ll&uvar _rho_=_rho_&uvar._ _tausq_=_tausq_&uvar._  
													 &emkfrename  /* eMKF: additional quadratic and cubic terms to rename */
											 ));;  
			  			run;

					%end;
					%else %do;

			  			data &out._pred;
			    		  set _KF_&uvar._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar 
		                                        	 prediction2=pred2_&uvar predMSE2=pred2MSE_&uvar predSE2=pred2SE_&uvar)
		                            		  drop=_2loglike _rho_ _rho2_ _tausq_ _tausq2_ 
												   &emkfkeep  /* eMKF: additional quadratic and cubic terms to drop */
										 	 );; 
			  			run;

			  			data &out._ests;
			    		  set _KF_&uvar._ests(rename=(_2loglike=_2ll&uvar _rho_=_rho_&uvar._ _tausq_=_tausq_&uvar._  
													 &emkfrename /* eMKF: additional quadratic and cubic terms to rename */
											 ));; 
			  			run;

					%end;

					/* eMKF: Keep the log likelihood estimate for the type of model into llike1, 2, 3, 4, 5, 6, and 7*/

					data  _loglike_; /*eMKF: modified to allow for quad and cubic trend models */
					   set _KF_&uvar._ests(keep= _rep _2loglike );
					   llike1=.;
					   llike2=.;
					   llike3=.;
					   llike4=.;
					   llike5=.;
					   llike6=.;
					   llike7=.;
					   %if %upcase(&uvar)=INDEP_CUBIC  	%then llike1 = _2loglike ;;
					   %if %upcase(&uvar)=INDEP_QUAD 	%then llike2 = _2loglike ;; 
					   %if %upcase(&uvar)=INDEP_LINEAR 	%then llike3 = _2loglike ;;
					   %if %upcase(&uvar)=COMMON_CUBIC 	%then llike4 = _2loglike ;;
					   %if %upcase(&uvar)=COMMON_QUAD  	%then llike5 = _2loglike ;;
					   %if %upcase(&uvar)=COMMON_LINEAR %then llike6 = _2loglike ;;
					   %if %upcase(&uvar)=DROPPED      	%then llike7 = _2loglike ;;
					   if _2loglike ne .;
					run;

					proc sort data=_loglike_ nodupkey;
					   by _rep _2loglike;
					run;

					data _loglike_;
					   set _loglike_;
					   drop _2loglike;
					run;

				%end;

				%if &ui > 1 %then %do;

					data  _loglike2_;
					run;

					%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;

			  			data &out._pred;
			    			merge &out._pred _KF_&uvar._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar )
		              							              drop=_2loglike _rho_ _tausq_ 
																   &emkfkeep    /* eMKF: additional quadratic and cubic terms to drop */
															 );; 
							by _rep _group_ _time;
			  			run;

			  			data &out._ests;
			    			merge &out._ests _KF_&uvar._ests(rename=(_2loglike=_2ll&uvar _rho_=_rho_&uvar._ _tausq_=_tausq_&uvar._  
																	 &emkfrename /* eMKF: additional quadratic and cubic terms to rename */
															));; 
							by _rep _group_;
			  			run;

					%end;
					%else %do;

			  			data &out._pred;
			    			merge &out._pred _KF_&uvar._pred(rename=(prediction=pred_&uvar predMSE=predMSE_&uvar predSE=predSE_&uvar 
		         					                                 prediction2=pred2_&uvar predMSE2=pred2MSE_&uvar predSE2=pred2SE_&uvar )
		              							              drop= _2loglike _rho_ _rho2_ _tausq_ _tausq2_ 
																	&emkfkeep    /* eMKF: additional quadratic and cubic terms to drop */ 
															 );;
							by _rep _group_ _time;
			  			run;

			  			data &out._ests;
			    			merge &out._ests _KF_&uvar._ests(rename=(_2loglike=_2ll&uvar _rho_=_rho_&uvar._ _tausq_=_tausq_&uvar._  
																	 &emkfrename /* eMKF: additional quadratic and cubic terms to rename */
               												));; 
							by _rep _group_;
			  			run;

					%end;

					/* eMKF: Keep the log likelihood estimate for the type of model into llike1, 2, 3, 4, 5, 6, and 7*/

			  		data  _loglike2_; /* eMKF: modified to allow for quad and cubic trend models */
			   			set _KF_&uvar._ests(keep= _rep _2loglike);
			   			%if %upcase(&uvar)=INDEP_CUBIC   %then llike1 = _2loglike ;; 
			   			%if %upcase(&uvar)=INDEP_QUAD    %then llike2 = _2loglike ;;
						%if %upcase(&uvar)=INDEP_LINEAR  %then llike3 = _2loglike ;;
						%if %upcase(&uvar)=COMMON_CUBIC  %then llike4 = _2loglike ;;
						%if %upcase(&uvar)=COMMON_QUAD   %then llike5 = _2loglike ;;
						%if %upcase(&uvar)=COMMON_LINEAR %then llike6 = _2loglike ;;
						%if %upcase(&uvar)=DROPPED       %then llike7 = _2loglike ;;
			   			if _2loglike ne .;
			  		run;

			  		proc sort data=_loglike2_ nodupkey;
			   			by _rep _2loglike;
			  		run;

			  		data _loglike2_;
			   			set _loglike2_;
			   			drop _2loglike;
			  		run;

		      		data _loglike_;
			   			merge _loglike_ _loglike2_;
			   			by _rep;
			  		run;

				%end;

 			%end; /*End of do ui loop */

		    data _loglike_; /*eMKF: modified to allow for quad and cubic trend models */
				 set _loglike_;

				 df1 = 3 * &ug; 
				 df2 = 2 * &ug;
				 df3 = 1 * &ug; 
				 df4 = 3; 
				 df5 = 2;
				 df6 = 1; 
				 df7 = 0; 

				 ntime=&un;

				 difll67=llike7 - llike6;
				 difll57=llike7 - llike5;
				 difll47=llike7 - llike4;
				 difll37=llike7 - llike3;
				 difll27=llike7 - llike2;
				 difll17=llike7 - llike1;
				 difll56=llike6 - llike5;
				 difll46=llike6 - llike4;
				 difll36=llike6 - llike3;
				 difll26=llike6 - llike2;
				 difll16=llike6 - llike1;
				 difll45=llike5 - llike4;						 
				 difll35=llike5 - llike3;
				 difll25=llike5 - llike2;
				 difll15=llike5 - llike1;
				 difll34=llike4 - llike3;
				 difll24=llike4 - llike2;
				 difll14=llike4 - llike1;
				 difll23=llike3 - llike2;
				 difll13=llike3 - llike1;
				 difll12=llike2 - llike1;
				 
				 if difll67 < 0 and difll67 ne . then difll67=0;
				 if difll57 < 0 and difll57 ne . then difll57=0;
				 if difll47 < 0 and difll47 ne . then difll47=0;
				 if difll37 < 0 and difll37 ne . then difll37=0;
				 if difll27 < 0 and difll27 ne . then difll27=0;
				 if difll17 < 0 and difll17 ne . then difll17=0;
				 if difll56 < 0 and difll56 ne . then difll56=0;
				 if difll46 < 0 and difll46 ne . then difll46=0;
				 if difll36 < 0 and difll36 ne . then difll36=0;
				 if difll26 < 0 and difll26 ne . then difll26=0;
				 if difll16 < 0 and difll16 ne . then difll16=0;
				 if difll45 < 0 and difll45 ne . then difll45=0;
				 if difll35 < 0 and difll35 ne . then difll35=0; 
				 if difll25 < 0 and difll25 ne . then difll25=0;
				 if difll15 < 0 and difll15 ne . then difll15=0;
				 if difll34 < 0 and difll34 ne . then difll34=0; 
				 if difll24 < 0 and difll24 ne . then difll24=0;
				 if difll14 < 0 and difll14 ne . then difll14=0;
				 if difll23 < 0 and difll23 ne . then difll23=0;
				 if difll13 < 0 and difll13 ne . then difll13=0;
				 if difll12 < 0 and difll12 ne . then difll12=0;

				 df67=df6-df7;
				 df57=df5-df7;
				 df47=df4-df7;
				 df37=df3-df7;
				 df27=df2-df7;
				 df17=df1-df7;
				 df56=df5-df6;
				 df46=df4-df6;
				 df36=df3-df6;
				 df26=df2-df6;
				 df16=df1-df6;
				 df45=df4-df5;
				 df35=df3-df5;
				 df25=df2-df5;
				 df15=df1-df5;
				 df34=df3-df4;
				 df24=df2-df4;
				 df14=df1-df4;
				 df23=df2-df3;
				 df13=df1-df3;
				 df12=df1-df2;

				 bic67= df67*log(&ug * &un) - difll67;
				 bic57= df57*log(&ug * &un) - difll57; 
				 bic47= df47*log(&ug * &un) - difll47;
				 bic37= df37*log(&ug * &un) - difll37;
				 bic27= df27*log(&ug * &un) - difll27;
				 bic17= df17*log(&ug * &un) - difll17;
				 bic56= df56*log(&ug * &un) - difll56;
				 bic46= df46*log(&ug * &un) - difll46;
				 bic36= df36*log(&ug * &un) - difll36;
				 bic26= df26*log(&ug * &un) - difll26;
				 bic16= df16*log(&ug * &un) - difll16;
				 bic45= df45*log(&ug * &un) - difll45;
				 bic35= df35*log(&ug * &un) - difll35;
				 bic25= df25*log(&ug * &un) - difll25;
				 bic15= df15*log(&ug * &un) - difll15;
				 bic34= df34*log(&ug * &un) - difll34;
				 bic24= df24*log(&ug * &un) - difll24;
				 bic14= df14*log(&ug * &un) - difll14;
				 bic23= df23*log(&ug * &un) - difll23;
				 bic13= df13*log(&ug * &un) - difll13;
				 bic12= df12*log(&ug * &un) - difll12;

				 /* eMKF: Use coalesce to force the Bayes factor to 0 when missing */
				 bf67 = coalesce(exp(-0.5*bic67), 0);
				 bf57 = coalesce(exp(-0.5*bic57), 0);
				 bf47 = coalesce(exp(-0.5*bic47), 0);
				 bf37 = coalesce(exp(-0.5*bic37), 0);
				 bf27 = coalesce(exp(-0.5*bic27), 0);
				 bf17 = coalesce(exp(-0.5*bic17), 0);
				 bf56 = coalesce(exp(-0.5*bic56), 0);
				 bf46 = coalesce(exp(-0.5*bic46), 0);
				 bf36 = coalesce(exp(-0.5*bic36), 0);
				 bf26 = coalesce(exp(-0.5*bic26), 0);
				 bf16 = coalesce(exp(-0.5*bic16), 0);
				 bf45 = coalesce(exp(-0.5*bic45), 0);
				 bf35 = coalesce(exp(-0.5*bic35), 0); 
				 bf25 = coalesce(exp(-0.5*bic25), 0);
				 bf15 = coalesce(exp(-0.5*bic15), 0);
				 bf34 = coalesce(exp(-0.5*bic34), 0);
				 bf24 = coalesce(exp(-0.5*bic24), 0);
				 bf14 = coalesce(exp(-0.5*bic14), 0);
				 bf23 = coalesce(exp(-0.5*bic23), 0);
				 bf13 = coalesce(exp(-0.5*bic13), 0);
				 bf12 = coalesce(exp(-0.5*bic12), 0);

				 /* Let's setup something in case one of the models did not converge */
				 /* eMKF: Expanded to account for various combinations up to 7 models */
				 /* eMKF: Instead of working through all the combinatorial sequences, BF was set to 0 when missing */

				 /* eMKF: Identify simplest available model and use it as reference for the Bayes factor-based weight calculations. */ 

				 /* eMKF: Model 7 was requested, is the simplest, and at least one other model was requested and converged */
				 if difll67 ne . or difll57 ne . or difll47 ne . or difll37 ne . or difll27 ne . or difll17 ne . then do; 
				 	p1 = bf17 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 	p2 = bf27 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 	p3 = bf37 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 	p4 = bf47 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 	p5 = bf57 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 	p6 = bf67 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
					p7 =    1 /(1+bf17+bf27+bf37+bf47+bf57+bf67);
				 end;
				 else do;
				 	/* eMKF: Model 6 was requested, is the simplest, and at least one other model was requested and converged */
					if difll56 ne . or difll46 ne . or difll36 ne . or difll26 ne . or difll16 ne . then do; 
						p1 = bf16 /(1+bf16+bf26+bf36+bf46+bf56);
					 	p2 = bf26 /(1+bf16+bf26+bf36+bf46+bf56);
					 	p3 = bf36 /(1+bf16+bf26+bf36+bf46+bf56);
					 	p4 = bf46 /(1+bf16+bf26+bf36+bf46+bf56);
					 	p5 = bf56 /(1+bf16+bf26+bf36+bf46+bf56);
					 	p6 =    1 /(1+bf16+bf26+bf36+bf46+bf56);
						p7 =    0;
					end;
					else do;
						/* eMKF: Model 5 was requested, is the simplest, and at least one other model was requested and converged.
					     * This will be superseded by the previous clauses when both models 3 and 5 are selected 
						 * because the code ensure models 6 or 7 are included 
						 */
						if difll45 ne . or difll35 ne . or difll25 ne . or difll15 ne . then do; 
							p1 = bf15 /(1+bf15+bf25+bf35+bf45);
							p2 = bf25 /(1+bf15+bf25+bf35+bf45);
						 	p3 = bf35 /(1+bf15+bf25+bf35+bf45);
						 	p4 = bf45 /(1+bf15+bf25+bf35+bf45);
						 	p5 =    1 /(1+bf15+bf25+bf35+bf45);
						 	p6 =    0;
							p7 =    0;
						end;
				 		else do;
							/* eMKF: Model 4 was requested, is the simplest, and at least one other model was requested and converged.
						 	 * This will be superseded by the previous clauses when both models 2 and 4 are selected 
							 *  because the code will ensure models 5, 6 or 7 is included.
						 	 * Similarly if both models 3 and 4 had been selected 
						 	 */
							if difll34 ne . or difll24 ne . or difll14 ne . then do; 
							 	p1 = bf14 /(1+bf14+bf24+bf34);
							 	p2 = bf24 /(1+bf14+bf24+bf34);
							 	p3 = bf34 /(1+bf14+bf24+bf34);
							 	p4 = 	1 /(1+bf14+bf24+bf34);
							 	p5 =    0;
							 	p6 =    0;
								p7 =    0;
							end;
							else do;
						 		/* eMKF: Model 3 was requested, is the simplest, and at least one other model was requested and converged */
							 	if difll23 ne . or difll13 ne . then do; 
							 		p1 = bf13 /(1+bf13+bf23);
							 		p2 = bf23 /(1+bf13+bf23);
							 		p3 = 	1 /(1+bf13+bf23);
							 		p4 = 	0;
							 		p5 =    0;
							 		p6 =    0;
									p7 =    0;
							 	end;
						 	 	else do;
						 	 		/* eMKF: Models 2 and 1 were requested and converged */
								 	if difll12 ne . then do; 
								 		p1 = bf12 /(1+bf12);
								 		p2 =    1 /(1+bf12);
								 		p3 = 	0;
								 		p4 = 	0;
								 		p5 =    0;
								 		p6 =    0;
										p7 =    0;
									end;
				 				end;
							end;
						end;
					end;
				end;
		    run; /*eMKF: runs data step */

			data _KF_covY;
			run;

			data _KF_covY;
			  merge %if %upcase(&flag1)=YES %then _KF_INDEP_CUBIC_covY(keep=_rep _group_ _time)  ;
                    %if %upcase(&flag2)=YES %then _KF_INDEP_QUAD_covY(keep=_rep _group_ _time)  ;
                    %if %upcase(&flag3)=YES %then _KF_INDEP_LINEAR_covY(keep=_rep _group_ _time)  ;
			        %if %upcase(&flag4)=YES %then _KF_COMMON_CUBIC_covY(keep=_rep _group_ _time); 
                    %if %upcase(&flag5)=YES %then _KF_COMMON_QUAD_covY(keep=_rep _group_ _time); 
                    %if %upcase(&flag6)=YES %then _KF_COMMON_LINEAR_covY(keep=_rep _group_ _time);  
			        %if %upcase(&flag7)=YES %then _KF_DROPPED_covY(keep=_rep _group_ _time); ;;
			  by _rep _group_ _time ;
			  if _time=. then delete;
			run;

			proc sort data=_loglike_;
			  by _rep;
			run;

			proc sort data= &out._pred;
			  by _rep;
			run;

			/* eMKF: revised averaging to account for re-numbering of models and up to 7 models requested */
			/* eMKF: used coalesce to force to 0 when missing instead of dealing with combinatorial number of model combinations */

			data &out._pred;		
			  merge &out._pred _loglike_(keep= _rep p1 p2 p3 p4 p5 p6 p7);
			  by _rep;

              pred_INDEP_CUBIC  	= coalesce(pred_INDEP_CUBIC, 0); 
			  pred_INDEP_QUAD   	= coalesce(pred_INDEP_QUAD, 0);
			  pred_INDEP_LINEAR  	= coalesce(pred_INDEP_LINEAR, 0);
              pred_COMMON_CUBIC 	= coalesce(pred_COMMON_CUBIC, 0);
			  pred_COMMON_QUAD  	= coalesce(pred_COMMON_QUAD, 0);
			  pred_COMMON_LINEAR 	= coalesce(pred_COMMON_LINEAR, 0);
			  pred_DROPPED      	= coalesce(pred_DROPPED, 0);

              predMSE_INDEP_CUBIC  	= coalesce(predMSE_INDEP_CUBIC, 0);
			  predMSE_INDEP_QUAD   	= coalesce(predMSE_INDEP_QUAD, 0);
			  predMSE_INDEP_LINEAR  = coalesce(predMSE_INDEP_LINEAR, 0);
              predMSE_COMMON_CUBIC 	= coalesce(predMSE_COMMON_CUBIC, 0);
			  predMSE_COMMON_QUAD  	= coalesce(predMSE_COMMON_QUAD, 0);
			  predMSE_COMMON_LINEAR = coalesce(predMSE_COMMON_LINEAR, 0);
			  predMSE_DROPPED      	= coalesce(predMSE_DROPPED, 0);

              pred_ModelAvg    = p1*pred_INDEP_CUBIC  + p2*pred_INDEP_QUAD  + p3*pred_INDEP_LINEAR  + 
		 						 p4*pred_COMMON_CUBIC + p5*pred_COMMON_QUAD + p6*pred_COMMON_LINEAR + p7*pred_DROPPED;

              predMSE_ModelAvg = p1*predMSE_INDEP_CUBIC 	+ p2*predMSE_INDEP_QUAD  + p3*predMSE_INDEP_LINEAR 	+ 
								 p4*predMSE_COMMON_CUBIC + p5*predMSE_COMMON_QUAD + p6*predMSE_COMMON_LINEAR + p7*predMSE_DROPPED;

              predSE_ModelAvg= sqrt(predMSE_ModelAvg);

			run;

			%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
				data &out._pred;
				  set &out._pred;

				  pred2_INDEP_CUBIC  	 = coalesce(pred2_INDEP_CUBIC, 0);
			 	  pred2_INDEP_QUAD   	 = coalesce(pred2_INDEP_QUAD, 0);
			 	  pred2_INDEP_LINEAR  	 = coalesce(pred2_INDEP_LINEAR, 0);
                  pred2_COMMON_CUBIC 	 = coalesce(pred2_COMMON_CUBIC, 0);
			 	  pred2_COMMON_QUAD  	 = coalesce(pred2_COMMON_QUAD, 0);
			 	  pred2_COMMON_LINEAR    = coalesce(pred2_COMMON_LINEAR, 0);
			 	  pred2_DROPPED      	 = coalesce(pred2_DROPPED, 0);

                  pred2MSE_INDEP_CUBIC   = coalesce(pred2MSE_INDEP_CUBIC, 0);
			 	  pred2MSE_INDEP_QUAD    = coalesce(pred2MSE_INDEP_QUAD, 0);
			 	  pred2MSE_INDEP_LINEAR  = coalesce(pred2MSE_INDEP_LINEAR, 0);
                  pred2MSE_COMMON_CUBIC  = coalesce(pred2MSE_COMMON_CUBIC, 0);
			 	  pred2MSE_COMMON_QUAD   = coalesce(pred2MSE_COMMON_QUAD, 0);
			 	  pred2MSE_COMMON_LINEAR = coalesce(pred2MSE_COMMON_LINEAR, 0);
			 	  pred2MSE_DROPPED       = coalesce(pred2MSE_DROPPED, 0);

				  pred2_ModelAvg    = p1*pred2_INDEP_CUBIC  + p2*pred2_INDEP_QUAD  + p3*pred2_INDEP_LINEAR  + 
				   				      p4*pred2_COMMON_CUBIC + p5*pred2_COMMON_QUAD + p6*pred2_COMMON_LINEAR + p7*pred2_DROPPED;

                  pred2MSE_ModelAvg = p1*pred2MSE_INDEP_CUBIC  + p2*pred2MSE_INDEP_QUAD  + p3*pred2MSE_INDEP_LINEAR  + 
								      p4*pred2MSE_COMMON_CUBIC + p5*pred2MSE_COMMON_QUAD + p6*pred2MSE_COMMON_LINEAR + p7*pred2MSE_DROPPED;

				  pred2SE_ModelAvg= sqrt(pred2MSE_ModelAvg);

				run;
			%end;

			/* eMKF: clean-up */
			proc datasets nolist;
			  delete _junk_ _loglike_ _loglike2_ _KF_covY
			  		 /* eMKF: Deleting model-specific results but those could have useful information for the advanced user */
					 %if %upcase(&flag1)=YES %then _KF_INDEP_CUBIC_: ;
		             %if %upcase(&flag2)=YES %then _KF_INDEP_QUAD_: ;
	                 %if %upcase(&flag3)=YES %then _KF_INDEP_LINEAR_: ;
	                 %if %upcase(&flag4)=YES %then _KF_COMMON_CUBIC_: ;
					 %if %upcase(&flag5)=YES %then _KF_COMMON_QUAD_: ;
					 %if %upcase(&flag6)=YES %then _KF_COMMON_LINEAR_: ;
					 %if %upcase(&flag7)=YES %then _KF_DROPPED_: ;
			        ;
			run ;
			quit;

		%end; /*End of condition if uj>1 */

	%end;  /*End of slopes estimations using MLE */

    /*********************************/
	/* Estimate the Bayesian model(s)*/
	/*********************************/

	%if %upcase(&Bayesian) = YES %then %do;

		/* eMKF: Total number of Bayesian models requested */
		%let uii = %_counts_(&Bayesmodel);

		/* eMKF: Guard against too many groups/timepoints: code for current implementation cannot handle more than 204 groups or more than 5508 data points */
		%if (&ug > 204 or %eval(&ug*&un) > 5508) %then %do;
			%put ERROR: Code for current implementation of this macro cannot handle more than 204 groups or more than 5508 data points (group by time) per stratum;
			%put ERROR- Please reduce the number of groups and/or timepoints;
			proc iml;
			  print "  Error Note:";
			  print "  Code for current implementation of this macro cannot handle more than 204 groups or more than 5508 data points (group by time) per stratum. ";
			  print "  Please reduce the number of groups and/or timepoints. ";
			quit;
			%return;
		%end;

		/* eMKF: Make sure there are enough time points to fit the desired trends (up to cubic) */

		%if (&un < 2) %then %do; /*eMKF: minimal check for 2+ points */
			%put ERROR: There are not enough timepoints to work with! At least 2 timepoints are required (4+ recommended) to use this macro;
			proc iml;
			  print "  Error Note:";
			  print "  There are not enough timepoints to work with! At least 2 timepoints are required (4+ recommended) to use this macro. ";
			quit;
			%return;
		%end;

		%if %upcase(&checkSampleSize) = YES %then %do;

			%if (&un < 4) %then %do;
				%put ERROR: There are not enough timepoints to work with! At least 4 timepoints are recommended to use this macro;
				proc iml;
				   print "  Error Note:";
				   print "  There are not enough timepoints to work with! At least 4 timepoints are recommended to use this macro. ";
				quit;
				%return;
			%end;

			%let ui=0; 	        
			%do ui=1 %to &uii;
				%let uvar=;
				%let uvar=%scan(&Bayesmodel, &ui);

				%if (%upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR or 
					 %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = BMA_LINEAR) and (&un < 5) %then %do;
					  %put ERROR: There are not enough time points to fit linear trends: You may try dropping the time trends instead;
					  proc iml;
					   print "  Error Note:";
					   print "  At least 5 timepoints are recommended to fit linear trends. You may try dropping the time trends instead. ";
					  quit;
					  %return;
				%end;
				%if (%upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = COMMON_QUAD or 
					 %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = BMA_QUAD) and (&un < 6) %then %do;
					  %put ERROR: At least 6 timepoints are recommended to fit quadratic trends: You may try linear trends instead;
					  proc iml;
					   print "  Error Note:";
					   print "  At least 6 timepoints are recommended to fit quadratic trends. You may try linear trends instead. ";
					  quit;
	                  %return;
				%end;
				%if (%upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = COMMON_CUBIC or 
					 %upcase(&uvar) = FULL_CUBIC or %upcase(&uvar) = BMA_CUBIC) and (&un < 7) %then %do;
					%put ERROR: At least 7 timepoints are recommended to fit cubic trends: You may try quadratic or linear trends instead;
				  	proc iml;
				   	 print "  Error Note:";
				   	 print "  At least 7 timepoints are recommended to fit cubic trends. You may try quadratic or linear trends instead. ";
				  	quit;
		  			%return;
				%end;
			%end;

			%let uvar =;
		%end;

		/* eMKF: Loop through sequence of requested model(s) to set model flag(s) */
		%let flag1a =0; %let flag2a =0; %let flag3a =0;
		%let flag1f =0; %let flag2f =0; %let flag3f =0;
		%let flag1  =0; %let flag2  =0; %let flag3  =0; 
		%let flag4  =0; %let flag5  =0; %let flag6  =0; %let flag7=0;
		%let ui = 0; 
		%do ui=1 %to &uii;
			%let uvar=;
			%let uvar=%scan(&Bayesmodel, &ui);
			%if %upcase(&uvar) = BMA_CUBIC     %then %let flag1a = 1;
			%if %upcase(&uvar) = BMA_QUAD      %then %let flag2a = 1;
			%if %upcase(&uvar) = BMA_LINEAR    %then %let flag3a = 1; 
			%if %upcase(&uvar) = FULL_CUBIC    %then %let flag1f = 1;
			%if %upcase(&uvar) = FULL_QUAD     %then %let flag2f = 1;
			%if %upcase(&uvar) = FULL_LINEAR   %then %let flag3f = 1; 
			%if %upcase(&uvar) = INDEP_CUBIC   %then %let flag1 = 1;
			%if %upcase(&uvar) = INDEP_QUAD    %then %let flag2 = 1;
			%if %upcase(&uvar) = INDEP_LINEAR  %then %let flag3 = 1; 
			%if %upcase(&uvar) = COMMON_CUBIC  %then %let flag4 = 1;
			%if %upcase(&uvar) = COMMON_QUAD   %then %let flag5 = 1;
			%if %upcase(&uvar) = COMMON_LINEAR %then %let flag6 = 1; 
			%if %upcase(&uvar) = DROPPED 	   %then %let flag7 = 1; 
        %end;
		%let uvar=;

		/* eMKF: If there are multiple Bayesian models requested and BayesmodelAvg = YES, 
		         or if bma_*** models are requested, then conduct a model averaging estimation.
				 Also over-rides BMA if fully Bayesian models are listed */

		%let _BMAmodel = ;
		%if &BayesmodelAvg = %str() %then %let BayesmodelAvg = NO;
		%if %upcase(&BayesmodelAvg) = YES and (&flag1f = 1 or &flag2f = 1 or &flag3f = 1) %then %do;
		    %put ;
			%put Because a fully Bayesian model was requested, no model averaging will be applied;
			%let BayesmodelAvg = NO;
		%end;
		%if %upcase(&BayesmodelAvg) = YES and &flag1a ^= 1 and &flag2a ^= 1 and &flag3a ^= 1 and &uii = 1 %then %do;
		    %put ;
			%put Because only one Bayesian model was specified, no model averaging will be applied;
			%let BayesmodelAvg = NO;
		%end;
		%if &flag1a = 1 or ((&flag1 = 1 or &flag4 = 1) and (&BayesmodelAvg = YES)) %then %do;
		    %put ;	
			%if %upcase(&BayesmodelAvg) = YES %then %put Because a cubic trend model was requested, all possible models up to cubic will be included in Bayesian model averaging;;
			%let _BMAmodel = BMA_CUBIC;
			%let flag1a = 1;
			%let BayesmodelAvg = YES;
		%end;
		%if &flag1a ^= 1 and &flag1f ^= 1 and &flag1 ^= 1 and &flag4 ^= 1 and 
			(&flag2a = 1 or ((&flag2 = 1 or &flag5 = 1) and (&BayesmodelAvg = YES))) %then %do;
		    %put ;	
			%if %upcase(&BayesmodelAvg) = YES %then %put Because a quadratic trend model was requested, all possible models up to quadratic will be included in Bayesian model averaging;;
			%let _BMAmodel = BMA_QUAD;
			%let flag2a = 1;
			%let BayesmodelAvg = YES;
		%end;
		%if &flag1a ^= 1 and &flag1f ^= 1 and &flag1 ^= 1 and &flag4 ^= 1 and 
			&flag2a ^= 1 and &flag2f ^= 1 and &flag2 ^= 1 and &flag5 ^= 1 and 
			(&flag3a = 1 or ((&flag3 = 1 or &flag6 = 1) and (&BayesmodelAvg = YES))) %then %do;		
		    %put ;	
			%if %upcase(&BayesmodelAvg) = YES %then %put Because a linear trend model was requested, all possible models up to linear will be included in Bayesian model averaging;;
			%let _BMAmodel = BMA_LINEAR;
			%let flag3a = 1;
			%let BayesmodelAvg = YES;
		%end;
		%if &flag1a ^= 1 and &flag1f ^= 1 and &flag1 ^= 1 and &flag4 ^= 1 and 
			&flag2a ^= 1 and &flag2f ^= 1 and &flag2 ^= 1 and &flag5 ^= 1 and 
			&flag3a ^= 1 and &flag3f ^= 1 and &flag3 ^= 1 and &flag6 ^= 1 %then %do;		
		    %put ;	
			%if %upcase(&BayesmodelAvg) = YES %then %put Because no trend model was requested, an intercepts-only model will be fit;;
			%let _BMAmodel = DROPPED;
			%let flag7 = 1;
			%let BayesmodelAvg = NO;
		%end;

		%let toprint=;

		/*eMKF: Find the number of replications (from stratification variable _rep) */
		%let crep=0;
		data _freq_;
		run;
		proc freq data=_bayesdata_ noprint;
		  tables _rep /list out=_freq_;
		run;
		data _null_;
		  set _freq_;
		  _crep_ +1;
		  call symput('crep',_crep_);
		  keep _rep _crep_;
		run;
		%let crep= %eval(0 + &crep);

		proc datasets nolist;
  			delete _freq_;
		run ;

		/* eMKF: nmc and thin macro variables are needed below outside of call to PROC MCMC */
		%if &thin = %str() %then %let thin = 1; 	 /* set to proc mcmc default if missing */ 
		%if &nbi = %str() %then %let nbi = 1000; 	 /* set to proc mcmc default if missing */
		%if &nmc = %str() %then %let nmc = 1000; 	 /* set to proc mcmc default if missing */

		/***********************************************************************************/
		/* eMKF: Pre-compile appropriate Gibbs sampler subroutine(s) for selected model(s) */
 		/***********************************************************************************/

		/* eMKF: Set CMP library if not provided. Could also use sasuser.funcs if user has write-permission */
		%if &cmploc = %str() %then %let cmploc = work.funcs; 

		/* eMKF: true states */
		%gibbs_uds_compile_EP(g=&ug, n=&un, loc=&cmploc);

		/* eMKF: variance parameters */
		%if %upcase(&randomVars) = YES %then %gibbs_uds_compile_RP(g=&ug, n=&un, loc=&cmploc);; 

		/* eMKF: mean hyperparameters + mbetag/Dbetag arrays in the fully Bayesian models */
		%if &flag1f = 1 %then %gibbs_uds_compile_MP(uvar=full_cubic,  g=&ug, loc=&cmploc);;
		%if &flag2f = 1 %then %gibbs_uds_compile_MP(uvar=full_quad,   g=&ug, loc=&cmploc);;
		%if &flag3f = 1 %then %gibbs_uds_compile_MP(uvar=full_linear, g=&ug, loc=&cmploc);;

		/* eMKF: regression coefficients and predictions */
		%if &flag1a = 1 			  %then %gibbs_uds_compile_CP(uvar=bma_cubic,     g=&ug, n=&un, loc=&cmploc);; 
		%if &flag2a = 1 			  %then %gibbs_uds_compile_CP(uvar=bma_quad,      g=&ug, n=&un, loc=&cmploc);; 
		%if &flag3a = 1 			  %then %gibbs_uds_compile_CP(uvar=bma_linear,    g=&ug, n=&un, loc=&cmploc);; 
		%if &flag1 = 1 or &flag1f = 1 %then %gibbs_uds_compile_CP(uvar=indep_cubic,   g=&ug, n=&un, loc=&cmploc);; 
		%if &flag2 = 1 or &flag2f = 1 %then %gibbs_uds_compile_CP(uvar=indep_quad,    g=&ug, n=&un, loc=&cmploc);; 
		%if &flag3 = 1 or &flag3f = 1 %then %gibbs_uds_compile_CP(uvar=indep_linear,  g=&ug, n=&un, loc=&cmploc);; 
		%if &flag4 = 1 				  %then %gibbs_uds_compile_CP(uvar=common_cubic,  g=&ug, n=&un, loc=&cmploc);; 
		%if &flag5 = 1 				  %then %gibbs_uds_compile_CP(uvar=common_quad,   g=&ug, n=&un, loc=&cmploc);; 
		%if &flag6 = 1 				  %then %gibbs_uds_compile_CP(uvar=common_linear, g=&ug, n=&un, loc=&cmploc);; 
		%if &flag7 = 1 				  %then %gibbs_uds_compile_CP(uvar=dropped, 	  g=&ug, n=&un, loc=&cmploc);; 

		/* eMKF: model flags in Bayesian model averaging */
		%if &flag1a = 1 			  %then %gibbs_uds_compile_FP(uvar=bma_cubic,     g=&ug, n=&un, loc=&cmploc);; 
		%if &flag2a = 1 			  %then %gibbs_uds_compile_FP(uvar=bma_quad,      g=&ug, n=&un, loc=&cmploc);; 
		%if &flag3a = 1 			  %then %gibbs_uds_compile_FP(uvar=bma_linear,    g=&ug, n=&un, loc=&cmploc);; 

		/* eMKF: Start model fitting loop(s) (tuning loop from original MKF is now incorporated in the call to proc mcmc) */
		%put Start Bayesian model fitting loop(s); 

		/* eMKF: to avoid larger dataset size than necessary, initialize one posterior log file per replication */
		data _bayesdata1_ &out._bayes &out._bayesparm &out._bayeslogGR ;
		run;
		%let uj=0;
		%do uj=1 %to &crep; 
			data &out._bayeslog_rep&uj; 
			run;
		%end; 

		/* eMKF: additional log file to hold posterior model weights in Bayesian model averaging */
		%if &flag1a = 1 or &flag2a = 1 or &flag3a = 1 %then %do;
			data &out._bayeslogmod;
			run;
		%end;

		/* eMKF: Total number of models to loop through */
		%let uii=0;
		%if %upcase(&BayesmodelAvg) ^= YES %then %let uii = %_counts_(&Bayesmodel);;
		%if %upcase(&BayesmodelAvg) = YES  %then %let uii = %_counts_(&_BMAmodel);;

		/* Loop over replications */
		%let uj=0;
		%do uj=1 %to &crep;
			data _bayesdata1_;
			  set _bayesdata_;
			  if _rep= &uj;
			run;

			data _bayesfit2_ _bayesparam2_ _bayeslog2_ _bayeslogGR2_ ;
			run;

			%if &flag1a = 1 or &flag2a = 1 or &flag3a = 1 %then %do;
				data _bayeslogmod2_;
				run;
			%end;

			/* Loop over the sequence of specified trend models */
			%let ui=0; 
			%do ui=1 %to &uii;

				/* eMKF: Disabled option from original MKF that allowed user to include F, say, instead of FULL_LINEAR in Bayesmodel */
				%let uvar=;
				%if %upcase(&BayesmodelAvg) ^= YES %then %let uvar= %scan(&Bayesmodel, &ui);;
				%if %upcase(&BayesmodelAvg) = YES %then %let uvar= %scan(&_BMAmodel, &ui);;

				/* eMKF: Simplified toprint suffixes, assuming input is in lowercase (lowcase macro function could be used as an alternative) */
				%let toprint = &uvar;

				/* eMKF: Simplified suffix construction for _thekeep* variables */
				%let newuvar = ;
				%if %upcase(&uvar) = BMA_CUBIC     %then %let newuvar = _BMAC;
				%if %upcase(&uvar) = FULL_CUBIC    %then %let newuvar = _BFC;
				%if %upcase(&uvar) = INDEP_CUBIC   %then %let newuvar = _BGC;
				%if %upcase(&uvar) = COMMON_CUBIC  %then %let newuvar = _B1C;
				%if %upcase(&uvar) = BMA_QUAD      %then %let newuvar = _BMAQ;
				%if %upcase(&uvar) = FULL_QUAD     %then %let newuvar = _BFQ;
				%if %upcase(&uvar) = INDEP_QUAD    %then %let newuvar = _BGQ;
				%if %upcase(&uvar) = COMMON_QUAD   %then %let newuvar = _B1Q;
				%if %upcase(&uvar) = BMA_LINEAR    %then %let newuvar = _BMAL;
				%if %upcase(&uvar) = FULL_LINEAR   %then %let newuvar = _BFL;
				%if %upcase(&uvar) = INDEP_LINEAR  %then %let newuvar = _BGL;
				%if %upcase(&uvar) = COMMON_LINEAR %then %let newuvar = _B1L;
				%if %upcase(&uvar) = DROPPED	   %then %let newuvar =  _B0;

				%if &uj = 1 %then %do;
					%let _thekeep1 = &_thekeep1 pred_Bayes_&uvar =pred&newuvar predse_Bayes_&uvar=rmse&newuvar;
					%let _thekeep1b = &_thekeep1b pred_Bayes_&uvar =pred&newuvar predse_Bayes_&uvar=rmse&newuvar;
					%let _thekeeps = &_thekeeps  pred&newuvar=&_oo1_._pred&newuvar  rmse&newuvar=&_oo1_._se&newuvar;
					%let _thekeepsb = &_thekeepsb  pred&newuvar=&_oo1_._pred&newuvar  rmse&newuvar=&_oo1_._se&newuvar;
					%let _thekeep2 = &_thekeep2 pred_Bayes_&uvar predse_Bayes_&uvar;
					%let _thekeep3 = &_thekeep3 &out.(keep= pred_Bayes_&uvar) &out.(keep= predse_Bayes_&uvar);
				%end;

				/* Run the Bayesian model for each chain */
				data _bayesfit_ _bayesparam_ _bayeslog_ ;
				run;

				%let uk=0;
				%do uk=1 %to &chains;

					%let _chainseed=;
					%let _chainseed=%eval(&seed + &uk - 1); /*eMKF: start seed numbering from &seed instead of &seed + 1 */

					data _bayeslogc_;
					run;

					/* eMKF: modified put statements */
					%put Replication _rep = &uj of &crep;
					%put Bayesian model suffix is &newuvar;
					%put Chain &uk of &chains;

					/* eMKF: workhorse for Bayesian inference -- completely revamped to use proc mcmc instead of pre-compiled C code */
					%if %upcase(&uvar) ^= BMA_CUBIC and %upcase(&uvar) ^= BMA_QUAD and %upcase(&uvar) ^= BMA_LINEAR %then 
						%bayesfit(	bdata = _bayesdata1_, blog = _bayeslogc_, btype = &uvar, 
									bgroup = _group_, btime = _time, boutcome = _y, bse = _se, 
									bn = &_ssn, brndvars = &randomVars, bARmodel = &ARmodel,
									bplot = &mcmcplot, bprint = &modelprint, bslicesampler = &slicesampler, 
									binit = &init, bprcov = &propcov, bmaxt = &maxtune, bttol= &targetaccept, batol= &accepttol,
									bseed = &_chainseed, btune = &ntu, bburn = &nbi, biter = &nmc, bthin = &thin, borpoly=&orpoly,
									bmalpha = &malpha, bpalpha = &palpha,
									bmbeta1 = &mbeta1, bpbeta1 = &pbeta1, bbeta1l = &beta1l, bbeta1u = &beta1u,
									bmbeta2 = &mbeta2, bpbeta2 = &pbeta2, bbeta2l = &beta2l, bbeta2u = &beta2u,
									bmbeta3 = &mbeta3, bpbeta3 = &pbeta3, bbeta3l = &beta3l, bbeta3u = &beta3u,
						    		bmrho = &mrho, bprho = &prho, btaul = &taul, btauu = &tauu,
									bvshape = &vshape, bvscale = &vscale, bcmploc = &cmploc
						         );;

					/* eMKF: MCMC for Bayesian model averaging uses mixture priors */
					%if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then 
						%bayesBMA(	bdata = _bayesdata1_, blog = _bayeslogc_, btype = &uvar, 
									bgroup = _group_, btime = _time, boutcome = _y, bse = _se, 
									bn = &_ssn, brndvars = &randomVars, bARmodel = &ARmodel,
									bplot = &mcmcplot, bprint = &modelprint, bslicesampler = &slicesampler,
									binit = &init, bprcov = &propcov, bmaxt = &maxtune, bttol= &targetaccept, batol= &accepttol,
									bseed = &_chainseed, btune = &ntu, bburn = &nbi, biter = &nmc, bthin = &thin, borpoly=&orpoly,
									bmalpha = &malpha, bpalpha = &palpha,
									bmbeta1 = &mbeta1, bpbeta1 = &pbeta1, bbeta1l = &beta1l, bbeta1u = &beta1u,
									bmbeta2 = &mbeta2, bpbeta2 = &pbeta2, bbeta2l = &beta2l, bbeta2u = &beta2u,
									bmbeta3 = &mbeta3, bpbeta3 = &pbeta3, bbeta3l = &beta3l, bbeta3u = &beta3u,
						    		bmrho = &mrho, bprho = &prho, btaul = &taul, btauu = &tauu,
									bvshape = &vshape, bvscale = &vscale, bwshape = &wshape, bcmploc = &cmploc
						         );;

					/* eMKF: fitc dataset holds the posterior samples: start by splitting each chain in half */
					data _bayeslogc_; 
					  set _bayeslogc_;
					  chain = &uk;
					  half = 1;
					  if mod(floor(&nmc/&thin), 2) ne 0 then do;
						if _n_ = 1 then delete; 			  /* eMKF: drop first row if chain length is not even */
					  	if _n_ > (1 + floor(&nmc/&thin))/2 then half = 2;
					  end;
					  else do;
					  	if _n_ > floor(&nmc/&thin)/2 then half = 2;
					  end;
					run;

					data _bayeslog_;
					  set _bayeslog_ _bayeslogc_;
					  if chain ne .;
					run;
					    
					/* eMKF: local cleanup */
					proc datasets nolist;
					  delete _bayeslogc_  ;
					run ;
					quit;
 					
				%end; /*End of &uk the chains */

				%put Start post-processing calculations across chains for _rep = &uj of &crep and Bayesian model suffix &newuvar;

				/**********************************************************************************************/
				/* eMKF: Post-processing calculations of group differences/disparities at the last time point */ 
				/**********************************************************************************************/

				%if &comparedto ^=%str() %then %do;

					%let _idf=0; %let _jdf=0;

					data _bayeslog_;
					  set _bayeslog_;

					  /* eMKF: summary measures */
					  groupmin = min(of eta&un._1-eta&un._&ug) ; 						/* smallest rate */
					  groupmax = max(of eta&un._1-eta&un._&ug) ; 						/* largest rate */
					  groupmn1 = (sum(of eta&un._1-eta&un._&ug) - groupmin)/(&ug - 1);	/* average of all but smallest rate */
					  groupmn2 = (sum(of eta&un._1-eta&un._&ug) - groupmax)/(&ug - 1);  /* average of all but largest rate */
					  groupdif_MRD = groupmax - groupmin;								/* maximal rate difference */
					  if abs(groupmin) > 0 then grouprat_MRR = groupmax/groupmin;		/* maximal rate ratio */
					  else grouprat_MRR = .;
					  groupdif_SRD1 = groupmn1 - groupmin;								/* summary rate difference using min as reference */
					  if abs(groupmin) > 0 then grouprat_SRR1=groupmn1/groupmin;		/* summary rate ratio using min as reference */
					  else grouprat_SRR1 = .;
					  groupdif_SRD2 = groupmax - groupmn2;								/* summary rate difference using max as reference */
					  if abs(groupmn2) > 0 then grouprat_SRR2=groupmax/groupmn2;		/* summary rate ratio using max as reference */
					  else grouprat_SRR2 = .;

					  /* eMKF: all pairwise differences -- this differs from original MKF where only the unique pairs were tracked */
					  /*       this change is made, here, for consistency with the new ratio calculations, below */
					  %do _idf=1 %to &ug;   
					  	  %do _jdf=1 %to &ug;
						  	  groupdif_&_idf._&_jdf = eta&un._&_idf - eta&un._&_jdf ;
						  %end;
					  %end;

					  /* eMKF: all pairwise ratios */
					  %do _idf=1 %to &ug;   
					  	  %do _jdf=1 %to &ug; 
							  if abs(eta&un._&_jdf) > 0 then grouprat_&_idf._&_jdf = eta&un._&_idf / eta&un._&_jdf ;
							  else grouprat_&_idf._&_jdf = . ;
						  %end;
					  %end;

					  /* eMKF: differences relative to min */
					  %do _idf=1 %to &ug;   
						  groupdif_&_idf._MIN = eta&un._&_idf - groupmin ;
					  %end;

					  /* eMKF: ratios relative to min */
					  %do _idf=1 %to &ug;   
						  if abs(groupmin) > 0 then grouprat_&_idf._MIN = eta&un._&_idf / groupmin ;
						  else grouprat_&_idf._MIN =. ;
					  %end;

					  /* eMKF: differences relative to max */
					  %do _idf=1 %to &ug;   
						  groupdif_MAX_&_idf = groupmax - eta&un._&_idf ;
					  %end;

					  /* eMKF: ratios relative to max */
					  %do _idf=1 %to &ug;   
						  if abs(eta&un._&_idf) > 0 then grouprat_MAX_&_idf = groupmax / eta&un._&_idf ;
						  else grouprat_MAX_&_idf =. ;
					  %end;

					run;
	
				%end;

				/***************************************************************************************************/
				/* eMKF: Rank-normalized and folded Gelman-Rubin split-Rhat, aka. potential scale reduction factor */ 
				/*       (see Vehtari et al 2021; DOI 10.1214/20-BA1221)                                           */
				/***************************************************************************************************/

				data _bayeslogfl_ _bayeslogflrk_ _bayeslogflrkmn_ _bayeslogflrkvr_ _bayeslogflrkbv_ _bayeslogflrkwv_ _bayeslogflrkpv_
					 			  _bayeslogrk_   _bayeslogrkmn_   _bayeslogrkvr_   _bayeslogrkbv_   _bayeslogrkwv_   _bayeslogrkpv_
					 		      			     _bayeslogmn_     _bayeslogvr_     _bayeslogbv_     _bayeslogwv_     _bayeslogpv_
					 _bayeslogGR_ 					
				  ;
				run;

				/************************* eMKF: Posterior model weights in Bayesian model averaging ****************/

				%if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then %do;

				    data _bayeslogmod_ _bayeslogmodb_ _bayeslogmodbm_ _bayeslogmodbv_ _bayeslogmodbw_;
					run;

					proc freq data=_bayeslog_ noprint;
					  by chain half;
					  tables flg /list out=_bayeslogmod_ ;
					run;

					/* eMKF: within chain */
					data _bayeslogmod_;
					  set _bayeslogmod_;
					  pweight = percent/100;
					  if mod(floor(&nmc/&thin), 2) ne 0 then vweight = (percent/100)*(1-percent/100)/((floor(&nmc/&thin)-1)/2);
					  else vweight = (percent/100)*(1-percent/100)/(floor(&nmc/&thin)/2);
					  drop count percent;
					run;

					proc sort data=_bayeslogmod_;
 				 		by flg;
					run;

					/* eMKF: mean across chains */
					proc means data=_bayeslogmod_ noprint;
					  by flg;
					  var pweight;
					  output out=_bayeslogmodbm_ mean=bpweight;
					run;

					/* eMKF: variance across chains */
					proc means data=_bayeslogmod_ noprint;
					  by flg;
					  var pweight;
					  output out=_bayeslogmodbv_ var=bvweight;
					run;

					/* eMKF: average variance within chains */
					proc means data=_bayeslogmod_ noprint;
					  by flg;
					  var vweight;
					  output out=_bayeslogmodbw_ mean=wvweight;
					run;

					/* eMKF: pooled posterior variances */
					data _bayeslogmodb_;
					  merge _bayeslogmodbm_ _bayeslogmodbv_ _bayeslogmodbw_;
					  by flg;
					  drop _TYPE_ _FREQ_;
					  if mod(floor(&nmc/&thin), 2) ne 0 then pvrweight = bvweight + (1-1/((floor(&nmc/&thin)-1)/2))*wvweight;
					  else pvrweight = bvweight + (1-1/(floor(&nmc/&thin)/2))*wvweight;
					run;

				    data _bayeslogmod_;
					  set _bayeslogmodb_;
					  SD = sqrt(pvrweight);
					  rename bpweight=weight;
					  drop bvweight wvweight pvrweight;
					run;

					/* eMKF: clean-up */
					proc datasets nolist;
					  delete _bayeslogmodb_ _bayeslogmodbm_ _bayeslogmodbv_ _bayeslogmodbw_;
					run ;
					quit;

				%end;

				/*******************eMKF: Posterior means and variances and traditional split-Rhat *****************/

				/*eMKF: Means within (split) chains */
				proc means data=_bayeslog_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
					  %if &comparedto ^=%str() %then group: ;
				  ;
				  output out=_bayeslogmn_ mean= ;
				run;

				/*eMKF: Variances within (split) chains */
				proc means data=_bayeslog_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
					  %if &comparedto ^=%str() %then group: ;
				  ;
				  output out=_bayeslogvr_ var= ;
				run;

				/*eMKF: Within-chain variances */
				proc means data=_bayeslogvr_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
					  %if &comparedto ^=%str() %then group: ;
				  ;
				  output out=_bayeslogwv_ mean= ;
				run;

				/*eMKF: Between-chain variances */
				proc means data=_bayeslogmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
				      %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
					  %if &comparedto ^=%str() %then group: ;
				  ;
				  output out=_bayeslogbv_ var= ;
				run;

				/*eMKF: Means across chains */
				proc means data=_bayeslogmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
					  %if &comparedto ^=%str() %then group: ;
				  ;
				  output out=_bayeslogmn_ mean= ;
				run;

				/*eMKF: Pooled posterior variances */
				data _bayeslogpv_;
				  set _bayeslogmn_ _bayeslogbv_ _bayeslogwv_ ;
				  drop _TYPE_ _FREQ_;
				run;
				proc transpose data=_bayeslogpv_ out = _bayeslogpv_;
				run;
				data _bayeslogpv_;
				  set  _bayeslogpv_(rename=(col1=pmn col2=bv col3=wv));
				  if mod(floor(&nmc/&thin), 2) ne 0 then pvr = bv + (1-1/((floor(&nmc/&thin)-1)/2))*wv;
				  else pvr = bv + (1-1/(floor(&nmc/&thin)/2))*wv;
				run;

				/*************eMKF: Repeat for the rank-normalized split-Rhat (with only the model parameters) ************/

				proc rank data=_bayeslog_(drop = %if &comparedto ^=%str() %then group: ;) /*eMKF: exclude disparities */
						  out=_bayeslogrk_ ties=mean normal=blom; 
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ;
					  eta:
				  ;
				run;

				proc means data=_bayeslogrk_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogrkmn_ mean= ;
				run;

				proc means data=_bayeslogrk_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogrkvr_ var= ;
				run;

				proc means data=_bayeslogrkvr_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogrkwv_ mean= ;
				run;

				proc means data=_bayeslogrkmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogrkbv_ var= ;
				run;

				proc means data=_bayeslogrkmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogrkmn_ mean= ;
				run;

				data _bayeslogrkpv_;
				  set _bayeslogrkmn_ _bayeslogrkbv_ _bayeslogrkwv_;
				  drop _TYPE_ _FREQ_;
				run;
				proc transpose data=_bayeslogrkpv_ out = _bayeslogrkpv_;
				run;
				data _bayeslogrkpv_;
				  set _bayeslogrkpv_(rename=(col1=rkpmn col2=rkbv col3=rkwv) drop=_LABEL_);
				  if mod(floor(&nmc/&thin), 2) ne 0 then rkpvr = rkbv + (1-1/((floor(&nmc/&thin)-1)/2))*rkwv;
				  else rkpvr = rkbv + (1-1/(floor(&nmc/&thin)/2))*rkwv;
				run;

				/**************eMKF: Repeat with rank-normalized split-Rhat computed on the folded draws *****************/

				/*eMKF: First, standardize selected variables relative to medians */
				proc stdize data=_bayeslog_(drop = %if &comparedto ^=%str() %then group: ;) /*eMKF: exclude disparities */
							out=_bayeslogfl_ method=median; 
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				run;

			   /*eMKF: Next, fold by applying absolute value on all numeric variables */
				data _bayeslogfl_;
				  set _bayeslogfl_;
				  array Nums[*] _numeric_;
				  do i = 1 to dim(Nums);
				    Nums[i] = abs(Nums[i]);
				  end;
				  drop i;
				run;

				/*eMKF: Now proceed as before using the folded draws */
				proc rank data=_bayeslogfl_ out=_bayeslogflrk_ ties=mean normal=blom;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				run;
				proc means data=_bayeslogflrk_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogflrkmn_ mean= ;
				run;
				proc means data=_bayeslogflrk_ noprint;
				  by chain half;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogflrkvr_ var= ;
				run;
				proc means data=_bayeslogflrkvr_ noprint; 
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogflrkwv_ mean= ;
				run;
				proc means data=_bayeslogflrkmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogflrkbv_ var= ;
				run;
				proc means data=_bayeslogflrkmn_ noprint;
				  var %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then flg;
					  %if %upcase(&randomVars) = YES %then varr: ;
					  %if %upcase(&ARmodel) = INDEP_AR %then spsi mpsi ;
					  rho: tausq: 
					  %if %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_CUBIC %then sb: mb: ;
					  a: 
					  %if %upcase(&uvar) ^= DROPPED %then b: ; 
					  eta:
				  ;
				  output out=_bayeslogflrkmn_ mean= ;
				run;

				data _bayeslogflrkpv_;
				  set _bayeslogflrkmn_ _bayeslogflrkbv_ _bayeslogflrkwv_;
				  drop _TYPE_ _FREQ_;
				run;
				proc transpose data=_bayeslogflrkpv_ out = _bayeslogflrkpv_;
				run;
				data _bayeslogflrkpv_;
				  set _bayeslogflrkpv_(rename=(col1=flrkpmn col2=flrkbv col3=flrkwv) drop=_LABEL_);
				  if mod(floor(&nmc/&thin), 2) ne 0 then flrkpvr = flrkbv + (1-1/((floor(&nmc/&thin)-1)/2))*flrkwv;
				  else flrkpvr = flrkbv + (1-1/(floor(&nmc/&thin)/2))*flrkwv;
				run;

				/******** eMKF: Combine different split-Rhats, calculate diagnostic, and create flag for poor mixing ******/

				%let uk=1;
				%let GRthreshold = %sysevalf(&GRthreshold + 0);
				data _bayeslogGR_; 
				  merge %if &comparedto ^=%str() %then _bayeslogpv_(where=(substr(_NAME_, 1, 5) ne "group")); /*eMKF: exclude disparities */
						%if &comparedto =%str() %then _bayeslogpv_; 
					    _bayeslogrkpv_ 
						_bayeslogflrkpv_;
				  if wv > 0 then splitRhat = sqrt(pvr/wv);
				  if rkwv > 0 then rksplitRhat = sqrt(rkpvr/rkwv);
				  if flrkwv > 0 then flrksplitRhat = sqrt(flrkpvr/flrkwv);
				  mixed = 1;
				  if rksplitRhat = . or flrksplitRhat = . or max(rksplitRhat, flrksplitRhat) ge &GRthreshold then do;
				  	mixed = 0;
				  	call symput("uk" , mixed);
				  end;
				  drop pmn rkpmn flrkpmn bv rkbv flrkbv;
				run;
				%let uk = %eval(&uk + 0);

				%put End post-processing calculations across chains for _rep = &uj of &crep and model suffix &newuvar;

				/* eMKF: Issue warning message for poor mixing */
				%if &uk = 0 %then %do;
					%put WARNING: Gelman-Rubin diagnostics at the threshold &GRthreshold suggest poor mixing with &chains chains;
					%put WARNING- Model predictions may be unreliable based on this threshold value ;

					proc iml;
					  print " Warning: Gelman-Rubin diagnostics at the threshold &GRthreshold suggest poor mixing with &chains chains";
					  print "  		   Model predictions may unreliable based on this threshold value: see dataset &out._bayeslogGR_";
					  print "          Consider investigating the chain-specific diagnostic plots and modifying the MCMC options ";
					quit;

			    %end;

				/* eMKF: Final MCMC estimates are the posterior means and corrected posterior variances across chains */
				data pmns pvrs;
				run;
				proc transpose data=_bayeslogpv_(drop = pvr) suffix=pmn out=pmns;
				  var pmn;
  				  ID _NAME_;
				run;
				proc transpose data=_bayeslogpv_(drop = pmn) suffix=pvr out=pvrs;
				  var pvr;
  				  ID _NAME_;
				run;
				data _bayesfit_;
				  merge pmns(drop=_NAME_) pvrs(drop=_NAME_);
				run;

				/*eMKF: param dataset will hold the model parameters + etas */
				data _bayesparam_; 
				  set _bayesfit_;
				  %if &comparedto ^=%str() %then drop group: ;
				;
				run; 

				/*eMKF: fitc dataset will hold the etas and group differences/disparities at the latest timepoint */
				data _bayesfit_; 
				  set _bayesfit_;
				  keep eta: %if &comparedto ^=%str() %then group: ;
 				;
				run;   

				/* eMKF: local cleanup */
				proc datasets nolist;
				  delete _bayeslogfl_ _bayeslogflrk_ _bayeslogflrkmn_ _bayeslogflrkvr_ _bayeslogflrkbv_ _bayeslogflrkwv_ _bayeslogflrkpv_
					     			  _bayeslogrk_   _bayeslogrkmn_   _bayeslogrkvr_   _bayeslogrkbv_   _bayeslogrkwv_   _bayeslogrkpv_
					 	 		      				 _bayeslogmn_     _bayeslogvr_     _bayeslogbv_     _bayeslogwv_     _bayeslogpv_
						 pmns pvrs
					  ;
				run ;
				quit;

				/*********************************************************************************************************/
				/* eMKF: Revert to long format instead of wide format for compatibility with remainder of code from RAND */
				/*********************************************************************************************************/

				data _bayesparamt_ _bayesfitt_;
				run;

				/*eMKF: transposed parameter data set */

				%let _igrp_ = 0; %let _t_ = 0;
				data _bayesparamt_ ;
				  set _bayesparam_;				  	  
				  %do _t_ = 1 %to &un;
				  	  %do _igrp_ = 1 %to &ug;

					  	  _group_ = &_igrp_;
					  	  _time = &_t_;

						  %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then %do;
						      flgpmn&newuvar = flgpmn;
						      flgpvr&newuvar = flgpvr;
							  drop flgpmn flgpvr;
						  %end;

						  %if %upcase(&randomVars) = YES %then %do;
						      varrpmn&newuvar = varr&_igrp_.pmn;
						      varrpvr&newuvar = varr&_igrp_.pvr;
							  drop varr&_igrp_.pmn varr&_igrp_.pvr;
						  %end;

						  %if %upcase(&ARmodel) = INDEP_AR %then %do;
							  spsipmn&newuvar = spsipmn;
							  spsipvr&newuvar = spsipvr;
							  mpsipmn&newuvar = mpsipmn;
							  mpsipvr&newuvar = mpsipvr;
							  rhopmn&newuvar = rho&_igrp_.pmn; 
							  rhopvr&newuvar = rho&_igrp_.pvr; 
							  tausqpmn&newuvar = tausq&_igrp_.pmn; 
							  tausqpvr&newuvar = tausq&_igrp_.pvr;
							  drop rho&_igrp_.pmn rho&_igrp_.pvr tausq&_igrp_.pmn tausq&_igrp_.pvr
                                   mpsipmn mpsipvr spsipmn spsipvr ;
						  %end;
	
						  %if %upcase(&ARmodel) = COMMON_AR %then %do;
							  rhopmn&newuvar = rhopmn; 
							  rhopvr&newuvar = rhopvr; 
							  tausqpmn&newuvar = tausqpmn; 
							  tausqpvr&newuvar = tausqpvr;
							  drop rhopmn rhopvr tausqpmn tausqpvr;
						  %end;

				  		  %if %upcase(&uvar) = FULL_CUBIC or %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = FULL_LINEAR %then %do;
							  sb1pmn&newuvar = sb1pmn; 
							  sb1pvr&newuvar = sb1pvr;
						  	  mb1pmn&newuvar = mb1pmn; 
							  mb1pvr&newuvar = mb1pvr; 
							  drop mb1pmn mb1pvr sb1pmn sb1pvr;
						  %end;

						  %if %upcase(&uvar) = FULL_CUBIC or %upcase(&uvar) = FULL_QUAD %then %do;
							  sb2pmn&newuvar = sb2pmn; 
							  sb2pvr&newuvar = sb2pvr;
						  	  mb2pmn&newuvar = mb2pmn; 
							  mb2pvr&newuvar = mb2pvr; 
							  drop mb2pmn mb2pvr sb2pmn sb2pvr;
						  %end;

						  %if %upcase(&uvar) = FULL_CUBIC %then %do;
							  sb3pmn&newuvar = sb3pmn; 
							  sb3pvr&newuvar = sb3pvr;
						  	  mb3pmn&newuvar = mb3pmn; 
							  mb3pvr&newuvar = mb3pvr; 
							  drop mb3pmn mb3pvr sb3pmn sb3pvr;
						  %end;

					  	  apmn&newuvar = a&_igrp_.pmn;
						  apvr&newuvar = a&_igrp_.pvr;
						  drop a&_igrp_.pmn a&_igrp_.pvr;

						  %if %upcase(&uvar) = COMMON_LINEAR or %upcase(&uvar) = COMMON_QUAD or %upcase(&uvar) = COMMON_CUBIC %then %do;
					  	  	  b1pmn&newuvar = b1pmn;
						  	  b1pvr&newuvar = b1pvr;
							  drop b1pmn b1pvr;
						  %end;

						  %if %upcase(&uvar) ^= DROPPED and %upcase(&uvar) ^= COMMON_LINEAR and 
							  %upcase(&uvar) ^= COMMON_QUAD and %upcase(&uvar) ^= COMMON_CUBIC %then %do;
					  	  	  b1arrpmn&newuvar = b1arr&_igrp_.pmn;
						  	  b1arrpvr&newuvar = b1arr&_igrp_.pvr;
							  drop b1arr&_igrp_.pmn b1arr&_igrp_.pvr;
						  %end;

						  %if %upcase(&uvar) = COMMON_QUAD or %upcase(&uvar) = COMMON_CUBIC %then %do;
					  	  	  b2pmn&newuvar = b2pmn;
						  	  b2pvr&newuvar = b2pvr;
							  drop b2pmn b2pvr;
						  %end;

						  %if %upcase(&uvar) = FULL_QUAD or %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = BMA_QUAD or
							  %upcase(&uvar) = FULL_CUBIC or %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = BMA_CUBIC %then %do;
					  	  	  b2arrpmn&newuvar = b2arr&_igrp_.pmn;
						  	  b2arrpvr&newuvar = b2arr&_igrp_.pvr;
							  drop b2arr&_igrp_.pmn b2arr&_igrp_.pvr;
						  %end;

						  %if %upcase(&uvar) = COMMON_CUBIC %then %do;
					  	  	  b3pmn&newuvar = b3pmn;
						  	  b3pvr&newuvar = b3pvr;
							  drop b3pmn b3pvr;
						  %end;

						  %if %upcase(&uvar) = FULL_CUBIC or %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = BMA_CUBIC %then %do;
					  	  	  b3arrpmn&newuvar = b3arr&_igrp_.pmn;
						  	  b3arrpvr&newuvar = b3arr&_igrp_.pvr;
							  drop b3arr&_igrp_.pmn b3arr&_igrp_.pvr;
						  %end;

						  etapmn&newuvar = eta&_t_._&_igrp_.pmn;
						  etapvr&newuvar = eta&_t_._&_igrp_.pvr;
						  drop eta&_t_._&_igrp_.pmn eta&_t_._&_igrp_.pvr;

						  output;

					  %end;
				  %end;
			    run; 

				proc sort data=_bayesparamt_ out=_bayesparam_ ;
				  by _group_ _time;
				run;

				/*eMKF: transposed predictions data set */

				%if &comparedto ^=%str() %then %do;

					proc transpose data=_bayesfit_(drop = eta:) out=_bayesfitt_;
					run;

					data _bayesfitt_;
					  set _bayesfitt_;
					  if scan(_name_,1,"_") in("groupdif", "grouprat") then do;
						  if scan(_name_,3,"_") ^= "" then do;
						  	 if scan(_name_,2,"_") ^= "MAX" then _diffgrp1_ = 1*scan(_name_,2,"_");
							 else _diffgrp1_ = - &ug;
							 if substr(scan(_name_,3,"_"), 1, length(scan(_name_,3,"_")) - 3) ^= "MIN" then 
						  	 	_diffgrp2_ = 1*substr(scan(_name_,3,"_"), 1, length(scan(_name_,3,"_")) - 3);
							 else _diffgrp2_ = 0;
						  end;
						  else do;
						  	 if substr(scan(_name_,2,"_"), 1, length(scan(_name_,2,"_")) - 3) in("MRD", "MRR") then do;
							 	_diffgrp1_ = - &ug;
							    _diffgrp2_ = 0;
							 end;
						  	 if substr(scan(_name_,2,"_"), 1, length(scan(_name_,2,"_")) - 3) in("SRD1", "SRR1") then do;
							 	_diffgrp1_ = - (&ug-1)/2;
							    _diffgrp2_ = 0;
							 end;
						  	 if substr(scan(_name_,2,"_"), 1, length(scan(_name_,2,"_")) - 3) in("SRD2", "SRR2") then do;
							 	_diffgrp1_ = - &ug;
							    _diffgrp2_ = -(&ug-1)/2;
							 end;
						  end;
					  end;

					  if substr(_name_, 6, 3) = "dif" then _disparity = "difference";
					  if substr(_name_, 6, 3) = "rat" then _disparity = "ratio";
					  if substr(_name_, 6, 3) = "min" then _disparity = "min";
					  if substr(_name_, 6, 3) = "max" then _disparity = "max";
					  if substr(_name_, 6, 3) = "mn1" then _disparity = "avgexclmin";
					  if substr(_name_, 6, 3) = "mn2" then _disparity = "avgexclmax";
					  if substr(_name_, length(_name_)- 2, 3) = "pvr" then type="pvr";
					  if substr(_name_, length(_name_)- 2, 3) = "pmn" then type="pmn";
					run;

					proc sort data=_bayesfitt_;
	 				  by _disparity _diffgrp1_ _diffgrp2_;
					run;

					data _bayesfitt_;
					  merge _bayesfitt_(where=(type="pmn") rename=(col1=pred_Bayes_&uvar))
					        _bayesfitt_(where=(type="pvr") rename=(col1=predVar_Bayes_&uvar))
					        ;
					  predSE_Bayes_&uvar = sqrt(predVar_Bayes_&uvar);
					  rename _diffgrp1_ = _group_;
					  drop type _name_;
					run;

					data _bayesfitt_;
	 				  merge _bayesfitt_(keep= _diffgrp2_ _disparity _group_) _bayesfitt_(drop= _diffgrp2_ _disparity _group_);
					  _time = &un;
					  _oldtime = _time;
					  if _diffgrp2_ > 0 then _newtime = _time + _diffgrp2_;
					  if _diffgrp2_ le 0 then _newtime = _time + &ug + 1 - _diffgrp2_;
					  if _diffgrp2_ =.  then _newtime = _time + &ug + 2 ;
					  _time = _newtime;
					run;
				
				%end;

				%let _igrp_ = 0; %let _t_ = 0;
				data _bayesfit_ ;
				  set _bayesfit_;
				  %if &comparedto ^= %str() %then drop group: ;;
				  %do _t_ = 1 %to &un;
				  	  %do _igrp_ = 1 %to &ug;
					  	  _group_ = &_igrp_;
						  _time = &_t_;
						  _oldtime = &_t_;
					  	  _newtime = &_t_;
						  pred_Bayes_&uvar = eta&_t_._&_igrp_.pmn;
						  predVar_Bayes_&uvar = eta&_t_._&_igrp_.pvr;
						  predSE_Bayes_&uvar = sqrt(eta&_t_._&_igrp_.pvr);
						  drop eta&_t_._&_igrp_.pmn eta&_t_._&_igrp_.pvr;
						  output;
					  %end;
				  %end;
			    run; 

				data _bayesfit_;
				  set _bayesfit_ %if &comparedto ^=%str() %then _bayesfitt_;;
				run;

				proc sort data=_bayesfit_;
				  by _group_ _time %if &comparedto ^=%str() %then _disparity;;
				run;

				/* eMKF: finalize datasets  */

				%if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then %do;
					data _bayeslogmod_;
					  set _bayeslogmod_;
					  Model = resolve('&newuvar'); 
					run;

					proc sort data=_bayeslogmod_ out=_bayeslogmod_;
					  by Model ;
					run;
				%end;

				data _bayeslogGR_;
				  set _bayeslogGR_;
				  Model = resolve('&newuvar'); 
				run;

				proc sort data=_bayeslogGR_ out=_bayeslogGR_;
				  by Model ;
				run;

				/*eMKF: remove differences and ratios from logfile to minimize dataset size (user could recalculate those if needed) */
				data _bayeslog_;
				  set _bayeslog_;
				  Model = resolve('&newuvar'); 
				  %if &comparedto ^=%str() %then drop group: ;; 
				run;

				proc sort data=_bayeslog_ out=_bayeslog_;
				  by Model Iteration ;
				run;

				/* eMKF: append datasets by model */
				data _bayesfit2_;
				  %if &ui = 1 %then set _bayesfit_;;
				  %if &ui > 1 %then merge _bayesfit2_ _bayesfit_;;
				  %if &ui > 1 %then by _group_ _time %if &comparedto ^=%str() %then _disparity;;;
				run;
				data _bayesparam2_;
				  %if &ui = 1 %then set _bayesparam_;;
				  %if &ui > 1 %then merge _bayesparam2_ _bayesparam_;;
				  %if &ui > 1 %then by _group_ _time;;
				run;
				data _bayeslogGR2_;
				  %if &ui = 1 %then set _bayeslogGR_;;
				  %if &ui > 1 %then merge _bayeslogGR2_ _bayeslogGR_;;
				  %if &ui > 1 %then by Model ;;
				run;
				data _bayeslog2_;
				  %if &ui = 1 %then set _bayeslog_;;
				  %if &ui > 1 %then merge _bayeslog2_ _bayeslog_;;
				  %if &ui > 1 %then by Model Iteration;;
				run;
				%if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then %do;
					data _bayeslogmod2_;
					  %if &ui = 1 %then set _bayeslogmod_;;
					  %if &ui > 1 %then merge _bayeslogmod2_ _bayeslogmod_;;
					  %if &ui > 1 %then by Model ;;
					run;
				%end;

				/* eMKF: local cleanup */
				proc datasets nolist;
				  delete _bayesparamt_ _bayesfitt_ _bayesparam_ _bayesfit_ _bayeslogGR_ _bayeslog_ 
				  		 %if %upcase(&uvar) = BMA_CUBIC or %upcase(&uvar) = BMA_QUAD or %upcase(&uvar) = BMA_LINEAR %then _bayeslogmod_;
				  ;
				run ;
				quit;

			%end; /*End of the Bayesian %do &ui*/

			%put End of Bayesian model fitting for _rep = &uj of &crep ;

			/*eMKF: Track the stratum identifier _rep */
			%if &flag1a = 1 or &flag2a = 1 or &flag3a = 1 %then %do;
				data _bayeslogmod2_;
				  set _bayeslogmod2_;
				  _rep = &uj;
				run;
			%end;
			data _bayeslogGR2_;
			  set _bayeslogGR2_;
			  _rep = &uj;
			run;
			data _bayeslog2_;
			  set _bayeslog2_;
			  _rep = &uj;
			run;
			data _bayesparam2_;
			  set _bayesparam2_;
			  _rep = &uj;
			run;
			data _bayesfit2_;
			  set _bayesfit2_;
			  _rep = &uj;
			run;

			/*eMKF: Prepare output dataset &out._bayes with differences and ratios */

			data _bayesfit2_;
			  set _bayesfit2_;
			  _time =_oldtime;
			  _idd+1;
			run;

			%if &comparedto ^=%str() %then %do;
				proc sort data=_bayesfit2_;
			  	  by _diffgrp2_;
				run;
			%end;

			data _junk_;
			run;

			proc sort data=_bayesdata1_ out= _junk_(keep= _group_ &group _rep) nodupkey;
			  by _group_;
			run;

			data _bayesfit2_;
			  merge _bayesfit2_ 
			  		%if &comparedto =%str() %then _junk_;
			  		%if &comparedto ^=%str() %then _junk_(rename=(_group_=_diffgrp2_ &group=diff_&group));
			  ;
			  %if &comparedto =%str() %then by _group_;;
			  %if &comparedto ^=%str() %then by _diffgrp2_;;
			run;

			proc sort data=_bayesfit2_;
			  by _idd;
			run;

			data _bayesfit2_;
			  set _bayesfit2_;
			  drop _idd;
			run;

			data &out._bayes;
			  set &out._bayes _bayesfit2_;
			  if _time ne .;
			  %if &comparedto ^=%str() %then if _diffgrp2_ = -(&ug-1)/2 then diff_&group = "AVGEXCLMAX";;;
			  %if &comparedto ^=%str() %then if _diffgrp2_ = 0 then diff_&group = "MIN";;;
			run;

			/*eMKF: Append remaining datasets */

			data &out._bayesparm;
			  set &out._bayesparm _bayesparam2_;
			  if _time ne .;
			run;

			data &out._bayeslogGR;
			  set &out._bayeslogGR _bayeslogGR2_;
			  if mixed ne .;
			run;

			%if &flag1a = 1 or &flag2a = 1 or &flag3a = 1 %then %do;
				data &out._bayeslogmod;
				  set &out._bayeslogmod _bayeslogmod2_;
				  if flg ne .;
				run;
			%end;

			/*eMKF: To keep file size manageable, use separate dataset for each replication instead of appending */
			%if %upcase(&mcmclog) = YES %then %do;
				data &out._bayeslog_rep&uj;
				  set _bayeslog2_;
				  if Iteration ne .;
				run;
			%end;

		%end; /* End of replications %do uj*/

		%put End of all Bayesian model fitting loops ;

		/* eMKF: local cleanup */

		%if %upcase(&mcmclog) ^= YES %then %do;
			%let uj=0;
			%do uj=1 %to &crep; 
				proc datasets nolist;
					delete &out._bayeslog_rep&uj; 
				run;
				quit;
			%end;
		%end;

		proc datasets nolist;
		  delete _junk_ _bayesdata1_ _bayesfit2_ _bayesparam2_ _bayeslogGR2_ _bayeslog2_ 
				 %if &flag1a = 1 or &flag2a = 1 or &flag3a = 1 %then _bayeslogmod2_;
		  ;
		run;
		quit;

		/* eMKF: Adding original data back for reference */

		data _junk_;
		run;
		proc sort data= _bayesdata_ out=_junk_(drop= &outcome2 &se2) nodupkey;
		  by _rep _group_ _time ;
		run;

		data _junkN_;
		run;
		proc sort data= _bayesdata_(where = (_group_ = 1 and _time = &un)) out=_junkN_(keep= &time &by _time _rep) nodupkey;
		  by _rep ;
		run;

		data _outbs1 _outbs2 _outbs3;
		run;

		data _outbs3;
		  merge _junk_ &out._bayes(where = (_newtime le &un));
		  by _rep _group_ _time;
		  %if &comparedto ^=%str() %then if _group_ = -&ug then &group = "MAX";;;
		  %if &comparedto ^=%str() %then if _group_ = -(&ug-1)/2 then &group = "AVGEXCLMIN";;;
		run;

		%if &comparedto ^=%str() %then %do;

			data _outbs1;
			  merge _junkN_(drop=_time) &out._bayes(where = (_newtime > &un));
			  by _rep;
			run;

			proc sort data= _outbs1;
			  by _rep _group_ _newtime ;
			run;

			data _outbs2;
			  merge _junk_(where=(_time = &un)) _outbs1;
			  by _rep _group_ ;
			  if _group_ = -&ug then &group = "MAX";
			  if _group_ = -(&ug-1)/2 then &group = "AVGEXCLMIN";
			run;

			proc sort data= _outbs2;
			  by _rep _group_ _newtime ;
			run;

		%end;

		data &out._bayes;
		    set _outbs2 _outbs3;
		run;

		proc sort data= &out._bayes;
		  by _rep _group_ _time _newtime;
		run;

		data &out._bayesparm;
		  merge _junk_ &out._bayesparm;
		  by _rep _group_ _time;
		run;

		/* eMKF: Set up comparisons data set if requested */

		%let _comp2=; 

		%if &comparedto ^=%str() %then %do;

			data &comparedata;
			run;

			data &comparedata;
			  set &out._bayes;
			  if _group_ = . or _diffgrp2_ ne . ;
			run;

			/*eMKF: find label of comparison group - modified to allow for min/max used as reference  */
			data _junk_ _junk0_;
			run;
			proc freq data=&comparedata noprint;
			  tables &group/ list out=_junk_;
			run;
			proc freq data=&comparedata noprint;
			  tables diff_&group/ list out=_junk0_;
			run;
			data _junk_;
			  set _junk_ _junk0_;
			  _comp_ = upcase(compress("&comparedto"));
			  if upcase(compress(&group)) = _comp_ then call symput("_comp2", &group);
			  else if upcase(compress(diff_&group)) = _comp_ then call symput("_comp2", diff_&group);;
			run;

			%if &_comp2 = %str() and %upcase(&_comp2) ^= MIN and %upcase(&_comp2) ^= MAX %then %do;
			    /*eMKF: Shortened warning message in log*/
				%put WARNING: The comparison group &comparedto is not a &group or a recognized reference value;
				%put WARNING- No comparison will be printed at this point;

				proc iml; /*eMKF: Added warning message to HTML output*/
				    print " Warning: The comparison group &comparedto is not a &group or a recognized reference value";
					print "  		 Check to make sure the value &comparedto is correct";
					print "          No comparison will be printed at this point";
                    print "          All comparisons could be found in the &out._bayes data ";
				quit;

			%end;

			proc datasets nolist;
		  		delete _junk0_ ;
			run ;
			quit;

		%end;

		%let comparedto = &_comp2;

		%if &comparedto ^=%str() %then %do;

			data &comparedata;
			run;

			data &comparedata;
			  set &out._bayes;
			  if _diffgrp2_ ne . ;
			  if upcase(compress(&group)) = upcase(compress("&comparedto")) or upcase(compress(diff_&group)) = upcase(compress("&comparedto"));
			  if upcase(compress(&group)) = upcase(compress(diff_&group)) then delete;
			  if upcase(compress("&comparedto")) ^= "MAX" and upcase(compress("&comparedto")) ^= "MIN" then do;
			  	if upcase(compress(&group)) = "MAX" or upcase(compress(diff_&group)) = "MIN" then delete;
			  end;
			  &group._1 = &group;
			  &group._2 = diff_&group;
			  _thekey = 1;
			  rename &_thekeep1b;;
			  drop _avgse _y _se &outcome &se &outcome2 &se2 impute _diffgrp2_ _oldtime _newtime predvar_Bayes_: 
			       %if &by ^=%str() %then _avgseb imputeb;
				   %if %upcase(&randomVars) = YES %then _avgn _n;
				   %if %upcase(&randomVars) = YES and &by ^=%str() %then _avgnb;
              ;
			run;

			data _junk_;
			run;

			proc sort data= &out._bayes out= _junk_(keep= &by inputorder ) nodupkey;
			  by &by descending inputorder;
			run;

			data _junk_;
			  set _junk_;
			  _thekey =1;
			run;

			proc sort data= _junk_ out= _junk_ nodupkey;
			  by _thekey &by ;
			run;

			data _junk_;
			  set _junk_;
			  inputorder = inputorder + 0.01;
			run;

			proc sort data=&comparedata;
			  by _thekey &by inputorder;
			run;

			data &comparedata;
			  merge &comparedata(drop=inputorder) _junk_;
			  by _thekey &by;
			  drop _thekey;
			run;

			proc sort data=&comparedata out=&comparedata;
			  by _disparity %if &by ^=%str() %then &by ;;
			run;

			data &comparedata;
			  set &comparedata;
			  by _disparity %if &by ^=%str() %then &by ;;
			  _mid + 1;
			  if first._disparity then _mid = 1;
			  %if &by ^=%str() %then if first.&by then _mid = 1; ;;
			  if _mid=1 then firstp=1;
			run;

		%end;

		/* eMKF: more local cleanup */
		proc datasets nolist;
		  delete _junk_ _junkN_ _outbs1 _outbs2 _outbs3;
		run ;
		quit; 

		%let uloc = &cmploc..uds;
		proc fcmp outlib=&uloc; 
			deletesubr EP;
			%if %upcase(&randomVars) = YES %then deletesubr RP;; 
			%if &flag1f = 1 %then deletesubr MP_bfc;;
			%if &flag2f = 1 %then deletesubr MP_bfq;;
			%if &flag3f = 1 %then deletesubr MP_bfl;;
			%if &flag1a = 1 			  %then deletesubr CP_bmac;; 
			%if &flag2a = 1 			  %then deletesubr CP_bmaq;; 
			%if &flag3a = 1 			  %then deletesubr CP_bmal;; 
			%if &flag1 = 1 or &flag1f = 1 %then deletesubr CP_bgc;; 
			%if &flag2 = 1 or &flag2f = 1 %then deletesubr CP_bgq;; 
			%if &flag3 = 1 or &flag3f = 1 %then deletesubr CP_bgl;; 
			%if &flag4 = 1 				  %then deletesubr CP_b1c;; 
			%if &flag5 = 1 				  %then deletesubr CP_b1q;; 
			%if &flag6 = 1 				  %then deletesubr CP_b1l;; 
			%if &flag7 = 1 				  %then deletesubr CP_b0;; 
			%if &flag1a = 1 %then deletesubr FP_bmac;; 
			%if &flag2a = 1 %then deletesubr FP_bmaq;; 
			%if &flag3a = 1 %then deletesubr FP_bmal;; 
		run;
		quit;	
		options cmplib = _null_;

	%end; /* End of the Bayesian model(s) fitting*/

    /************************************************************************************************/
	/* Finalizing the results into a single dataset (eMKF: modified to account for added variables) */
    /************************************************************************************************/

	%if &slopes =%str() and %upcase(&Bayesian)=YES %then %do;
		data &out._pred;
		  set &out._bayes;
		  %if &comparedto ^=%str() %then if _group_ ^= . and _diffgrp2_ = . ;;;
		run;
	%end;

	%if &slopes ^=%str() and %upcase(&Bayesian)=YES %then %do;
		proc sort data=&out._pred;
		  by _rep &group _time;
		run;
		proc sort data=&out._bayes;
		  by _rep &group _time;
		run;
		data &out._pred;
		  merge &out._pred 
				%if &comparedto =%str() %then &out._bayes;
				%if &comparedto ^=%str() %then &out._bayes(where=(_group_ ^= . and _diffgrp2_ = . ));
			;
		  by _rep &group _time;
		run;
		proc sort data=&out._pred;
		  by &by &group _time;
		run;
	%end;

	data &out._pred;
	  merge &out._pred(keep = &by &group &time &outcome &se &neff &outcome2 &se2 &neff2)
			&out._pred(keep= _rep)
		    &out._pred(keep= _group_ _time _rtime) 
	        &out._pred(keep= _y _se %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _y2 _se2; ) 
	        &out._pred(keep= _avgse 
                             %if &by ^=%str() %then _avgseb; 
                             %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _avgse2; 
                             %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() and &by ^=%str() %then _avgse2b; )
			%if %upcase(&randomVars) = YES and %upcase(&Bayesian)=YES %then %do;
				&out._pred(keep= _n _avgn 
                                 %if &by ^=%str() %then _avgnb;
                                 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _n2 _avgn2;
                                 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() and &by ^=%str() %then _avgn2b; )
			%end;
	        &out._pred(keep= impute inputorder %if &by ^=%str() %then imputeb;)
			&out._pred(keep= pred:) 
			%if %upcase(&toprint2) = MODELAVG %then &out._pred(keep= p:);
		;
	run;

	/* eMKF: added functionality to remind the user about the predictions dataset if not tabulated */
	%if %upcase(&finalprint) ^= YES %then %do;
	    %put ;
		%put Tabulated printout of MKF predictions for the last time point was turned off by the user;
		%put Model predictions are in dataset &out._pred;

		proc iml;
		    print " Note: The tabulated printout of MKF predictions for the last time point was turned off by the user";
			print "  	  Model predictions can always be found in dataset &out._pred";
		quit;

	%end;
	%else %do;

		data _junk_;
		run;

		%if %upcase(&Bayesian)= YES %then %do;
			%if &comparedto ^=%str() %then %do;
				proc sort data=&out._bayes(where=(_group_ ^= . and _diffgrp2_ = . )) out= _junk_;
			  		by &by _time;
				run;
			%end;
			%else %do;
				proc sort data=&out._bayes out= _junk_;
			  		by &by _time;
				run;
			%end;
		%end;
		%else %do;
			proc sort data=&out._pred out= _junk_;
			  by &by _time;
			run;
		%end;

		proc sort data=_finalprint_;
		  by &by _time;
		run;

		data _junk_;
		  merge _junk_(in=b) _finalprint_(in=a);
		  by &by _time;
		  if a and b;
		run;

		data _junk_;
		  merge _junk_(drop=pred:) _junk_(keep=pred:) ;
		run;

		proc sort data=_junk_;
		  by &by &group _time;
		run;

		data _junk_;
		  set _junk_;
		  %if %upcase(&Bayesian)= YES %then final_pred=pred_Bayes_&toprint ;;
		  %if %upcase(&Bayesian)= YES %then final_se=predSE_Bayes_&toprint ;;
		  %if %upcase(&Bayesian)^= YES %then final_pred=pred_&toprint2 ;;
		  %if %upcase(&Bayesian)^= YES %then final_se=predSE_&toprint2 ;;
		  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then final_pred2=pred2_&toprint2 ;;
		  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then final_se2=pred2SE_&toprint2 ;;  
		  keep &by &group _group_ _time &time _y _se final_pred final_se inputorder impute 
               %if &by ^=%str() %then imputeb;
		       %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then  final_pred2 final_se2  _y2 _se2 ;
		  ;
		run;

		/*eMKF: Extended to account for additional options */
		%let _thet=;
		%if %upcase(&Bayesian)  = YES  %then %do;
			%if %upcase(&toprint)=FULL_CUBIC 	%then %let _thet= Fully Bayesian cubic trends ;
			%if %upcase(&toprint)=INDEP_CUBIC 	%then %let _thet= Independent Bayesian cubic trends ;
			%if %upcase(&toprint)=COMMON_CUBIC 	%then %let _thet= Common Bayesian cubic trend ;
			%if %upcase(&toprint)=FULL_QUAD 	%then %let _thet= Fully Bayesian quadratic trends ;
			%if %upcase(&toprint)=INDEP_QUAD 	%then %let _thet= Independent Bayesian quadratic trends ;
			%if %upcase(&toprint)=COMMON_QUAD 	%then %let _thet= Common Bayesian quadratic trend ;
			%if %upcase(&toprint)=FULL_LINEAR 	%then %let _thet= Fully Bayesian linear trends ;
			%if %upcase(&toprint)=INDEP_LINEAR 	%then %let _thet= Independent Bayesian linear trends ;
			%if %upcase(&toprint)=COMMON_LINEAR %then %let _thet= Common Bayesian linear trend ;
			%if %upcase(&toprint)=DROPPED 		%then %let _thet= Bayesian model with no trends ;
			%if %upcase(&toprint)=BMA_CUBIC 	%then %let _thet= Bayesian model average up to unconstrained cubic trend ;
			%if %upcase(&toprint)=BMA_QUAD 		%then %let _thet= Bayesian model average up to unconstrained quadratic trend ;
			%if %upcase(&toprint)=BMA_LINEAR 	%then %let _thet= Bayesian model average up to unconstrained linear trend ;
		%end;
		%if %upcase(&Bayesian)  ^= YES %then %do;
			%if %upcase(&toprint2)=INDEP_CUBIC 	 %then %let _thet= Independent ML-based cubic trends ;
			%if %upcase(&toprint2)=INDEP_QUAD 	 %then %let _thet= Independent ML-based quadratic trends ;
			%if %upcase(&toprint2)=INDEP_LINEAR  %then %let _thet= Independent ML-based linear trends ;
			%if %upcase(&toprint2)=COMMON_CUBIC  %then %let _thet= Common ML-based cubic trend ;
			%if %upcase(&toprint2)=COMMON_QUAD 	 %then %let _thet= Common ML-based quadratic trend ;
        	%if %upcase(&toprint2)=COMMON_LINEAR %then %let _thet= Common ML-based linear trend ;
			%if %upcase(&toprint2)=DROPPED 		 %then %let _thet= ML-based model with no trends ;
			%if %upcase(&toprint2)=MODELAVG 	 %then %let _thet= ML-based model average of selected trend models ;
		%end;
		%if %upcase(&randomVars) = YES  and %upcase(&ARmodel) = COMMON_AR %then %let _thet = &_thet with random sampling variances and common AR parameters ;
		%if %upcase(&randomVars) ^= YES and %upcase(&ARmodel) = COMMON_AR %then %let _thet = &_thet with fixed sampling variances and common AR parameters ;
		%if %upcase(&randomVars) = YES  and %upcase(&ARmodel) = INDEP_AR  %then %let _thet = &_thet with random sampling variances and independent AR parameters ;
		%if %upcase(&randomVars) ^= YES and %upcase(&ARmodel) = INDEP_AR  %then %let _thet = &_thet with fixed sampling variances and independent AR parameters ;
		%let _thet = &_thet across &group groups ;

		data _junk_;
		  set _junk_;
		  stddiff= (final_pred - _y)/_se;
		  _rse=_se/_y;
		  bayes_rse = final_pred / final_se;
		  ratio_se= final_se /_se;
		  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
			  stddiff2= (final_pred2 - _y2)/_se2;
			 _rse2=_se2/_y2;
			 bayes_rse2 = final_pred2 / final_se2;
			 ratio_se2= final_se2 /_se2;
		  %end;
		run;

		%if %scan(&outcome2,1) =%str() or %scan(&se2,1) =%str() %then %do;
			data _junk_;
			  set _junk_;
			  array Apred(1:2) _y final_pred;
			  array Ase(1:2) _se  final_se;
			  array Arse(1:2) _rse  bayes_rse;
			  do _i_ = 1 to 2;
			  	prediction = Apred[_i_];
			   	predse = Ase[_i_];
			   	predrse = Arse[_i_];
			   	if _i_=2 then pred_stddiff = stddiff;
			   	if _i_=2 then pred_ratiose = ratio_se;
			   	output;
			  end;
			  drop _y  _se final_pred final_se stddiff ratio_se _rse bayes_rse;
			run;
		%end;
		%else %do;
			data _junk_;
			  set _junk_;
			  array Apred(1:4) _y   final_pred _y2   final_pred2;
			  array Ase(1:4)   _se  final_se   _se2  final_se2;
			  array Arse(1:4)  _rse bayes_rse  _rse2 bayes_rse2;
			  do _i_ = 1 to 4;
			   	prediction = Apred[_i_];
			   	predse = Ase[_i_];
			   	predrse = Arse[_i_];
			   	if _i_=2 then pred_stddiff = stddiff;
			   	if _i_=2 then pred_ratiose = ratio_se;
			   	if _i_=3 then pred_stddiff = .;
			   	if _i_=3 then pred_ratiose = .;
			   	if _i_=4 then pred_stddiff = stddiff2;
			   	if _i_=4 then pred_ratiose = ratio_se2;
			   	output;
			  end;
			  drop _y  _se final_pred final_se stddiff ratio_se _rse bayes_rse;
			run;
		%end;

		%if &comparedto ^=%str() and %upcase(&Bayesian) = YES %then %do;

			data _junnk_;
			run;

			data _junnk_;
			  set &comparedata;
			  %if %upcase(&toprint)=FULL_CUBIC 	  %then rename pred_bfc= prediction rmse_bfc =predse ;;
			  %if %upcase(&toprint)=INDEP_CUBIC   %then rename pred_bgc= prediction rmse_bgc =predse ;;
			  %if %upcase(&toprint)=COMMON_CUBIC  %then rename pred_b1c= prediction rmse_b1c =predse ;;
			  %if %upcase(&toprint)=FULL_QUAD 	  %then rename pred_bfq= prediction rmse_bfq =predse ;;
			  %if %upcase(&toprint)=INDEP_QUAD    %then rename pred_bgq= prediction rmse_bgq =predse ;;
			  %if %upcase(&toprint)=COMMON_QUAD   %then rename pred_b1q= prediction rmse_b1q =predse ;;
			  %if %upcase(&toprint)=FULL_LINEAR   %then rename pred_bfl= prediction rmse_bfl =predse ;;
			  %if %upcase(&toprint)=INDEP_LINEAR  %then rename pred_bgl= prediction rmse_bgl =predse ;;
			  %if %upcase(&toprint)=COMMON_LINEAR %then rename pred_b1l= prediction rmse_b1l =predse ;;
			  %if %upcase(&toprint)=DROPPED 	  %then rename pred_b0= prediction rmse_b0 =predse ;;
			  %if %upcase(&toprint)=BMA_CUBIC     %then rename pred_bmac= prediction rmse_bmac =predse ;;
			  %if %upcase(&toprint)=BMA_QUAD      %then rename pred_bmaq= prediction rmse_bmaq =predse ;;
			  %if %upcase(&toprint)=BMA_LINEAR    %then rename pred_bmal= prediction rmse_bmal =predse ;;
			run;

			data _junnk_;
			  set _junnk_;
			  _i_=5;
			  if _disparity = "difference" then _groupdiff_ = compress(&group._1)||" - "|| &group._2;
			  if _disparity = "ratio" then _groupratio_ = compress(&group._1)||" / "|| &group._2;
			  keep &by &time &group &group._1 &group._2 prediction predse _i_ _groupdiff_ _groupratio_ _disparity inputorder;
			run;

			data _junk_;
			  set _junk_ _junnk_ ;
			run;

			/* eMKF: clean up */		
			proc datasets nolist;
			  delete _junnk_ ;
			run ;
			quit;

		%end;

		data _junk_;
		  set _junk_;

		  /*eMKF: CIs formatted to account for extra space for the minus sign */
		  %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do;
			  if _disparity = "difference" or _disparity = "" then do;
				  if prediction -(1.96*predse) ge 0 and prediction +(1.96*predse) ge 0 
					then pred_ci = "[ "||compress(put(prediction -(1.96*predse), 16.&pdigit))||",  "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
				  if prediction -(1.96*predse)  < 0 and prediction +(1.96*predse) ge 0 
					then pred_ci = "["||compress(put(prediction -(1.96*predse), 16.&pdigit))||",  "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
				  if prediction -(1.96*predse) ge 0 and prediction +(1.96*predse)  < 0 
					then pred_ci = "[ "||compress(put(prediction -(1.96*predse), 16.&pdigit))||", "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
				  if prediction -(1.96*predse)  < 0 and prediction +(1.96*predse)  < 0 
					then pred_ci = "["||compress(put(prediction -(1.96*predse), 16.&pdigit))||", "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
			  end;
			  if _disparity = "ratio" and prediction > 0 then do; /* eMKF: Added for printing ratios and their lognormal CIs */
				  pred_ci = "[ "||compress(put(exp(log(prediction) -(1.96*predse/prediction)), 16.&pdigit))||",  "|| compress(put(exp(log(prediction) +(1.96*predse/prediction)), 16.&pdigit))||"]";
			  end;
		  %end;
		  %else %do;
			  if prediction -(1.96*predse) ge 0 and prediction +(1.96*predse) ge 0 
				then pred_ci = "[ "||compress(put(prediction -(1.96*predse), 16.&pdigit))||",  "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
			  if prediction -(1.96*predse)  < 0 and prediction +(1.96*predse) ge 0 
				then pred_ci = "["||compress(put(prediction -(1.96*predse), 16.&pdigit))||",  "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
			  if prediction -(1.96*predse) ge 0 and prediction +(1.96*predse)  < 0 
				then pred_ci = "[ "||compress(put(prediction -(1.96*predse), 16.&pdigit))||", "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
			  if prediction -(1.96*predse)  < 0 and prediction +(1.96*predse)  < 0 
				then pred_ci = "["||compress(put(prediction -(1.96*predse), 16.&pdigit))||", "|| compress(put(prediction +(1.96*predse), 16.&pdigit))||"]";
		  %end;

		  if _i_ in (1, 3) then label="Sample       ";
		  if _i_ in (2, 4) then label="MKF estimate   ";

		  if pred_stddiff ne . then stddiff= put(pred_stddiff, 16.&pdigit);
		  if pred_ratiose ne . then ratiose= put(pred_ratiose, 16.&pdigit);
		  if pred_stddiff = . then stddiff= "  ~~             ";
		  if pred_ratiose = . then ratiose=  "  ~~             ";
		  format prediction predse pred_stddiff pred_ratiose 16.&pdigit;
		run;

		data _junk_;
		  set _junk_;
		  ll1=0;
		  %if &by ^=%str() %then ll1 =length(compress(&by));;;
		  ll2=length(compress(&group));
		  ll3=0;
		  ll3=length(compress(_time));
		  %if %eval(0+ %_counts_(&time)) = 1 %then ll3=length(compress(&time));;;
		  ll4=length(compress(label));
		  ll5=length(compress(put(prediction, 16.&pdigit)));
		  ll6=length(compress(put(predse, 16.&pdigit)));
		  ll7=length(compress(pred_ci));
		  ll8=length(compress(stddiff));
		  ll9=length(compress(ratiose));
		run;

		data _freqg_;
		run;

		proc means data=_junk_ noprint;
		  var ll1-ll9;
		  output out=_freqg_ max=ll1-ll9;
		run;

		%local _flign;
		%let _flign=134; /*eMKF: increased line width from 94*/

		data _freqg_;
		  set _freqg_;
		  id+1;
		  space=5; /*eMKF: increased from 3*/
		  n1 = 4;
		  n2 = n1 + max(3, ll1) -1 + space;
		  n3 = n2 + max(5, ll2) -1 + space;
		  n4 = n3 + max(4, ll3) -1 + space;
		  n5 = n4 + max(12, ll4) -1 + space;
		  n6 = n5 + max(10, ll5) -1 + space;
		  n7 = n6 + max(5, ll6) -1 + space;
		  n8 = n7 + max(17, ll7) -1 + space;
		  n9 = n8 + max(7, ll8) -1 + space;
		  n0 = n9 + max(10, ll9) +1 + space;
		  if id=10 then call symput("_flign", col1);
		run;

		%let _flign =%eval(0 + &_flign);

		proc transpose data=_freqg_(keep=n0-n9) out=_freqg_;
		run;

		/* Let's set up the printing*/

		data _freqg_;
		 set _freqg_;
		 id +1;
		 %if &by ^=%str() %then if id=1 then name1=compress("@"||col1)||" &by                     ";;;
		 if id=2 then name1=compress("@"||col1)||" &group                     ";
		 %if &time  =%str() %then  if id=3 then name1=compress("@"||col1)||" _time                    ";;;
		 %if &time ^=%str() %then  if id=3 then name1=compress("@"||col1)||" &time                    ";;;
		 if id=4 then name1=compress("@"||col1)||" label                    ";
		 if id=5 then name1=compress("@"||col1)||" prediction                     ";
		 if id=6 then name1=compress("@"||col1)||" predse                     ";
		 if id=7 then name1=compress("@"||col1)||" pred_ci                     ";
		 if id=8 then name1=compress("@"||col1)||" stddiff                     ";
		 if id=9 then name1=compress("@"||col1)||" ratiose                     ";

		 %if &by ^=%str() %then if id=1 then name1a=compress("@"||col1)||" &by                     ";;;
		 if id=2 then name1a=compress("@"||col1)||" &group                     ";
		 %if &time  =%str() %then  if id=3 then name1a=compress("@"||col1)||" _time                    ";;;
		 %if &time ^=%str() %then  if id=3 then name1a=compress("@"||col1)||" &time                    ";;;
		 if id=4 then name1a=compress("@"||col1)||" '    &outcome Estimation:'                    ";

		 if id=4 then name1b=compress("@"||col1)||" '    &outcome2 Estimation:'                    ";

		 if id=4 then name2=compress("@"||col1)||" label                    ";
		 if id=5 then name2=compress("@"||col1)||" prediction                     ";
		 if id=6 then name2=compress("@"||col1)||" predse                     ";
		 if id=7 then name2=compress("@"||col1)||" pred_ci                     ";
		 if id=8 then name2=compress("@"||col1)||" stddiff                     ";
		 if id=9 then name2=compress("@"||col1)||" ratiose                     ";

		 %if &by ^=%str() %then if id=1 then name3=compress("@"||col1)||" '&by'                     ";;;
		 if id=2 then name3=compress("@"||col1)||" '   &group'                     ";
		 %if &time  =%str() %then  if id=3 then name3=compress("@"||col1)||" 'Time'                    ";;;
		 %if &time ^=%str() %then  if id=3 then name3=compress("@"||col1)||" '&time'                    ";;;
		 if id=4 then name3=compress("@"||col1)||" 'Estimation'                    ";
		 if id=5 then name3=compress("@"||col1)||" ' Point'                     ";
		 if id=6 then name3=compress("@"||col1)||" 'RMSE'                     "; /* eMKF: changed label from Std. Error to RMSE to avoid confusion */
		 if id=7 then name3=compress("@"||col1)||" '   Wald 95% CI'                     ";
		 if id=8 then name3=compress("@"||col1)||" 'Std.'                     ";
		 if id=9 then name3=compress("@"||col1)||" 'Rel.'                     ";

		 if id=4 then name4=compress("@"||col1)||" '   Type '                    ";
		 if id=5 then name4=compress("@"||col1)||" 'Estimate'                     ";
		 if id=8 then name4=compress("@"||col1)||" 'Diff'                     ";
		 if id=9 then name4=compress("@"||col1)||" 'RMSE'                     ";

		 /*eMKF: added if-clause to correct uninitialized variable &group._1 warning in original MKF macro */
		 %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do; 
			 %if &by ^=%str() %then if id=1 then name5=compress("@"||col1)||" &by                     ";;
			 if id=2 then name5=compress("@"||col1)||" &group._1                  ";
			 if id=3 then name5=compress("@"||col1)||" group0                     ";
			 if id=5 then name5=compress("@"||col1 -1)||" prediction                      ";
			 if id=6 then name5=compress("@"||col1)||" predse                     ";
			 if id=7 then name5=compress("@"||col1)||" pred_ci                     ";

			 if id=2 then name6=compress("@"||col1)||" &group._1                 ";
			 if id=3 then name6=compress("@"||col1)||" group0                    ";
			 if id=5 then name6=compress("@"||col1 -1)||" prediction                     ";
			 if id=6 then name6=compress("@"||col1)||" predse                     ";
			 if id=7 then name6=compress("@"||col1)||" pred_ci                     ";

			 %if &by ^=%str() %then if id=1 then name5a=compress("@"||col1)||" &by                     ";;
			 if id=2 then name5a=compress("@"||col1)||" &group._1                 ";
			 if id=3 then name5a=compress("@"||col1)||" group0                    ";
			 if id=5 then name5a=compress("@"||col1)||" prediction                     ";
			 if id=6 then name5a=compress("@"||col1)||" predse                     ";
			 if id=7 then name5a=compress("@"||col1)||" pred_ci                     ";

			 if id=2 then name6a=compress("@"||col1)||" &group._1                 ";
			 if id=3 then name6a=compress("@"||col1)||" group0                    ";
			 if id=5 then name6a=compress("@"||col1)||" prediction                     ";
			 if id=6 then name6a=compress("@"||col1)||" predse                     ";
			 if id=7 then name6a=compress("@"||col1)||" pred_ci                     ";

			 %if &by ^=%str() %then if id=1 then name7=compress("@"||col1)||" '&by'                     ";;
			 if id=2 then name7=compress("@"||col1)||" '    Disparity Measure'                     ";
			 if id=5 then name7=compress("@"||col1)||"'Estimate'                     "; /*eMKF: Modified column label*/
			 if id=6 then name7=compress("@"||col1)||" ' RMSE '                     ";
			 if id=7 then name7=compress("@"||col1)||" '      95% CI'                     ";

		 %end;

		run;

		%local nname1 nname1a nname1b nname2 nname3 nname4 nname5 nname6 nname5a nname6a nname7 ;

		%let nname1=; %let nname1a=; %let nname1b=; %let nname2=; %let nname3=; %let nname4=;
		%let nname5=; %let nname6=; %let nname5a=; %let nname6a=; %let nname7=;

		proc sql noprint;
		   select name1  into :nname1  separated by ' '  from _freqg_ ;
		   select name1a into :nname1a separated by ' '  from _freqg_ ;
		   select name1b into :nname1b separated by ' '  from _freqg_ ;
		   select name2  into :nname2  separated by ' '  from _freqg_ ;
		   select name3  into :nname3  separated by ' '  from _freqg_ ;
		   select name4  into :nname4  separated by ' '  from _freqg_ ;

		   /*eMKF: added if-clause to correct warnings in original MKF macro */
		   %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do;
			   select name5  into :nname5  separated by ' '  from _freqg_ ;
			   select name6  into :nname6  separated by ' '  from _freqg_ ;
			   select name5a into :nname5a separated by ' '  from _freqg_ ;
			   select name6a into :nname6a separated by ' '  from _freqg_ ;
			   select name7  into :nname7  separated by ' '  from _freqg_ ;
		  %end;
		quit;

		%if &_flign > 94 %then options  linesize=&_flign ;;;

		title "MKF &_thet for the outcome(s)"; /*eMKF: edited title for consistency */
		%if %scan(&outcome2,1)  =%str() or  %scan(&se2,1)  =%str() %then  title2 "&outcome";;;
		%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then  title2 "&outcome and &outcome2";;;

		proc sort data=_junk_;
		  by %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then _disparity; inputorder _i_;
		run;

		%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
			data _junk000_;
			run;
			data _junk000_;
			  set _junk_;
			  if _i_ in (1,3);
			  _i_=_i_-0.5;
			  keep &by &group _time _i_ inputorder &time;
			run;
			data _junk_;
			  set _junk_ _junk000_;
			run;
			proc sort data=_junk_;
			  by inputorder &by &group _time _i_;
			run;

			/* eMKF: clean up */		
			proc datasets nolist;
			  delete _junk000_ 
			        ;
			run ;
			quit;

		%end;

		/*eMKF: This is the data step used for printing the table */
		data _null_;
		   set _junk_;
		   by %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then _disparity; inputorder; /*eMKF: added _disparity to the by list */
		   file print header=newpage;
		   label  _group_="Group ID" _time="Time " 
		         label="Estimation Type" prediction="Prediction" predse="Std. Err" 
		         pred_ci="95% CI" stddiff="Standardized Difference" ratiose="Relative RMSE" ;
		   if first.inputorder then firstp=1;
		   %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do;
		   		if _disparity = "difference" then group0="- "||compress(&group._2);
				if _disparity = "ratio" then group0="/ "||compress(&group._2); 	    /*eMKF: Added display format for ratios */
		   %end;
		   if _i_=0.5 then put &nname1a;
		   %if %scan(&outcome2,1)  =%str() or  %scan(&se2,1)  =%str() %then if _i_=1 then put &nname1;;;
		   %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _i_=1 then put &nname2;;;
		   if _i_=2 then put &nname2;
		   if _i_=2.5 then put &nname1b;
		   if _i_=3 then put &nname2;
		   if _i_=4 then put &nname2;
		   %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do; /*eMKF: added this if-clause to correct warnings in original MKF macro */
			   if _i_=5 and first.inputorder and _disparity = "difference" 
				   then put "                         Differences between MKF point estimates by &group            ";
			   if _i_=5 and first.inputorder and _disparity = "ratio" 
				   then put "                         Ratios between MKF point estimates by &group             "; /*eMKF: Added for ratios */
			   if _i_=5 and first.inputorder then put "                        ";  
			   if _i_=5 and first.inputorder then put &nname7;   
			   if _i_=5 and firstp  = 1 and prediction < 0 then put &nname5;
			   if _i_=5 and firstp ne 1 and prediction < 0 then put &nname6;
			   if _i_=5 and firstp  = 1 and prediction ge 0 then put &nname5a;
			   if _i_=5 and firstp ne 1 and prediction ge 0 then put &nname6a;
		   %end;
		   if last.inputorder and impute=1 then put "   Warning:  For this group, user supplied SE=0 were set to average of nonzero values across timepoints";
		   %if &by ^=%str() %then if last.inputorder and imputeb=1 then put "   Warning:  For this group, user supplied SE=0 were set to average of nonzero values across strata";;;
		   if last.inputorder and _i_ ne 5 then put &_flign.*'-';;;;
		   %if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then if last.inputorder and _i_=5 then put &_flign.*'-';;;;
		   /* %if &by ^=%str() %then if last.&by or last._group_ then put &_flign.*'-';;; */
		   /* %if &by  =%str() %then if last._group_ then put &_flign.*'-';;; */
		   return;
		   newpage:
		      if _i_ ne 5 then do;
 			    put &nname3;
			    put &nname4;
		        put &_flign.*'#';
		        return;
			  end;
			  else do; 		/*eMKF: To do: would be nice to add by group label on page 2*/
				put &nname7;
				return;
			  end;
		run;

		title ;

		/* eMKF: clean up */		
		proc datasets nolist;
		  delete _junk_ _freqg_ _finalprint_ ;
		run ;
		quit;
	
	%end; /*eMKF: end if finalprint = YES */
	
%end; /*End of &run1 &run2 */

%if &comparedto ^=%str() and %upcase(&Bayesian)  = YES %then %do;

	data &comparedata; /*eMKF: Added ratios */
	  set &comparedata;
	  if _disparity = "difference" then _measure = compress(&group._1)|| " - " ||compress(&group._2);
	  if _disparity = "ratio" then _measure = compress(&group._1)|| " / " ||compress(&group._2); 
	  drop _group_ _rep _time inputorder _mid firstp;
	run;

	data &comparedata; /*eMKF: Kept ratios */
	  merge &comparedata(keep= _disparity _measure) &comparedata(keep= &group._1 &group._2) 
	        &comparedata(keep= &by &time) &comparedata(keep= pred_: rmse_:)
	  ;
	run;

%end;

data &out;
run;
data &out;
 set &out._pred;
 keep &_thekeep2;
run;
data &out;
 merge &_thekeep3;
run;
data &out;
 set &out;
 rename &_thekeep1;
run;
data &out;
 set &out;
 rename &_thekeeps;
run;

%if %upcase(&Bayesian)=YES %then %do;
	data &out._bayes;
	 set &out._bayes;
	 rename &_thekeep1b;
	 drop predVar_Bayes_: ;
	run;
	data &out._bayes;
	 set &out._bayes;
	 rename &_thekeepsb;
	run;
%end;

/* eMKF: clean up */		
proc datasets nolist;
  delete %if %upcase(&Bayesian) ^= YES %then _nlmixdata_; _bayesdata_ ;
run ;
quit;

%mend;



data _null_;
run;

/* eMKF: BAYESFIT completely overhauled to use SAS PROC MCMC for the Bayesian estimations instead of precompiled C code
 bdata              : Name of the data to be used
 blog               : Name of the output data containing full set of &biter/&bthin posterior draws
 btype              : full_cubic, full_quad, full_linear, indep_cubic, indep_quad, indep_linear, common_cubic, common_quad, common_linear, or dropped
 bgroup             : Group variable in the dataset 
 btime              : Time variable in the dataset 
 boutcome           : Outcome of interest variable in the dataset 
 bse                : Standard error variable in the dataset 
 bn				    : Effective sample size variable in the dataset (if applicable)
 brndvars			: YES if variances should be modeled; NO if variances should be assumed known
 bARmodel			: common_ar if AR parameters are common across groups; indep_ar if they are independently drawn from a common prior
 bslicesampler		: YES to use the slice sampler instead of MH algorithm for parameters that are not included in the Gibbs sampling step 
					  Default is NO due to heavier computational load.
 bseed              : random number generating seed that will allow the user to reproduce the same results in the Bayesian model
 bprcov				: method used in constructing initial covariance matrix for the MH algorithm (see proc mcmc documentation)
					  If empty, proc mcmc default of IND will be used.
 binit				: Option for generating initial values for the parameters (see documentation and leave empty to apply proc mcmc default)
					  eMKF default is REINIT to reset chains after tuning at the values set by the user
 bmaxt				: maximum number of proposal tuning loops (if empty, proc mcmc default of 24 is used; if 0, tuning will be skipped)
 batol				: Tolerance for acceptance probabilities (if empty, proc mcmc default of 0.075 is used in bttol +|- batol)
 bttol				: Target acceptance rate for random walk Metropolis. If empty, proc mcmc defaults are used, as follows: 
					  0.45 for models with 1 parameter, 0.35 for 2-4 parameters, and 0.234 for models with 5+ parameters.
 btune				: number of tuning iterations to use in each MCMC proposal tuning phase (if empty, proc mcmc default of 500 is used)
 bburn              : number of burn-in MCMC iterations (if empty, proc mcmc default of 1000 is used)
 biter              : number of post-burn-in MCMC iterations (if empty, proc mcmc default of 1000 is used)
 bthin				: controls thinning rate (if empty, proc mcmc default of 1 is used)
 borpoly  			: YES (default) for pre-transforming the design matrix using SAS IML orpol function. NO for "raw" polynomials.
          			  If YES, regression coefficients will be reverse-transformed prior to macro end. 
					  However, prior values below are assumed to be for the coefficients of the orthogonal polynomial regression if borpoly=YES.
 bmalpha , bpalpha 	: prior mean and precision for alphas
 bmbeta1			: prior mean for mean linear coefficient across groups -- used for cubic, quadratic, or linear trends (full_, indep_, and common_)
 bpbeta1 			: prior precision for mean linear coefficient across groups -- used only for SD hyperprior(s) in full_cubic, full_quad, or full_linear
 bmbeta2 			: prior mean for mean quadratic coefficient across groups -- used for cubic or quadratic trends (full_, indep_, and common_)
 bpbeta2			: prior precision for mean quadratic coefficient across groups -- used for SD hyperprior(s) in full_cubic or full_quad
 bmbeta3 			: prior mean for mean cubic coefficient across groups -- used for cubic trends (full_, indep_, and common_)
 bpbeta3			: prior precision for mean cubic coefficient across groups -- only used for SD hyperprior(s) in full_cubic
 bbeta1l, bbeta1u	: bounds for U(a,b) prior for SD of linear coefficients across groups -- only used for hyperprior(s) in full_cubic, full_quad, or full_linear
 bbeta2l, bbeta2u	: bounds for U(a,b) prior for SD of quadratic coefficients across groups -- only used for hyperprior(s) in full_cubic or full_quad
 bbeta3l, bbeta3u	: bounds for U(a,b) prior for SD of cubic coefficients across groups -- only used for hyperprior(s) in full_cubic
 bmrho, bprho		: prior mean and precision for transformed rho -- ie., psi = ln[(1-rho)/(1+rho)]
 btaul, btauu		: bounds for U(a,b) prior for tau (SD of innovation variance tausq)
 bvshape			: Shape parameter for inverse gamma prior distribution of the variance (when applicable) 
 bvscale			: Scale parameter for inverse gamma prior distribution of the variance (when applicable)
 bprint				: If YES, posterior parameter estimates and default chain-specific convergence diagnostics are printed (default is NO)
 bplot				: If YES, trace/diagnostics plots from proc mcmc will be included (default is NO)
 bcmploc			: location of the CMP library (usually set in parent macro mkf)

*/

%macro bayesfit(
             bdata	= , 
			 blog	= ,
			 btype	= full_linear, 
	   /* eMKF: Variable labels assumed to have been reformatted using macro reformat */
			 bgroup	= _group_, 
			 btime	= _time, 
			 boutcome= _y, 
			 bse	= _se,
			 bn 	= ,
			 brndvars = NO,
			 bARmodel = common_ar,
			 bslicesampler = NO,
	   /* eMKF: MCMC tuning parameters: if missing, proc mcmc defaults will be used */
			 bseed	= ,
			 bprcov = ,
			 binit  = reinit,
			 bmaxt  = ,
			 batol 	= ,	
			 bttol 	= ,
			 btune	= ,			
			 bburn  = ,
			 biter  = ,
			 bthin 	= ,
			 borpoly = YES,
	   /* eMKF: Model parameters: if missing, the data will be used to generate starting values*/
			 bmalpha  = , 		
			 bpalpha  = ,
			 bmbeta1  = 0,      /*eMKF: Constant c3 or c7 in RAND's MKF User's Guide */
			 bpbeta1  = ,
			 bmbeta2  = 0,
			 bpbeta2  = ,
			 bmbeta3  = 0,
			 bpbeta3  = ,
			 bbeta1l  = 0,		/*eMKF: Constant c5 in RAND's MKF User's Guide */
			 bbeta1u  = ,
			 bbeta2l  = 0,
			 bbeta2u  = ,
			 bbeta3l  = 0,
			 bbeta3u  = ,
             bmrho    = 0,		/*eMKF: Constant c9  in RAND's MKF User's Guide */
			 bprho    = 1,		/*eMKF: Constant c10 in RAND's MKF User's Guide */
			 btaul    = 0.0001,	/*eMKF: Constant c11 in RAND's MKF User's Guide */
			 btauu    = ,
			 bvshape   = ,
			 bvscale   = ,
	    /* eMKF: Printing and diagnostic plots are off by default */
			 bprint   = NO,
			 bplot 	  = NO,
			 bcmploc = work.funcs
               );
 
%local g n p brtm _brtimess brangeY bqrangeV bmedianV 
	   formatted dsop dscl _i _j oPPmat
       b1line b2line b3line vline etaarrline etamnarrline tauparline psiparline tausqparline rhoparline
       parline aparline vparline mbparline sbparline udsparline tauparline2 psiparline2 sbparline2
       plinea plineb1 plineb2 plineb3 plinev plinetau plinepsi 
       hplinemb1 hplinesb1 hplinemb2 hplinesb2 hplinemb3 hplinesb3 hplinempsi hplinespsi
	   initlinea initlineb1 initlineb2 initlineb3 initlinevarr initlinetau initlinepsi
	   initlinemb1 initlinemb2 initlinemb3 initlinesb1 initlinesb2 initlinesb3 
       monitorline optionline udsline rcXline rcNline initmbeta Narrline;

/* eMKF: Data assumed to have been pre-formatted using macro reformat: check and reformat if not */

%let formatted = 0;

%let dsop = %sysfunc(open(&bdata));
%if &dsop ne 0 %then %do;
	%if %sysfunc(varnum(&dsop, inputorder)) ne 0 and %sysfunc(varnum(&dsop, &btime)) ne 0 %then %let formatted = 1;
%end; 
%let dscl = %sysfunc(close(&dsop));

%let formatted = %eval(&formatted + 0);

data _bbdata_ _bbdata1_;
run;

%if &formatted = 1 %then %do;
	data _bbdata_;
	  set &bdata;
	run;
%end;
%else %do;
    %put ;
	%put Reformatting data prior to Bayesian estimation;
	%if %upcase(&brndvars) = YES and &bn = %str() %then %do;
		%put ERROR: (Effective) sample sizes bn must be specified to fit random sampling variances;
		proc iml;
			print "  Error Note:";
			print "  (Effective) sample sizes bn must be specified to fit random sampling variances. ";
		quit;
		%return;
	%end;
	%reformat(data=&bdata, 
		      outcome=&boutcome, se=&bse, neff=&bn, outcome2=, se2=, neff2=,
			  group=&bgroup, time=&btime, by=, randomVars = &brndvars,
 			  outformat= _bbdata_ );
%end;

/*eMKF: Sort by replications, group, and time */
proc sort data= _bbdata_;
  by _rep _group_ _time ;
run;

/* eMKF: Macro variable for the number of groups */
%let g=;
data _bfreqg_;
run;
proc freq data=_bbdata_ noprint;
 tables _group_ /list out=_bfreqg_;
run;
data _bfreqg_;
 set _bfreqg_;
 _grp_ +1;
 call symput('g',_grp_);
 keep _grp_ _group_;
run;

/* eMKF: Macro variable for the number of time points */
%let n=;
data _bfreqn_;
run;
proc freq data=_bbdata_ noprint;
 tables _rtime /list out=_bfreqn_;
run;
data _bfreqn_;
 set _bfreqn_;
 _tm +1;
 call symput('n',_tm);
 keep _tm _rtime;
run;

/*eMKF: Macro variable for the real times to use in calculations */
%let _brtimess = ;
data _bfreqn_;
  set _bfreqn_;
  retain _rts;
  if _n_= 1 then _rts = cat(_rtime);
  else _rts = catx(" ", _rts, _rtime);
  call symput('_brtimess', _rts);
  drop _rts;
run;

/* eMKF: variable that will be used for real time in case times are irregular */
%let brtm  = _rtime;

/* eMKF: Set numerical values to use in code */
%let n=%eval(0+&n);
%let g=%eval(0+&g);

/* eMKF: Compute variances */
data _bbdata_;
  set _bbdata_ ;
  _var = _se**2;
run;

/* eMKF: Modified to use orthogonal cubic polynomial design matrix */

data _oXmat_ _oPmat_;
run;

%if %upcase(&borpoly) = YES %then %do;

	proc iml;
	  x = { &_rtimess };	
	  x = T(x);							/* eMKF: column vector with real times */
	  oP = orpol(x, 3);					/* eMKF: orthonormal design matrix for cubic orthogonal polynomials */
	  x0 = { %cnstss(1, &n) };
	  x0 = T(x0);
	  x1 = x;
	  x2 = x#x;
	  x3 = x#x2;
	  uP = x0 || x1 || x2 || x3;		/* eMKF: raw/unstandardized design matrix */
	  oP1 = inv(T(uP)*uP)*T(uP)*oP[,1];
      oP2 = inv(T(uP)*uP)*T(uP)*oP[,2];
      oP3 = inv(T(uP)*uP)*T(uP)*oP[,3];
      oP4 = inv(T(uP)*uP)*T(uP)*oP[,4];
	  oPP = oP1 || oP2 || oP3 || oP4;	/* eMKF: right multiplication of raw uP with oPP produces orthonormal oP */
	  y = T(do(1, &n, 1));				/* eMKF: column vector of consecutive time indices */
	  yP = y || oP;
	  create _oXmat_ from yP [ colname = {"_time" "&brtm.0" "&brtm.1" "&brtm.2" "&brtm.3"} ] ;
	  append from yP; close _oXmat_;
	  create _oPmat_ from oPP [ colname = {"t0" "t1" "t2" "t3"} ] ;
	  append from oPP; close _oPmat_;
	quit;

	proc sort data=_bbdata_;
	  by _time;
	run;

	data _bbdata_;
	  merge _bbdata_ _oXmat_;
	  by _time;
	run;

	proc sort data= _bbdata_;
	  by _rep _group_ _time ;
	run;

%end;
%else %do;								/* eMKF: Add raw predictor variables */
	data _bbdata_;
	  set _bbdata_;
	  &brtm.0 = 1;
	  &brtm.1 = &brtm;
	  &brtm.2 = &brtm**2;
	  &brtm.3 = &brtm**3;
	run;
%end;

/* eMKF: Evaluate range of the data to use in setting prior parameters, as in MKF */
%let brangeY=;
data _bbjunk;
run;
proc means data=_bbdata_ noprint;
  var _y;
  output out=_bbjunk range=range;
run;
data _null_;
 set _bbjunk;
 call symput("brangeY", range);
run;
%let brangeY = %sysevalf(&brangeY + 0);

/*******************************************************************/
/* eMKF: Set any prior parameters not already provided by the user */
/*******************************************************************/

/* eMKF: c1 in RAND's MKF User's Guide */
%if &bmalpha = %str() %then %let bmalpha = %sysevalf(0.5 * &brangeY);;

/* eMKF: 1/c2 in RAND's MKF User's Guide */	
%if &bpalpha = %str() %then %let bpalpha = %sysevalf(0.000001/(&brangeY**2));; 

/* eMKF: c3 or c7 in RAND's MKF User's Guide */
%if &bmbeta1 = %str() %then %let bmbeta1 = %sysevalf(0);;

/* eMKF: 1/c4 in RAND's MKF User's Guide */ 
%if &bpbeta1  = %str() and (%upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR) 
	%then %let bpbeta1  = %sysevalf(10/(&brangeY**2));;

/* eMKF: 1/c8 in RAND's MKF User's Guide */	
%if &bpbeta1  = %str() and %upcase(&btype) ^= FULL_CUBIC and %upcase(&btype) ^=FULL_QUAD and %upcase(&btype) ^=FULL_LINEAR 
	%then %let bpbeta1  = %sysevalf(0.000001/(&brangeY**2));;

/* eMKF: c5 in RAND's MKF User's Guide  */
%if &bbeta1l  = %str() %then %let bbeta1l = %sysevalf(0);;						

/* eMKF: c6 in RAND's MKF User's Guide  */	
%if &bbeta1u  = %str() %then %let bbeta1u = %sysevalf(0.5 * &brangeY);;		

/* eMKF: c9 in RAND's MKF User's Guide  */	
%if &bmrho    = %str() %then %let  bmrho  = %sysevalf(0);;	

/* eMKF: c10 in RAND's MKF User's Guide  */	
%if &bprho    = %str() %then %let  bprho  = %sysevalf(1);;	

/* eMKF: c11 in RAND's MKF User's Guide  */	
%if &btaul    = %str() %then %let  btaul  = %sysevalf(0.0001);;	

/* eMKF: c12 in RAND's MKF User's Guide  */	
%if &btauu    = %str() %then %let  btauu  = %sysevalf(0.1 * &brangeY);;			

/* eMKF: Set cubic and quad precisions so that the coefficients tend to be smaller in magnitude as the degree increases */
%if &bmbeta2 = %str() %then %let bmbeta2 = %sysevalf(0);; 
%if &bpbeta2 = %str() %then %let bpbeta2 = %sysevalf(2.0 * &bpbeta1);; 
%if &bmbeta3 = %str() %then %let bmbeta3 = %sysevalf(0);; 	
%if &bpbeta3 = %str() %then %let bpbeta3 = %sysevalf(4.0 * &bpbeta1);; 		
%if &bbeta2l = %str() %then %let bbeta2l = %sysevalf(0);; 						
%if &bbeta2u = %str() %then %let bbeta2u = %sysevalf(1.5 * &bbeta1u);;
%if &bbeta3l = %str() %then %let bbeta3l = %sysevalf(0);; 	
%if &bbeta3u = %str() %then %let bbeta3u = %sysevalf(2.0 * &bbeta1u);; 

/***************************************************************************************/
/* eMKF: Use data to inform prior parameters for variances in the random variance case */
/***************************************************************************************/

%if %upcase(&brndvars) = YES %then %do;
	%let bqrangeV=0; %let bmedianV=0;
	data _bbjunk;
	run;
	proc means data=_bbdata_ noprint;
	  var _var;
	  output out=_bbjunk median=median qrange=qrange;
	run;
	data _null_;
	 set _bbjunk;
	 call symput("bqrangeV", qrange);
	 call symput("bmedianV", median);
	run;
	%let bqrangeV = %sysevalf(&bqrangeV + 0);
	%let bmedianV = %sysevalf(&bmedianV + 0);
	/* eMKF: Use median for mean and 10 times IQR for standard deviation of sampling variances (inverse gamma prior) */
	%if &bvshape = %str() %then %let bvshape = %sysevalf(2 + ( &bmedianV**2 / ((10 * &bqrangeV)**2) ) );;
	%if &bvscale = %str() %then %let bvscale = %sysevalf((&bvshape - 1)*&bmedianV);;
%end;
%else %do;
	%let bvshape =; 
	%let bvscale =;
%end;
 
/*************************************************************/
/* eMKF: Symbolic array declarations (resolved in proc mcmc) */
/*************************************************************/

/* eMKF: Named 1-dimensional arrays of regression parameters other than intercept (if any) */
%let b1line=; %let b2line=; %let b3line=; 
%if %upcase(&btype) ^= COMMON_CUBIC and %upcase(&btype) ^= COMMON_QUAD and %upcase(&btype) ^= COMMON_LINEAR and %upcase(&btype) ^=DROPPED 
	%then %let b1line = array b1arr[&g] b1arr1-b1arr&g ;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = INDEP_QUAD
	%then %let b2line = array b2arr[&g] b2arr1-b2arr&g ;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = INDEP_CUBIC
	%then %let b3line = array b3arr[&g] b3arr1-b3arr&g ;

/* eMKF: Constant arrays of hyperparameters to pass to UDS (fully Bayesian setting only) */
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR
	%then %let b1line = &b1line%str(;) array mb1hyp[2] (&bmbeta1 &bpbeta1) ;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD
	%then %let b2line = &b2line%str(;) array mb2hyp[2] (&bmbeta2 &bpbeta2) ;
%if %upcase(&btype) = FULL_CUBIC
	%then %let b3line = &b3line%str(;) array mb3hyp[2] (&bmbeta3 &bpbeta3) ;

/* eMKF: Named 1-dimensional arrays of unobserved true states and their means 
  (consistent with internal SAS names for random effects in proc mcmc) */
%let etamnarrline = array etamnarr[%eval(&g*&n)];
%let etaarrline   = array etaarr[%eval(&g*&n)];
%let _i = 0; %let _j = 0; 
%do _i = 1 %to &g;
   %do _j = 1 %to &n; 
		%let etamnarrline = &etamnarrline etamn&_j._&_i;
		%let etaarrline   = &etaarrline eta&_j._&_i;
   %end;
%end;

/* eMKF: Named 1-dimensional array of random sampling variances (if applicable) */
%let vline=;
%if %upcase(&brndvars) = YES %then %do;
	%let vline = array varr[&g] varr1-varr&g ; 
	%let vline = &vline%str(;) array vhyp[2] (&bvshape &bvscale) ; /* add array of hyperparameters to pass to UDS */
%end;

/* eMKF: Dynamic array of effective sample sizes (if applicable) to use with read_array */
%let Narrline=;
%if %upcase(&brndvars) = YES %then %let Narrline = array Narr[1] /nosymbols ;

/*****************************************************************/
/* eMKF: Symbolic parameter declarations (resolved in proc mcmc) */
/*****************************************************************/

/* eMKF: Slice sampler, if requested, would apply to selected parameters for which Gibbs sampling is not available */
%let bslice =%str(;) ;
%if %upcase(&bslicesampler) = YES %then %let bslice = %str(/slice ;);

/* eMKF: Group-specific AR parameters (if applicable) */
%let tauparline=; %let psiparline=; %let tausqparline=; %let rhoparline=; %let tauparline2 =; %let psiparline2 =;
%if %upcase(&bARmodel) = INDEP_AR %then %do;
  	%let psiparline2 = &psiparline2 parms spsi &bslice;			/* SD hyperparameter for mean of psi */
  	%let psiparline2 = &psiparline2 parms mpsi &bslice;			/* mean hyperparameter for mean of psi */
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let psiparline   = &psiparline psi&_i ; 
		%let tauparline   = &tauparline tau&_i ; 
		%let tausqparline = &tausqparline tausq&_i ; 
		%let rhoparline   = &rhoparline rho&_i ; 
    	%let psiparline2  = &psiparline2 parms psi&_i &bslice;	/* Group-specific psi1 through psi&g  */
    	%let tauparline2  = &tauparline2 parms tau&_i &bslice;	/* Group-specific innovation SDs tau1 through tau&g */
	%end;
%end;
%if %upcase(&bARmodel) = COMMON_AR %then %do;
    %let psiparline2 = parms psi &bslice;						/* Common psi = ln[(1-rho)/(1+rho)] */
    %let tauparline2 = parms tau &bslice;    					/* Common innovation SD tau */
%end;

/* eMKF: Hyper-parameters in the full Bayesian models */
%let mbparline=; %let sbparline=; %let sbparline2=;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR %then %do;
	%let mbparline = &mbparline mb1 ;
	%let sbparline = &sbparline sb1 ;
	%let sbparline2 = &sbparline2 parms sb1 &bslice;
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD %then %do;
	%let mbparline = &mbparline mb2 ;
	%let sbparline = &sbparline sb2 ;
	%let sbparline2 = &sbparline2 parms sb2 &bslice;
%end;
%if %upcase(&btype) = FULL_CUBIC %then %do;
	%let mbparline = &mbparline mb3 ;
	%let sbparline = &sbparline sb3 ;
	%let sbparline2 = &sbparline2 parms sb3 &bslice;
%end;

/* eMKF: Intercepts */
%let aparline=; %let _i = 0;
%do _i = 1 %to &g; 
	%let aparline = &aparline a&_i ; 
%end; 

/* eMKF: Linear, quadratic, and cubic coefficients, as needed */
%let parline = &aparline; %let _i = 0;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = INDEP_CUBIC %then %do;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i b2arr&_i b3arr&_i ;
	%end;
%end;
%if %upcase(&btype) = FULL_QUAD or %upcase(&btype) = INDEP_QUAD %then %do;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i b2arr&_i ;
	%end;
%end;
%if %upcase(&btype) = FULL_LINEAR or %upcase(&btype) = INDEP_LINEAR %then %do;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i ;
	%end;
%end;
%if %upcase(&btype) = COMMON_CUBIC   %then %let parline = &parline b1 b2 b3 ;
%if %upcase(&btype) = COMMON_QUAD    %then %let parline = &parline b1 b2 ;
%if %upcase(&btype) = COMMON_LINEAR  %then %let parline = &parline b1 ;

/* eMKF: Variance parameters (if applicable) */
%let vparline = ; %let _i = 0;
%if %upcase(&brndvars) = YES %then %do;
	%do _i = 1 %to &g; 
		%let vparline = &vparline varr&_i ;
	%end;
%end;

/*************************************/
/* eMKF: UDS parameters declarations */
/*************************************/
%let udsparline = ;

/* mbetag and Dbetag updated with the hyper-parameters in the fully Bayesian models */
%if &mbparline ^= %str() %then %let udsparline = &udsparline parms &mbparline mbetag Dbetag %str(/uds ;);

/* etamnarr updated with the regression coefficients */
%let udsparline = &udsparline parms &parline etamnarr %str(/uds ;);

/* true states updated in a separate UDS block */
%let udsparline = &udsparline parms etaarr %str(/uds ;);

/* variances updated in a separate UDS block (when applicable) */
%if &vparline ^= %str() %then %let udsparline = &udsparline parms &vparline %str(/uds ;);

/**************************************************************************/
/* eMKF: Symbolic prior/hyperprior specifications (resolved in proc mcmc) */
/**************************************************************************/

/* eMKF: Priors for AR parameters */
%let plinetau=; %let plinepsi=; %let hplinempsi=; %let hplinespsi=;
%if %upcase(&bARmodel) = COMMON_AR %then %do; 							/* common AR parameters */
	%let plinepsi = prior psi ~ normal(&bmrho, prec=&bprho); 	
	%let plinetau = prior tau ~ uniform(&btaul, &btauu);     	
%end;
%if %upcase(&bARmodel) = INDEP_AR %then %do; 							/* group-specific AR parameters */
	%let hplinespsi = hyperprior spsi ~ uniform(0.0001,sqrt(1/&bprho)); /* Keep away from zero */
	%let hplinempsi = hyperprior mpsi ~ normal(&bmrho, prec=&bprho);
	%let plinepsi = prior &psiparline ~ normal(mpsi, sd=spsi);
	%let plinetau = prior &tauparline ~ uniform(&btaul, &btauu);   	
%end;

/* eMKF: Hyper-prior specification in the full Bayesian models */
%let hplinemb1=; %let hplinemb2=; %let hplinemb3=; %let hplinesb1=; %let hplinesb2=; %let hplinesb3=;
%if %upcase(&btype) = FULL_CUBIC %then %do;
	%let hplinesb3 = hyperprior sb3 ~ uniform(&bbeta3l, &bbeta3u);
	%let hplinemb3 = hyperprior mb3 ~ normal(&bmbeta3, prec=&bpbeta3);
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD %then %do;
	%let hplinesb2 = hyperprior sb2 ~ uniform(&bbeta2l, &bbeta2u);
	%let hplinemb2 = hyperprior mb2 ~ normal(&bmbeta2, prec=&bpbeta2);
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR %then %do;
	%let hplinesb1 = hyperprior sb1 ~ uniform(&bbeta1l, &bbeta1u);
	%let hplinemb1 = hyperprior mb1 ~ normal(&bmbeta1, prec=&bpbeta1); 
%end;

/* eMKF: Prior for intercepts a1 through a&g */
%let plinea = prior &aparline ~ normal(&bmalpha, prec=&bpalpha);

/* eMKF: Priors for regression parameters other than intercepts */
%let plineb1=; %let plineb2=; %let plineb3=;
%if %upcase(&btype) = INDEP_CUBIC %then
	%let plineb3 = prior b3arr: ~ normal(&bmbeta3, prec=&bpbeta3);
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = INDEP_QUAD %then 
	%let plineb2 = prior b2arr: ~ normal(&bmbeta2, prec=&bpbeta2);
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = INDEP_QUAD or %upcase(&btype) = INDEP_LINEAR %then 
	%let plineb1 = prior b1arr: ~ normal(&bmbeta1, prec=&bpbeta1);
%if %upcase(&btype) = COMMON_CUBIC %then 
	%let plineb3 = prior b3 ~ normal(&bmbeta3, prec=&bpbeta3);
%if %upcase(&btype) = COMMON_CUBIC or %upcase(&btype) = COMMON_QUAD %then 
	%let plineb2 = prior b2 ~ normal(&bmbeta2, prec=&bpbeta2);
%if %upcase(&btype) = COMMON_CUBIC or %upcase(&btype) = COMMON_QUAD or %upcase(&btype) = COMMON_LINEAR %then 
	%let plineb1 = prior b1 ~ normal(&bmbeta1, prec=&bpbeta1);
%if %upcase(&btype) = FULL_CUBIC %then 
	%let plineb3 = prior b3arr: ~ normal(mb3, sd=sb3);
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD %then 
	%let plineb2 = prior b2arr: ~ normal(mb2, sd=sb2);
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR %then
	%let plineb1 = prior b1arr: ~ normal(mb1, sd=sb1);

/* eMKF: Prior for variance parameters */
%let plinev=;
%if %upcase(&brndvars) = YES %then
	%let plinev = prior varr: ~ igamma(&bvshape, scale=&bvscale);

/******************************************************************************/
/* eMKF: Symbolic initialization for model parameters (resolved in proc mcmc) */
/******************************************************************************/

/* eMKF: Initial values for AR parameters */
%let initlinetau = ; %let initlinepsi = ; %let _i = 0;
%if %upcase(&bARmodel) = COMMON_AR %then %do; 	/*common AR parameters */
	%let initlinepsi = psi = &bmrho + sqrt(1/&bprho)*rand('normal');
	%let initlinetau = tau = rand('uniform', &btaul, &btauu); 
%end;
%if %upcase(&bARmodel) = INDEP_AR %then %do; /* Group-specific AR parameters */
    %let initlinepsi = &initlinepsi spsi = rand('uniform', 0.0001, sqrt(1/&bprho))%str(;) ;		
    %let initlinepsi = &initlinepsi mpsi = &bmrho + sqrt(1/&bprho)*rand('normal')%str(;) ;
	%do _i = 1 %to &g;
		%let initlinepsi = &initlinepsi psi&_i=mpsi+spsi*rand('normal')%str(;) ;
		%let initlinetau = &initlinetau tau&_i=rand('uniform',&btaul,&btauu)%str(;) ;
	%end;	
%end;

/* eMKF: Initial values for regression hyper-parameters in the fully Bayesian models */
%let initlinemb1=; %let initlinemb2=; %let initlinemb3=; %let initlinesb1=; %let initlinesb2=; %let initlinesb3=;
%if %upcase(&btype) = FULL_CUBIC %then %do;
	%let initlinesb3 = sb3 = rand('uniform', &bbeta3l, &bbeta3u);
	%let initlinemb3 = mb3 = &bmbeta3 + sqrt(1/&bpbeta3)*rand('normal') ;
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD %then %do;
	%let initlinesb2 = sb2 = rand('uniform', &bbeta2l, &bbeta2u);
	%let initlinemb2 = mb2 = &bmbeta2 + sqrt(1/&bpbeta2)*rand('normal') ;
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR %then %do;
	%let initlinesb1 = sb1 = rand('uniform', &bbeta1l, &bbeta1u);
	%let initlinemb1 = mb1 = &bmbeta1 + sqrt(1/&bpbeta1)*rand('normal') ;
%end;

/* eMKF: Initial values for prior mean vector mbetag and precision matrix Dbetag for use with matrix operations */
%let initmbeta = call zeromatrix(Dbetag);	
%let initmbeta = &initmbeta%str(;) mbetag[1,1] = &bmalpha%str(;) Dbetag[1,1] = &bpalpha;

/* eMKF: In models other than the fully Bayesian trend models, mbetag and Dbetag are constants */
%if %upcase(&btype) ^= FULL_CUBIC and %upcase(&btype) ^= FULL_QUAD and %upcase(&btype) ^= FULL_LINEAR and %upcase(&btype) ^= DROPPED
	%then %let initmbeta = &initmbeta%str(;) mbetag[2,1] = &bmbeta1%str(;) Dbetag[2,2] = &bpbeta1;
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = INDEP_QUAD or %upcase(&btype) = COMMON_CUBIC or %upcase(&btype) = COMMON_QUAD 
	%then %let initmbeta = &initmbeta%str(;) mbetag[3,1] = &bmbeta2%str(;) Dbetag[3,3] = &bpbeta2;
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = COMMON_CUBIC 
	%then %let initmbeta = &initmbeta%str(;) mbetag[4,1] = &bmbeta3%str(;) Dbetag[4,4] = &bpbeta3;

/* eMKF: In fully Bayesian models, mbetag and Dbetag depend on model parameters and are updated accordingly */
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR 
	%then %let initmbeta = &initmbeta%str(;) mbetag[2,1] = mb1%str(;) Dbetag[2,2] = 1/sb1**2;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD
	%then %let initmbeta = &initmbeta%str(;) mbetag[3,1] = mb2%str(;) Dbetag[3,3] = 1/sb2**2;
%if %upcase(&btype) = FULL_CUBIC
	%then %let initmbeta = &initmbeta%str(;) mbetag[4,1] = mb3%str(;) Dbetag[4,4] = 1/sb3**2;

/* eMKF: Initial values for intercepts */
%let initlinea=; %let _i = 0;
%do _i = 1 %to &g; 
	%let initlinea = &initlinea a&_i = &bmalpha+sqrt(1/&bpalpha)*rand('normal')%str(;) ;
%end;

/* eMKF: Initial values for regression parameters */
%let initlineb1=; %let initlineb2=; %let initlineb3=; %let _i=0;
%if %upcase(&btype) = FULL_CUBIC %then %do;
	%do _i = 1 %to &g; 
		%let initlineb3 = &initlineb3 b3arr&_i=mb3+sb3*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD %then %do;
	%do _i = 1 %to &g; 
		%let initlineb2 = &initlineb2 b2arr&_i=mb2+sb2*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR %then %do;
	%do _i = 1 %to &g; 
		%let initlineb1 = &initlineb1 b1arr&_i=mb1+sb1*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = INDEP_CUBIC %then %do;
	%do _i = 1 %to &g; 
		%let initlineb3 = &initlineb3 b3arr&_i=&bmbeta3+sqrt(1/&bpbeta3)*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = INDEP_QUAD %then %do;
	%do _i = 1 %to &g; 
		%let initlineb2 = &initlineb2 b2arr&_i=&bmbeta2+sqrt(1/&bpbeta2)*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = INDEP_QUAD or %upcase(&btype) = INDEP_LINEAR %then %do;
	%do _i = 1 %to &g; 
		%let initlineb1 = &initlineb1 b1arr&_i=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
	%end;
%end;
%if %upcase(&btype) = COMMON_CUBIC %then
	%let initlineb3 = b3 = &bmbeta3 + sqrt(1/&bpbeta3)*rand('normal') ;
%if %upcase(&btype) = COMMON_CUBIC or %upcase(&btype) = COMMON_QUAD %then
	%let initlineb2 = b2 = &bmbeta2 + sqrt(1/&bpbeta2)*rand('normal') ;
%if %upcase(&btype) = COMMON_CUBIC or %upcase(&btype) = COMMON_QUAD or %upcase(&btype) = COMMON_LINEAR %then
	%let initlineb1 = b1 = &bmbeta1 + sqrt(1/&bpbeta1)*rand('normal') ;

/* eMKF: Initial values for unobserved true states predictions given regression parameters */
%let _i = 0; %let _j = 0; 
%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = INDEP_CUBIC %then %do;
    %do _i = 1 %to &g; 
	    %local initetamnarr&_i;   /*eMKF: broken up into one macro variable per group instead of single combined macro variable to avoid max length error (65534) */
  	    %do _j = 1 %to &n; 
	        %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1arr&_i+X[&_j,3]*b2arr&_i+X[&_j,4]*b3arr&_i%str(;) ;
	    %end;
    %end;
%end;
%if %upcase(&btype) = FULL_QUAD or %upcase(&btype) = INDEP_QUAD %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	  		%let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,1]*b1arr&_i+X[&_j,2]*b2arr&_i%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = FULL_LINEAR or %upcase(&btype) = INDEP_LINEAR %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
			%let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1arr&_i%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = COMMON_CUBIC %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	        %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1+X[&_j,3]*b2+X[&_j,4]*b3%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = COMMON_QUAD %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	        %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1+X[&_j,3]*b2%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = COMMON_LINEAR %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	        %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = DROPPED %then %do;
	%do _i = 1 %to &g; 
  	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	        %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i%str(;) ;
		%end;
	%end;
%end;

/* eMKF: Initial values for variance parameters from igamma(&bvshape, scale=&bvscale) (if applicable) */
%let initlinevarr = ; %let _i = 0;
%if %upcase(&brndvars) = YES %then %do; 
	%do _i = 1 %to &g;
		%let initlinevarr = &initlinevarr varr&_i=1/rand('gamma',&bvshape,1/&bvscale)%str(;) ;
	%end;
%end;

/* eMKF: temporary dataset for building group-specific design matrix */
data _bbdata1_;
  set _bbdata_(where=(_group_ = 1));
  keep &brtm.0 &brtm.1 &brtm.2 &brtm.3;
run;

/* eMKF: dimensionality */
%let p = 0;
%if %upcase(&btype) = DROPPED %then %let p = 1;
%if %upcase(&btype) = FULL_LINEAR or %upcase(&btype) = INDEP_LINEAR or %upcase(&btype) = COMMON_LINEAR %then %let p = 2;
%if %upcase(&btype) = FULL_QUAD   or %upcase(&btype) = INDEP_QUAD   or %upcase(&btype) = COMMON_QUAD   %then %let p = 3;
%if %upcase(&btype) = FULL_CUBIC  or %upcase(&btype) = INDEP_CUBIC  or %upcase(&btype) = COMMON_CUBIC  %then %let p = 4;
%let p = %eval(0+&p);

/* eMKF: applicable read_array statement for the design matrix */
%let rcXline = ;
%if &p = 1 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'));
%if &p = 2 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'));
%if &p = 3 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'), resolve('&brtm.2'));
%if &p = 4 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'), resolve('&brtm.2'), resolve('&brtm.3'));

/* eMKF: applicable read_array statement for the effective sample sizes */
%let rcNline = ;
%if %upcase(&brndvars) = YES %then %let	rcNline = rcN = read_array('_bbdata_', Narr, '_n');

/*************************************************************************************/
/* eMKF: Applicable UDS statements - see macros gibbs_uds_compile_** for definitions */
/*       These will be applied in the order defined here, and after any M-H samplers */
/*************************************************************************************/
%let udsline = ; 

/* eMKF: UDS statement for mean hyper-parameters in fully Bayesian models (if applicable) */
/*       The pseudo-parameters mbetag and Dbetag are also updated in those UDS subroutines */
/*       Note that the sb1-sb3 are updated via the M-H sampler, which is applied before the UDS per SAS documentation */
%if %upcase(&btype) = FULL_CUBIC %then
	%let udsline = &udsline uds MP_bfc(mb1, mb2, mb3, mbetag, Dbetag, b1arr, b2arr, b3arr, mb1hyp, mb2hyp, mb3hyp, sb1, sb2, sb3)%str(;) ;
%if %upcase(&btype) = FULL_QUAD %then
	%let udsline = &udsline uds MP_bfq(mb1, mb2, mbetag, Dbetag, b1arr, b2arr, mb1hyp, mb2hyp, sb1, sb2)%str(;) ;
%if %upcase(&btype) = FULL_LINEAR %then
	%let udsline = &udsline uds MP_bfl(mb1, mbetag, Dbetag, b1arr, mb1hyp, sb1)%str(;) ;

/* eMKF: UDS statement for regression coefficients  */
/*       The pseudo-parameter etamnarr is also updated in those subroutines to hold the updated regression predictions */
%if %upcase(&btype) = DROPPED %then
	%let udsline = &udsline uds CP_b0(a, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = COMMON_LINEAR %then
	%let udsline = &udsline uds CP_b1l(a, b1, etamnarr, ambetag, bmbetag, aDbetag, bDbetag, rhoarr, nuarr, rts, aX, bX, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = COMMON_QUAD %then 
	%let udsline = &udsline uds CP_b1q(a, b1, b2, etamnarr, ambetag, bmbetag, aDbetag, bDbetag, rhoarr, nuarr, rts, aX, bX, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = COMMON_CUBIC %then
	%let udsline = &udsline uds CP_b1c(a, b1, b2, b3, etamnarr, ambetag, bmbetag, aDbetag, bDbetag, rhoarr, nuarr, rts, aX, bX, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = INDEP_LINEAR or %upcase(&btype) = FULL_LINEAR %then
	%let udsline = &udsline uds CP_bgl(a, b1arr, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = INDEP_QUAD or %upcase(&btype) = FULL_QUAD %then
	%let udsline = &udsline uds CP_bgq(a, b1arr, b2arr, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = INDEP_CUBIC or %upcase(&btype) = FULL_CUBIC %then
	%let udsline = &udsline uds CP_bgc(a, b1arr, b2arr, b3arr, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;

/* eMKF: UDS statement for true states etaarr */
%let udsline = &udsline uds EP(etaarr, etamnarr, rhoarr, nuarr, rts, Yarr, Sarr)%str(;) ;

/* eMKF: UDS statement for variances (if applicable) */
%if %upcase(&brndvars) = YES %then 
	%let udsline = &udsline uds RP(varr, vhyp, Sarr, Narr)%str(;) ;

/* eMKF: library location for pre-compiled UDS subroutines */
options cmplib = &bcmploc;

/* eMKF: Options will be the proc mcmc defaults if not specified by the user */
%let optionline=;
%if &bseed ^= %str()  %then %let optionline = &optionline seed 		= %eval(0+&bseed);;
%if &bmaxt ^= %str()  %then %let optionline = &optionline maxtune 	= %eval(0+&bmaxt);;
%if &btune ^= %str()  %then %let optionline = &optionline ntu 		= %eval(0+&btune);;
%if &bburn ^= %str()  %then %let optionline = &optionline nbi 		= %eval(0+&bburn);;
%if &biter ^= %str()  %then %let optionline = &optionline nmc 		= %eval(0+&biter);;
%if &bthin ^= %str()  %then %let optionline = &optionline thin 		= %eval(0+&bthin);;
%if &batol ^= %str()  %then %let optionline = &optionline accepttol = %sysevalf(&batol);;
%if &bttol ^= %str()  %then %let optionline = &optionline targaccept = %sysevalf(&bttol);;
%if &bprcov ^= %str() %then %let optionline = &optionline propcov 	= &bprcov;
%if &binit ^= %str()  %then %let optionline = &optionline init 		= &binit;

/* eMKF: Disable summary statistics if not requested by the user */
%if %upcase(&bprint) ^= YES %then %let optionline = &optionline stats = none;

/* eMKF: Diagnostics plots and ODS graphics enabled if requested by the user */
%if %upcase(&bplot) = YES %then %do; 
	%let optionline = &optionline plots = all;
	ods graphics on;
%end;
%else %let optionline = &optionline plots = none;

/* eMKF: Add jointmodel option (log-likelihood constructed using stored arrays) */
%let optionline = &optionline jointmodel;

/* eMKF: Monitor selected model parameters */
%let monitorline = &sbparline &mbparline &parline etaarr &vparline;
%if %upcase(&bARmodel) = INDEP_AR %then %let monitorline = spsi mpsi &tausqparline &rhoparline &monitorline ;
%if %upcase(&bARmodel) = COMMON_AR %then %let monitorline = tausq rho &monitorline ;

/* eMKF: Empty dataset to pass to proc mcmc: data from _bbdata_ will be read directly into arrays */
data _bb_;
run;

/* eMKF: Call proc mcmc using the above customizations  */
%put Call to PROC MCMC initiated; %let _i = 0;
proc mcmc data=_bb_ outpost= &blog monitor = ( &monitorline ) &optionline;;	

	  %if %upcase(&bprint) ^=YES and %upcase(&bplot) ^=YES 	/* Disable output tables and plots as applicable */
		%then ods select none;;

	  /**********************/
	  /* Array declarations */
	  /**********************/
	  array rts[&n] (&_brtimess); 	 						/* constant array with real times */
	  array Xarr[1]						   	    /nosymbols;	/* dynamic array for predictors to read in from dataset */
	  array Yarr[1]			   	  			    /nosymbols;	/* dynamic array for _y from dataset */
	  array Sarr[1]			   	   			    /nosymbols;	/* dynamic array for _var from dataset */
	  &Narrline;;											/* dynamic array for _n from dataset (if applicable) */
	  array X[&n, &p];										/* design matrix to use in matrix multiplication */
	  array mbetag[&p, 1];									/* prior mean vector for betas (assumed common accross groups) */
	  array Dbetag[&p, &p];									/* diagonal prior precision matrix for betas (assumed common accross groups) */
	  %if %upcase(&btype) = COMMON_LINEAR or 
			%upcase(&btype) = COMMON_QUAD or 
				%upcase(&btype) = COMMON_CUBIC %then %do;	/* conformal arrays used in UDS subroutines for common trend models */
		 array aX[&n, 1];									/* design matrix (intercept only) */
		 array bX[&n, %eval(&p-1)];							/* design matrix (excl. intercept) */
	  	 array ambetag[1, 1];								/* hyper parameter vector (intercept only) */
		 array bmbetag[%eval(&p-1), 1];						/* hyper parameter vector (excl. intercept) */
		 array aDbetag[1, 1];								/* hyper parameter matrix (intercept only) */
		 array bDbetag[%eval(&p-1),%eval(&p-1)];			/* hyper parameter matrix (excl. intercept) */
	  %end;
	  %if %upcase(&bARmodel) = INDEP_AR %then %do;			/* AR-related parameters in the group-specific random effects model */
	 	  array psi[&g] psi1-psi&g;							/* group-specific psi = ln[(1-rho)/(1+rho)] */
	  	  array tau[&g] tau1-tau&g;				    		/* group-specific innovation SD tau */
	  	  array tausq[&g] tausq1-tausq&g;					/* squares of group-specific innovation SD tau */
	  	  array rho[&g] rho1-rho&g;							/* reverse-transformation for rho */
	  	  array nu[&g] nu1-nu&g;							/* innovation variance parameters under stationarity */
	  	  array dg[&g] dg1-dg&g;							/* determinants of AR variance-covariance matrices */
	  %end;
	  array rhoarr[&g]; 									/* temporary 1-dimensional array with group-specific parameters rho */
	  array nuarr[&g]; 										/* temporary 1-dimensional array with group-specific parameters nu */
	  array a[&g] a1-a&g;									/* named 1-dimensional array of group-specific intercepts */
	  &b1line;;												/* named 1-dimensional array of group-specific linear coefficients (if requested) */
 	  &b2line;;                     						/* named 1-dimensional array of group-specific quad coefficients (if requested) */
	  &b3line;; 											/* named 1-dimensional array of group-specific cubic coefficients (if requested) */
	  &etamnarrline;;										/* named 1-dimensional array etamnarr (gxn) for predictions from regression */
	  &vline;;												/* named 1-dimensional array group-specific variance parameters (if requested) */
	  &etaarrline;;											/* named 1-dimensional array etaarr (gxn) for unobserved true states */

	  begincnst;

		  /*****************/
	  	  /* Design matrix */
		  /*****************/
	  	  &rcXline;;										/* read in dynamic array of predictors Xarr */  
		  call zeromatrix(X);								/* initialize design matrix X to all zeroes */
		  if &p > 1 then do;								/* Xarr is a 2-dimensional array for p > 1 */
		 	  do i = 1 to &n;
			 	  do m = 1 to &p;
					  X[i,m] = Xarr[i,m];						
			  	  end;
		  	  end;
		  end;
		  else do;											/* !! Xarr collapses to a 1-dimensional array if p = 1 !! */
		 	  do i = 1 to &n;
			  	  X[i,1] = Xarr[i];						
		  	  end;
		  end;
		  %if %upcase(&btype) = COMMON_LINEAR or 
			    %upcase(&btype) = COMMON_QUAD or 
				  %upcase(&btype) = COMMON_CUBIC %then %do; /* conformal subarrays used in UDS subroutines for common trend models */
			  do i = 1 to &n;								
				  aX[i, 1] = X[i, 1];
				  do m = 2 to &p;
					  bX[i, m-1] = X[i, m];
				  end;
			  end;
		  %end;

		  /**********************/
		  /* Group sample means */
		  /**********************/
		  rcY = read_array('_bbdata_', Yarr, '_y');			/* read in 1-dimensional array of _y from dataset */

		  /**********************/
	  	  /* Sampling variances */
		  /**********************/
		  rcS = read_array('_bbdata_', Sarr, '_var');		/* read in 1-dimensional array of _var from dataset */

		  /******************************************/
		  /* Effective sample sizes (if applicable) */
		  /******************************************/
		  &rcNline;;										/* read in 1-dimensional array of _n from dataset (if applicable) */

		  /******************/
	 	  /* Initialization */
	  	  /******************/
		  call streaminit(%eval(0+&bseed));					/* set seed */

		  &initlinepsi;;									/* initialize psi = ln[(1-rho)/(1+rho)] */
		  &initlinetau;;									/* initialize innovation SD tau */
		  %if %upcase(&bARmodel) = COMMON_AR %then %do;		/* common AR parameters across groups */
			  tausq = tau**2;					 		 	/* track tau-squared */
	  		  rho = (1-exp(psi))/(1+exp(psi)); 		 		/* reverse-transformation for rho */
	  		  nu = tausq/(1-rho**2);			 		 	/* innovation variance parameter under stationarity */
			  dg = nu**&n;									/* recursive formula for determinant of Vgamma (assuming 2+ points) */
			  do i = 2 to &n;								
				  dg = dg*(1-(rho**(2*(rts[i]-rts[i-1])))); 
			  end;
			  if abs(rho) ge 1 or dg= . or dg le 0 then do; /* guard against numerical singularities */
				  rho = 0;
				  nu = tausq;
				  dg = nu**&n;
			  end;
			  do k=1 to &g;									/* temp parameter arrays (e.g., to pass to UDS subroutines) */
				  rhoarr[k] = rho;
			  	  nuarr[k] = nu;
			  end;
		  %end;
		  %if %upcase(&bARmodel) = INDEP_AR %then %do;		/* independent AR parameters across groups */
			  do k=1 to &g;
				  tausq[k] = tau[k]**2;		 	
		  		  rho[k] = (1-exp(psi[k]))/(1+exp(psi[k]));
		  		  nu[k] = tausq[k]/(1-rho[k]**2);
				  dg[k] = nu[k]**&n;
				  do i = 2 to &n;
					  dg[k] = dg[k]*(1-(rho[k]**(2*(rts[i] - rts[i-1]))));
				  end;
				  if abs(rho[k]) ge 1 or dg[k] = . or dg[k] le 0 then do;
					  rho[k] = 0;             
					  nu[k] = tausq[k]; 		 
				      dg[k] = nu[k]**&n;
				  end;
				  rhoarr[k] = rho[k];
			  	  nuarr[k] = nu[k];
			  end;
		  %end;

		  &initlinesb1;;									/* initialize SD hyperparameter sb1 (if applicable) */
	  	  &initlinesb2;;									/* initialize SD hyperparameter sb2 (if applicable) */
	  	  &initlinesb3;;									/* initialize SD hyperparameter sb3 (if applicable) */
		  &initlinemb1;;									/* initialize mean hyperparameter mb1 (if applicable) */
	  	  &initlinemb2;;									/* initialize mean hyperparameter mb2 (if applicable) */
	  	  &initlinemb3;;									/* initialize mean hyperparameter mb3 (if applicable) */
		  &initmbeta;;							 			/* initialize mbetag and Dbetag */
		  %if %upcase(&btype) = COMMON_LINEAR or 
				%upcase(&btype) = COMMON_QUAD or 
				  %upcase(&btype) = COMMON_CUBIC %then %do;	/* conformal subarrays used in UDS subroutines for common trend models */
			  call zeromatrix(aDbetag);
			  call zeromatrix(bDbetag); 
		  	  ambetag[1,1] = mbetag[1,1];
		  	  aDbetag[1,1] = Dbetag[1,1];	
			  do m = 2 to &p;
				  bmbetag[m-1, 1] = mbetag[m,1];
			      bDbetag[m-1, m-1] = Dbetag[m,m];	
			  end;
		  %end;

		  &initlinea;;										/* initialize intercepts */
		  &initlineb1;;										/* initialize linear coefficients (if applicable) */
		  &initlineb2;;										/* initialize quad coefficients (if applicable) */	
		  &initlineb3;;										/* initialize cubic coefficients (if applicable) */

          %do _i = 1 %to &g; 
		      &&initetamnarr&_i;;						    /* initialize conditional mean for true states  */ 
		  %end;

		  do k = 1 to &g; 			  						/* initialize etaarr using Markov property of AR process */
			  etaarr[(k-1)*&n+1] = etamnarr[(k-1)*&n+1] + 
							sqrt(nuarr[k])*rand('normal'); 	/* first timepoint from stationary distribution of AR process */

		  	  do i = 2 to &n; 								/* subsequent timepoints from implied conditional distributions */
			      etaarr[(k-1)*&n+i] = etamnarr[(k-1)*&n+i] + 
							((rhoarr[k]**(rts[i] - rts[i-1]))*(etaarr[(k-1)*&n+i-1] - etamnarr[(k-1)*&n+i-1])) + 
							sqrt(nuarr[k]*(1-(rhoarr[k]**(2*(rts[i] - rts[i-1])))))*rand('normal');
		  	  end;
		  end;

		  &initlinevarr;;									/* initialize sampling variances (if applicable) */

	  endcnst;

	  /********************/
	  /* UDS declarations */
	  /********************/
	  &udsline;;											/* Gibbs sampling done in the order specified in udsline */
	  														/* Per SAS documentation, parameters that use M-H will be sampled first */
	  /**************************/
	  /* Parameter declarations */
	  /**************************/
	  &psiparline2;;										/* psi = ln[(1-rho)/(1+rho)] and any hyperparameters */
	  &tauparline2;;										/* innovation SD tau */
	  &sbparline2;;											/* SD hyper-parameters (if any) */
	  &udsparline;;											/* UDS parameter blocks, one for each Gibbs sampler */

	  beginnodata;

	  	  /********************/
	  	  /* Prior statements */
	  	  /********************/
		  &hplinespsi;;										/* SD hyper-prior for mean of psi = ln[(1-rho)/(1+rho)] (if applicable) */
		  &hplinempsi;;										/* Mean hyper-prior for mean of psi = ln[(1-rho)/(1+rho)] (if applicable) */
		  &plinepsi;;										/* prior for psi = ln[(1-rho)/(1+rho)] */
		  &plinetau;;										/* prior for innovation SD tau */
		  %if %upcase(&bARmodel) = COMMON_AR %then %do;		/* AR parameters in the common case */
			  tausq = tau**2;					 			/* track tau-squared */
	  		  rho = (1-exp(psi))/(1+exp(psi)); 		 		/* reverse-transformation for rho */
	  		  nu = tausq/(1-rho**2);			 			/* innovation variance parameter under stationarity */
			  dg = nu**&n;									/* recursive formula for determinant of Vgamma (assuming 2+ points) */
			  do i = 2 to &n;								
				  dg = dg*(1-(rho**(2*(rts[i]-rts[i-1]))));
			  end;
			  if abs(rho) ge 1 or dg= . or dg le 0 then do; /* guard against numerical singularities */
				  rho = 0;
				  nu = tausq;
				  dg = nu**&n;
			  end;
			  do k=1 to &g;									/* parameter arrays to pass to UDS subroutines */
				  rhoarr[k] = rho;
			  	  nuarr[k] = nu;
			  end;
		  %end;
	      %if %upcase(&bARmodel) = INDEP_AR %then %do;	 	/* Group-specific AR parameters */
			 do k = 1 to &g; 
			    tausq[k] = tau[k]**2;
		  	    rho[k] = (1-exp(psi[k]))/(1+exp(psi[k]));
		  	    nu[k] = tausq[k]/(1-rho[k]**2);
			    dg[k] = nu[k]**&n;
			    do i = 2 to &n;
				    dg[k] = dg[k]*(1-(rho[k]**(2*(rts[i] - rts[i-1]))));
			    end;
			    if abs(rho[k]) ge 1 or dg[k] = . or dg[k] le 0 then do;
				    rho[k] = 0;             
				    nu[k] = tausq[k]; 		 
				    dg[k] = nu[k]**&n;
			    end;
				rhoarr[k] = rho[k];
			  	nuarr[k] = nu[k];		
			 end;
		  %end;

	  	  &hplinesb1;;						 	 			/* SD hyper-priors in the full_linear, full_quad, and full_cubic cases */
	  	  &hplinesb2;;						 	 			/* SD hyper-priors in the full_quad and full_cubic cases */
	  	  &hplinesb3;;						 	 			/* SD hyper-prior in the full_cubic case */
	      &hplinemb1;;						 	 			/* Mean hyper-priors in the full_linear, full_quad, and full_cubic cases */
	  	  &hplinemb2;;						 	 			/* Mean hyper-priors in the full_quad and full_cubic cases */	
      	  &hplinemb3;; 						 	 			/* Mean hyper-prior in the full_cubic case */
		  %if %upcase(&btype) = FULL_LINEAR or 
				%upcase(&btype) = FULL_QUAD or 
				  %upcase(&btype) = FULL_CUBIC %then	
		  	  prior mbetag Dbetag ~ general(0);;			/* pseudo-parameters mbetag and Dbetag do no contribute to prior */
	  	  &plinea;; 							 			/* prior for intercepts */
	  	  &plineb1;;							 			/* prior for linear coefficients (if those were requested) */
	  	  &plineb2;;							 			/* prior for quadratic coefficients (if those were requested) */
	  	  &plineb3;;							 			/* prior for cubic coefficients (if those were requested) */

		  prior etamnarr ~ general(0);						/* pseudo-parameters etamnarr do no contribute to prior */

		  lpr = 0; 											/* calculation of log-prior for etaarr from univariate conditionals */
		  do k = 1 to &g;
		  	  lpr = lpr + lpdfnorm(etaarr[(k-1)*&n+1], 		/* etaarr is updated in the UDS call for the true states */
								   etamnarr[(k-1)*&n+1], 	/* etamnarr is updated in the UDS call for the regression coefficients */
								   sqrt(nuarr[k]));			/* first timepoint from stationary distribution of AR process */
		  	  do i = 2 to &n; 								/* subsequent timepoints from implied conditional distributions */
			      lpr = lpr + lpdfnorm(etaarr[(k-1)*&n+i], 
									   etamnarr[(k-1)*&n+i] + 
									    (rhoarr[k]**(rts[i] - rts[i-1]))*(etaarr[(k-1)*&n+i-1] - etamnarr[(k-1)*&n+i-1]), 
									   sqrt(nuarr[k]*(1-(rhoarr[k]**(2*(rts[i] - rts[i-1])))))); 
		  	  end;
		  end;
		  prior etaarr ~ general(lpr);						/* prior for unobserved true states */

		  &plinev;;								 			/* Inverse gamma prior for sampling variances (if applicable) */

		  /********************************/
	  	  /* Loglikelihood calculation(s) */
	  	  /********************************/
	  	  lp = 0;		  									/* log of joint distribution of sample means */
		  do k = 1 to &g*&n;
		  	  lp = lp + lpdfnorm(Yarr[k], etaarr[k], sqrt(Sarr[k]));
		  end;

		  %if %upcase(&brndvars) = YES %then %do;		  	/* log of joint distribution of sample variances (if applicable) */
			  do k = 1 to &g;
				  do j = 1 to &n;
					  lp = lp + lpdfgamma(Sarr[(k-1)*&n+j],
										  (Narr[(k-1)*&n+j]-1)/2,
										  (2*varr[k])/(Narr[(k-1)*&n+j]-1)); /* varr is updated in the UDS call for the variances */
			 	  end;
			  end;
		  %end;

	  endnodata;

	  /*******************/
	  /* Model statement */
	  /*******************/
	  model general(lp);

run;

/* eMKF: Re-enable output tables and plots */
%if %upcase(&bprint) ^= YES and %upcase(&bplot) ^= YES %then ods select all;;

/* eMKF: Disable ODS graphics */
%if %upcase(&bplot) = YES %then ods graphics off;;

%put Call to PROC MCMC concluded;

/*eMKF: Keep only the desired columns in the posterior log dataset */
%if %upcase(&bARmodel) = INDEP_AR %then %do;
	data &blog;
	  merge &blog(drop= etamn: Log: spsi mpsi &tauparline &psiparline
				    	%if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or 
					        %upcase(&btype) = FULL_LINEAR %then mbetag: Dbetag: ; ) 
			&blog(keep = spsi mpsi);
	run;
%end;
%if %upcase(&bARmodel) = COMMON_AR %then %do;
	data &blog;
		set &blog(drop= etamn: Log: tau psi
				       %if %upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or 
					    	%upcase(&btype) = FULL_LINEAR %then mbetag: Dbetag: ; )
		;
	run;
%end;

/****************************************************/
/* eMKF: Reverse-transform regression coefficients  */
/****************************************************/

data _blogc2_ _tblogc2_ _tblogc_  ;
run;

%let _i = 0; 

%if %upcase(&borpoly) = YES %then %do;

	/* eMKF: order columns by group */
	data _blogc2_;
  	  retain  Iteration 
			  %do _i=1 %to &g;
				  a&_i
				  %if %upcase(&btype) = FULL_LINEAR or %upcase(&btype) = INDEP_LINEAR %then b1arr&_i;
				  %if %upcase(&btype) = FULL_QUAD   or %upcase(&btype) = INDEP_QUAD   %then b1arr&_i b2arr&_i;
				  %if %upcase(&btype) = FULL_CUBIC  or %upcase(&btype) = INDEP_CUBIC  %then b1arr&_i b2arr&_i b3arr&_i;
			  %end;
			  %if %upcase(&btype) = COMMON_LINEAR %then b1 ; 
			  %if %upcase(&btype) = COMMON_QUAD   %then b1 b2 ; 
			  %if %upcase(&btype) = COMMON_CUBIC  %then b1 b2 b3 ; 
	  ;
	  set &blog(keep = Iteration a: %if &p > 1 %then b: ;);
	run;

	/* eMKF: block diagonal by group */
	%let oPPmat = ; %let _i = 0;
	%do _i = 1 %to &g; 
		%if &_i = 1 %then %let oPPmat = block( oP ;
		%if &_i > 1 and &_i < &g %then %let oPPmat = &oPPmat , block ( oP ;
		%if &_i = &g and &g > 1  %then %let oPPmat = &oPPmat , oP %sysfunc(repeat( %str(%)), &g-2));
		%if &_i = &g and &g = 1  %then %let oPPmat = &oPPmat );
	%end;

	%let _i = 0;

	/* eMKF: call proc iml to perform matrix multiplication */
	proc iml;

		use _oPmat_;
		read all into oP; close _oPmat_;
		oP = oP[1:&p, 1:&p];
		oPP = &oPPmat;;

		/* eMKF: re-structure block matrix in the common trend cases (where &p > 1) */
		%if %upcase(&btype) = COMMON_LINEAR or %upcase(&btype) = COMMON_QUAD or %upcase(&btype) = COMMON_CUBIC %then %do;
			oPP1 = oPP[do(1, &g*&p, &p), do(1, &g*&p, &p)];
			oPP2 = vecdiag(oPP[do(1, &g*&p, &p), do(2, &g*&p, &p)]);
			ToPP2 = vecdiag(oPP[do(2, &g*&p, &p), do(1, &g*&p, &p)]);
			oPP1 = oPP1 // T(ToPP2);
			oPP0 = oPP2;
			%if &p > 2 %then %do;
				oPP3 = vecdiag(oPP[do(1, &g*&p, &p), do(3, &g*&p, &p)]);
				ToPP3 = vecdiag(oPP[do(3, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP3);
				oPP0 = oPP0 || oPP3;
			%end;
			%if &p > 3 %then %do;
				oPP4 = vecdiag(oPP[do(1, &g*&p, &p), do(4, &g*&p, &p)]);
				ToPP4 = vecdiag(oPP[do(4, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP4);
				oPP0 = oPP0 || oPP4;
			%end;
			oPP0 = oPP0 // oPP[2:&p, 2:&p];
			oPP = oPP1 || oPP0;
		%end;

		varNames = {"Iteration"};
		%do _i = 1 %to &g;
		  varNames = varNames || {"a&_i"};
		  %if %upcase(&btype)= FULL_LINEAR or %upcase(&btype)= INDEP_LINEAR %then varNames = varNames || {"b1arr&_i"};;
		  %if %upcase(&btype)= FULL_QUAD   or %upcase(&btype)= INDEP_QUAD   %then varNames = varNames || {"b1arr&_i"} || {"b2arr&_i"};;
		  %if %upcase(&btype)= FULL_CUBIC  or %upcase(&btype)= INDEP_CUBIC  %then varNames = varNames || {"b1arr&_i"} || {"b2arr&_i"} || {"b3arr&_i"};;
		%end;
		%if %upcase(&btype) = COMMON_LINEAR %then varNames = varNames || {"b1"};;
		%if %upcase(&btype) = COMMON_QUAD   %then varNames = varNames || {"b1" "b2"};;
		%if %upcase(&btype) = COMMON_CUBIC  %then varNames = varNames || {"b1" "b2" "b3"};;

		use _blogc2_;
		read all into oB;
		close _blogc2_;

		oB1 = oB[,1];
		oB = T(oB[,2:ncol(oB)]);
		oBB = oPP * oB;
		oBB = oB1 || T(oBB);

		create _tblogc2_ var varNames;
		append from oBB;
		close _tblogc2_;

	quit;

	/* eMKF: re-order columns as they were initially from PROC MCMC */
	data _tblogc_;
  	  retain  Iteration a1-a&g 
			  %if %upcase(&btype) = FULL_LINEAR or %upcase(&btype) = INDEP_LINEAR  %then b1arr1-b1arr&g ; 
			  %if %upcase(&btype) = FULL_QUAD   or %upcase(&btype) = INDEP_QUAD    %then b1arr1-b1arr&g b2arr1-b2arr&g ; 
			  %if %upcase(&btype) = FULL_CUBIC  or %upcase(&btype) = INDEP_CUBIC   %then b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g ; 
			  %if %upcase(&btype) = COMMON_LINEAR %then b1 ; 
			  %if %upcase(&btype) = COMMON_QUAD   %then b1 b2 ; 
			  %if %upcase(&btype) = COMMON_CUBIC  %then b1 b2 b3 ; 
	  ;
	  set _tblogc2_;
	run;

	/* eMKF: merge into &blog */
	data &blog;
	  merge &blog(keep = Iteration %if %upcase(&btype)=FULL_LINEAR or %upcase(&btype)=FULL_QUAD or %upcase(&btype)=FULL_CUBIC %then sb: mb: ;)
	        _tblogc_
			&blog(drop = %if %upcase(&btype)=FULL_LINEAR or %upcase(&btype)=FULL_QUAD or %upcase(&btype)=FULL_CUBIC %then sb: mb: ;
						 a: %if &p > 1 %then b: ; )
	  ;
	  by Iteration;
	run;

%end;

/* eMKF: clean-up */
proc datasets nolist;
 delete _bbdata_ _bbdata1_ _bb_ _bfreqg_ _bfreqn_ _bbjunk _oXmat_ _oPmat_ _tblogc_ _tblogc2_ _blogc2_;
run ;
quit;

%mend;

data _null_;
run;




/* eMKF: BAYESBMA -- Implements Bayesian model averaging via mixture prior approach -- 2023 Q2 by M. Talih
 bdata              : Name of the data to be used
 blog               : Name of the output data containing full set of &biter/&bthin posterior draws
 btype              : bma_cubic, bma_quad, or bma_linear
 bgroup             : Group variable in the dataset 
 btime              : Time variable in the dataset 
 boutcome           : Outcome of interest variable in the dataset 
 bse                : Standard error variable in the dataset 
 bn				    : Effective sample size variable in the dataset (if applicable)
 brndvars			: YES if variances should be modeled; NO if variances should be assumed known
 bARmodel			: common_ar if AR parameters are common across groups; indep_ar if they are independently drawn from a common prior
 bslicesampler		: YES to use the slice sampler instead of MH algorithm for parameters that are not included in the Gibbs sampling step 
					  Default is NO due to heavier computational load.
 bseed              : random number generating seed that will allow the user to reproduce the same results in the Bayesian model
 bprcov				: method used in constructing initial covariance matrix for the MH algorithm (see proc mcmc documentation)
					  If empty, proc mcmc default of IND will be used.
 binit				: Option for generating initial values for the parameters (see documentation and leave empty to apply proc mcmc default)
					  eMKF default is REINIT to reset chains after tuning at the values set by the user
 bmaxt				: maximum number of proposal tuning loops (if empty, proc mcmc default of 24 is used; if 0, tuning will be skipped)
 batol				: Tolerance for acceptance probabilities (if empty, proc mcmc default of 0.075 is used in bttol +|- batol)
 bttol				: Target acceptance rate for random walk Metropolis. If empty, proc mcmc defaults are used, as follows: 
					  0.45 for models with 1 parameter, 0.35 for 2-4 parameters, and 0.234 for models with 5+ parameters.
 btune				: number of tuning iterations to use in each MCMC proposal tuning phase (if empty, proc mcmc default of 500 is used)
 bburn              : number of burn-in MCMC iterations (if empty, proc mcmc default of 1000 is used)
 biter              : number of post-burn-in MCMC iterations (if empty, proc mcmc default of 1000 is used)
 bthin				: controls thinning rate (if empty, proc mcmc default of 1 is used)
 borpoly  			: YES (default) for pre-transforming the design matrix using SAS IML orpol function. NO for "raw" polynomials.
          			  If YES, regression coefficients will be reverse-transformed prior to macro end. 
					  However, prior values below are assumed to be for the coefficients of the orthogonal polynomial regression if borpoly=YES.
 bmalpha , bpalpha 	: prior mean and precision for alphas
 bmbeta1			: prior mean for mean linear coefficient across groups -- used for cubic, quadratic, or linear trends (full_, indep_, and common_)
 bpbeta1 			: prior precision for mean linear coefficient across groups -- used for SD hyperprior(s) in all 3 full_ models
 bmbeta2 			: prior mean for mean quadratic coefficient across groups -- used for cubic or quadratic trends (full_, indep_, and common_)
 bpbeta2			: prior precision for mean quadratic coefficient across groups -- used for SD hyperprior(s) in full_cubic or full_quad
 bmbeta3 			: prior mean for mean cubic coefficient across groups -- used for cubic trends (full_, indep_, and common_)
 bpbeta3			: prior precision for mean cubic coefficient across groups -- used for SD hyperprior(s) in full_cubic
 bbeta1l, bbeta1u	: bounds for U(a,b) prior for SD of linear coefficients across groups -- used for hyperprior(s) in all 3 full_ models
 bbeta2l, bbeta2u	: bounds for U(a,b) prior for SD of quadratic coefficients across groups -- used for hyperprior(s) in full_cubic or full_quad
 bbeta3l, bbeta3u	: bounds for U(a,b) prior for SD of cubic coefficients across groups -- used for hyperprior(s) in full_cubic
 bmrho, bprho		: prior mean and precision for transformed rho -- ie., psi = ln[(1-rho)/(1+rho)]
 btaul, btauu		: bounds for U(a,b) prior for tau (SD of innovation variance tausq)
 bvshape			: Shape parameter for inverse gamma prior distribution of the variance (when applicable) 
 bvscale			: Scale parameter for inverse gamma prior distribution of the variance (when applicable)
 vwshape			: Common shape parameter to use for Dirichlet prior on model indicators in mixture prior
 bprint				: If YES, posterior parameter estimates and default chain-specific convergence diagnostics are printed (default is NO)
 bplot				: If YES, trace/diagnostics plots from proc mcmc will be included (default is NO)
 bcmploc			: location of CMP library (usually set in parent macro mkf)

*/

%macro bayesBMA(
             bdata	= , 
			 blog	= ,
			 btype	= bma_linear, 
	   /* eMKF: Variable labels assumed to have been reformatted using macro reformat */
			 bgroup	= _group_, 
			 btime	= _time, 
			 boutcome= _y, 
			 bse	= _se,
			 bn 	= ,
			 brndvars = NO,
			 bARmodel = common_ar,
			 bslicesampler = NO,
	   /* eMKF: MCMC tuning parameters: if missing, proc mcmc defaults will be used */
			 bseed	= ,
			 bprcov = ,
			 binit  = reinit,
			 bmaxt  = ,
			 batol 	= ,	
			 bttol 	= ,
			 btune	= ,			
			 bburn  = ,
			 biter  = ,
			 bthin 	= ,
			 borpoly = YES,
	   /* eMKF: Model parameters: if missing, the data will be used to generate starting values*/
			 bmalpha  = , 		
			 bpalpha  = ,
			 bmbeta1  = 0,      /*eMKF: Constant c3 or c7 in RAND's MKF User's Guide */
			 bpbeta1  = ,
			 bmbeta2  = 0,
			 bpbeta2  = ,
			 bmbeta3  = 0,
			 bpbeta3  = ,
			 bbeta1l  = 0,		/*eMKF: Constant c5 in RAND's MKF User's Guide */
			 bbeta1u  = ,
			 bbeta2l  = 0,
			 bbeta2u  = ,
			 bbeta3l  = 0,
			 bbeta3u  = ,
             bmrho    = 0,		/*eMKF: Constant c9  in RAND's MKF User's Guide */
			 bprho    = 1,		/*eMKF: Constant c10 in RAND's MKF User's Guide */
			 btaul    = 0.0001,	/*eMKF: Constant c11 in RAND's MKF User's Guide */
			 btauu    = ,
			 bvshape  = ,
			 bvscale  = ,
			 bwshape  = 2,
	    /* eMKF: Printing and diagnostic plots are off by default */
			 bprint   = NO,
			 bplot 	  = NO,
			 bcmploc  = work.funcs
               );
 
%local g n p brtm _brtimess brangeY bqrangeV bmedianV
	   formatted dsop dscl _i _j _l _ll oPPmat
       b1line b2line b3line vline etaarrline etamnarrline tauparline psiparline tausqparline rhoparline wtsline
       parline parline2 aparline vparline mbparline sbparline udsparline tauparline2 psiparline2 sbparline2
       plinea plineb1 plineb2 plineb3 plinev plinetau plinepsi plinewts plineflg
       hplinemb1 hplinesb1 hplinemb2 hplinesb2 hplinemb3 hplinesb3 hplinempsi hplinespsi
	   initlinea initlineb1 initlineb2 initlineb3 initlinevarr initlinetau initlinepsi initlinewts initlineflg
	   initlinemb1 initlinemb2 initlinemb3 initlinesb1 initlinesb2 initlinesb3 
       monitorline optionline udsline rcXline rcNline initmbeta Narrline;

/* eMKF: Data assumed to have been pre-formatted using macro reformat: check and reformat if not */

%let formatted = 0;

%let dsop = %sysfunc(open(&bdata));
%if &dsop ne 0 %then %do;
	%if %sysfunc(varnum(&dsop, inputorder)) ne 0 and %sysfunc(varnum(&dsop, &btime)) ne 0 %then %let formatted = 1;
%end; 
%let dscl = %sysfunc(close(&dsop));

%let formatted = %eval(&formatted + 0);

data _bbdata_ _bbdata1_;
run;

%if &formatted = 1 %then %do;
	data _bbdata_;
	  set &bdata;
	run;
%end;
%else %do;
    %put ;
	%put Reformatting data prior to Bayesian estimation;
	%if %upcase(&brndvars) = YES and &bn = %str() %then %do;
		%put ERROR: (Effective) sample sizes bn must be specified to fit random sampling variances;
		proc iml;
			print "  Error Note:";
			print "  (Effective) sample sizes bn must be specified to fit random sampling variances. ";
		quit;
		%return;
	%end;
	%reformat(data=&bdata, 
		      outcome=&boutcome, se=&bse, neff=&bn, outcome2=, se2=, neff2=,
			  group=&bgroup, time=&btime, by=, randomVars = &brndvars,
 			  outformat= _bbdata_ );
%end;

/*eMKF: Sort by replications, group, and time */
proc sort data= _bbdata_;
  by _rep _group_ _time ;
run;

/* eMKF: Macro variable for the number of groups */
%let g=;
data _bfreqg_;
run;
proc freq data=_bbdata_ noprint;
 tables _group_ /list out=_bfreqg_;
run;
data _bfreqg_;
 set _bfreqg_;
 _grp_ +1;
 call symput('g',_grp_);
 keep _grp_ _group_;
run;

/* eMKF: Macro variable for the number of time points */
%let n=;
data _bfreqn_;
run;
proc freq data=_bbdata_ noprint;
 tables _rtime /list out=_bfreqn_;
run;
data _bfreqn_;
 set _bfreqn_;
 _tm +1;
 call symput('n',_tm);
 keep _tm _rtime;
run;

/*eMKF: Macro variable for the real times to use in calculations */
%let _brtimess = ;
data _bfreqn_;
  set _bfreqn_;
  retain _rts;
  if _n_= 1 then _rts = cat(_rtime);
  else _rts = catx(" ", _rts, _rtime);
  call symput('_brtimess', _rts);
  drop _rts;
run;

/* eMKF: variable that will be used for real time in case times are irregular */
%let brtm  = _rtime;

/* eMKF: Set numerical values to use in code */
%let n=%eval(0+&n);
%let g=%eval(0+&g);

/* eMKF: Compute variances */
data _bbdata_;
  set _bbdata_ ;
  _var = _se**2;
run;

/* eMKF: Modified to use orthogonal cubic polynomial design matrix */

data _oXmat_ _oPmat_;
run;

%if %upcase(&borpoly) = YES %then %do;

	proc iml;
	  x = { &_rtimess };	
	  x = T(x);							/* eMKF: column vector with real times */
	  oP = orpol(x, 3);					/* eMKF: orthonormal design matrix for cubic orthogonal polynomials */
	  x0 = { %cnstss(1, &n) };
	  x0 = T(x0);
	  x1 = x;
	  x2 = x#x;
	  x3 = x#x2;
	  uP = x0 || x1 || x2 || x3;		/* eMKF: raw/unstandardized design matrix */
	  oP1 = inv(T(uP)*uP)*T(uP)*oP[,1];
      oP2 = inv(T(uP)*uP)*T(uP)*oP[,2];
      oP3 = inv(T(uP)*uP)*T(uP)*oP[,3];
      oP4 = inv(T(uP)*uP)*T(uP)*oP[,4];
	  oPP = oP1 || oP2 || oP3 || oP4;	/* eMKF: right multiplication of raw uP with oPP produces orthonormal oP */
	  y = T(do(1, &n, 1));				/* eMKF: column vector of consecutive time indices */
	  yP = y || oP;
	  create _oXmat_ from yP [ colname = {"_time" "&brtm.0" "&brtm.1" "&brtm.2" "&brtm.3"} ] ;
	  append from yP; close _oXmat_;
	  create _oPmat_ from oPP [ colname = {"t0" "t1" "t2" "t3"} ] ;
	  append from oPP; close _oPmat_;
	quit;

	proc sort data=_bbdata_;
	  by _time;
	run;

	data _bbdata_;
	  merge _bbdata_ _oXmat_;
	  by _time;
	run;

	proc sort data= _bbdata_;
	  by _rep _group_ _time ;
	run;

%end;
%else %do; 								/* eMKF: Add raw predictor variables */
	data _bbdata_;
	  set _bbdata_;
	  &brtm.0 = 1;
	  &brtm.1 = &brtm;
	  &brtm.2 = &brtm**2;
	  &brtm.3 = &brtm**3;
	run;
%end;

/* eMKF: Evaluate range of the data to use in setting prior parameters, as in MKF */
%let brangeY=;
data _bbjunk;
run;
proc means data=_bbdata_ noprint;
  var _y;
  output out=_bbjunk range=range;
run;
data _null_;
 set _bbjunk;
 call symput("brangeY", range);
run;
%let brangeY = %sysevalf(&brangeY + 0);

/*******************************************************************/
/* eMKF: Set any prior parameters not already provided by the user */
/*******************************************************************/

/* eMKF: c1 in RAND's MKF User's Guide */
%if &bmalpha = %str() %then %let bmalpha = %sysevalf(0.5 * &brangeY);;

/* eMKF: 1/c2 in RAND's MKF User's Guide */	
%if &bpalpha = %str() %then %let bpalpha = %sysevalf(0.000001/(&brangeY**2));; 

/* eMKF: c3 or c7 in RAND's MKF User's Guide */
%if &bmbeta1 = %str() %then %let bmbeta1 = %sysevalf(0);;

/* eMKF: 1/c4 in RAND's MKF User's Guide */ 
%if &bpbeta1  = %str() and (%upcase(&btype) = FULL_CUBIC or %upcase(&btype) = FULL_QUAD or %upcase(&btype) = FULL_LINEAR) 
	%then %let bpbeta1  = %sysevalf(10/(&brangeY**2));;

/* eMKF: 1/c8 in RAND's MKF User's Guide */	
%if &bpbeta1  = %str() and %upcase(&btype) ^= FULL_CUBIC and %upcase(&btype) ^=FULL_QUAD and %upcase(&btype) ^=FULL_LINEAR 
	%then %let bpbeta1  = %sysevalf(0.000001/(&brangeY**2));;

/* eMKF: c5 in RAND's MKF User's Guide  */
%if &bbeta1l  = %str() %then %let bbeta1l = %sysevalf(0);;						

/* eMKF: c6 in RAND's MKF User's Guide  */	
%if &bbeta1u  = %str() %then %let bbeta1u = %sysevalf(0.5 * &brangeY);;		

/* eMKF: c9 in RAND's MKF User's Guide  */	
%if &bmrho    = %str() %then %let  bmrho  = %sysevalf(0);;	

/* eMKF: c10 in RAND's MKF User's Guide  */	
%if &bprho    = %str() %then %let  bprho  = %sysevalf(1);;	

/* eMKF: c11 in RAND's MKF User's Guide  */	
%if &btaul    = %str() %then %let  btaul  = %sysevalf(0.0001);;	

/* eMKF: c12 in RAND's MKF User's Guide  */	
%if &btauu    = %str() %then %let  btauu  = %sysevalf(0.1 * &brangeY);;			

/* eMKF: Set cubic and quad precisions so that the coefficients tend to be smaller in magnitude as the degree increases */
%if &bmbeta2 = %str() %then %let bmbeta2 = %sysevalf(0);; 
%if &bpbeta2 = %str() %then %let bpbeta2 = %sysevalf(2.0 * &bpbeta1);; 
%if &bmbeta3 = %str() %then %let bmbeta3 = %sysevalf(0);; 	
%if &bpbeta3 = %str() %then %let bpbeta3 = %sysevalf(4.0 * &bpbeta1);; 		
%if &bbeta2l = %str() %then %let bbeta2l = %sysevalf(0);; 						
%if &bbeta2u = %str() %then %let bbeta2u = %sysevalf(1.5 * &bbeta1u);;
%if &bbeta3l = %str() %then %let bbeta3l = %sysevalf(0);; 	
%if &bbeta3u = %str() %then %let bbeta3u = %sysevalf(2.0 * &bbeta1u);; 

/***************************************************************************************/
/* eMKF: Use data to inform prior parameters for variances in the random variance case */
/***************************************************************************************/

%if %upcase(&brndvars) = YES %then %do;
	%let bqrangeV=0; %let bmedianV=0;
	data _bbjunk;
	run;
	proc means data=_bbdata_ noprint;
	  var _var;
	  output out=_bbjunk median=median qrange=qrange;
	run;
	data _null_;
	 set _bbjunk;
	 call symput("bqrangeV", qrange);
	 call symput("bmedianV", median);
	run;
	%let bqrangeV = %sysevalf(&bqrangeV + 0);
	%let bmedianV = %sysevalf(&bmedianV + 0);
	/* eMKF: Use median for mean and 10 times IQR for standard deviation of sampling variances (inverse gamma prior) */
	%if &bvshape = %str() %then %let bvshape = %sysevalf(2 + ( &bmedianV**2 / ((10 * &bqrangeV)**2) ) );;
	%if &bvscale = %str() %then %let bvscale = %sysevalf((&bvshape - 1)*&bmedianV);;
%end;
%else %do;
	%let bvshape =; 
	%let bvscale =;
%end;
 
/*************************************************************/
/* eMKF: Symbolic array declarations (resolved in proc mcmc) */
/*************************************************************/

/* eMKF: Array structures for BMA weights and prior mixtures */
/* Recall: model 1=indep_cubic, 2=indep_quad, 3=indep_linear, 4=common_cubic, 5=common_quad, 6=common_linear, 7=dropped */
%let wtsline = ; 
%if %upcase(&btype) = BMA_CUBIC %then %do;		  /*  constant shape parameters for Dirichlet */
	%let wtsline = &wtsline array wtsshape[7] (&wshape &wshape &wshape &wshape &wshape &wshape &wshape)%str(;) ; 
	%let wtsline = &wtsline array wts[7]%str(;);  /*  model weights */
%end;
%if %upcase(&btype) = BMA_QUAD %then %do;
	%let wtsline = &wtsline array wtsshape[5] (&wshape &wshape &wshape &wshape &wshape)%str(;) ; 
	%let wtsline = &wtsline array wts[5]%str(;);
%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;
	%let wtsline = &wtsline array wtsshape[3] (&wshape &wshape &wshape)%str(;) ;
	%let wtsline = &wtsline array wts[3]%str(;); 
%end;

/* eMKF: Named 1-dimensional arrays of regression parameters other than intercept */
%let b1line=; %let b2line=; %let b3line=; 
%if %upcase(&btype) = BMA_CUBIC or %upcase(&btype) = BMA_QUAD or %upcase(&btype) = BMA_LINEAR
	%then %let b1line = array b1arr[&g] b1arr1-b1arr&g ;
%if %upcase(&btype) = BMA_CUBIC or %upcase(&btype) = BMA_QUAD
	%then %let b2line = array b2arr[&g] b2arr1-b2arr&g ;
%if %upcase(&btype) = BMA_CUBIC
	%then %let b3line = array b3arr[&g] b3arr1-b3arr&g ;

/* eMKF: Named 1-dimensional arrays of unobserved true states and their means
  (consistent with internal SAS names for random effects in proc mcmc) */
%let etamnarrline = array etamnarr[%eval(&g*&n)];
%let etaarrline   = array etaarr[%eval(&g*&n)];
%let _i = 0; %let _j = 0; 
%do _i = 1 %to &g;
   %do _j = 1 %to &n; 
		%let etamnarrline = &etamnarrline etamn&_j._&_i;
		%let etaarrline   = &etaarrline eta&_j._&_i;
   %end;
%end;

/* eMKF: Named 1-dimensional array of random sampling variances (if applicable) */
%let vline=;
%if %upcase(&brndvars) = YES %then %do;
	%let vline = array varr[&g] varr1-varr&g ; 
	%let vline = &vline%str(;) array vhyp[2] (&bvshape &bvscale) ; /* add array of hyperparameters to pass to UDS */
%end;

/* eMKF: Dynamic array of effective sample sizes (if applicable) to use with read_array */
%let Narrline=;
%if %upcase(&brndvars) = YES %then %let Narrline = array Narr[1] /nosymbols ;

/*****************************************************************/
/* eMKF: Symbolic parameter declarations (resolved in proc mcmc) */
/*****************************************************************/

/* eMKF: Slice sampler, if requested, would apply to parameters for which Gibbs sampling is not available */
%let bslice =%str(;) ;
%if %upcase(&bslicesampler) = YES %then %let bslice = %str(/slice ;);

/* eMKF: Group-specific AR parameters (if applicable) */
%let tauparline=; %let psiparline=; %let tausqparline=; %let rhoparline=; %let tauparline2 =; %let psiparline2 =;
%if %upcase(&bARmodel) = INDEP_AR %then %do;
  	%let psiparline2 = &psiparline2 parms spsi &bslice;			/* SD hyperparameter for mean of psi */
  	%let psiparline2 = &psiparline2 parms mpsi &bslice;			/* mean hyperparameter for mean of psi */
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let psiparline  = &psiparline psi&_i ; 
		%let tauparline  = &tauparline tau&_i ; 
		%let tausqparline= &tausqparline tausq&_i ; 
		%let rhoparline  = &rhoparline rho&_i ; 
    	%let psiparline2 = &psiparline2 parms psi&_i &bslice;	/* Group-specific psi1 through psi&g  */
    	%let tauparline2 = &tauparline2 parms tau&_i &bslice;	/* Group-specific innovation SDs tau1 through tau&g */
	%end; 
%end;
%if %upcase(&bARmodel) = COMMON_AR %then %do;
    %let psiparline2 = parms psi &bslice;						/* Common psi = ln[(1-rho)/(1+rho)] */
    %let tauparline2 = parms tau &bslice;    					/* Common innovation SD tau */
%end;

/* eMKF: Intercepts */
%let aparline=; %let _i = 0;
%do _i = 1 %to &g; 
	%let aparline = &aparline a&_i ; 
%end; 

/* eMKF: Linear, quadratic, and cubic coefficients, as needed */
%let parline = &aparline; %let parline2 = &aparline; %let _i = 0;
%if %upcase(&btype) = BMA_CUBIC %then %do;
	%let parline = &parline b1 b2 b3;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i b2arr&_i b3arr&_i ;
		%let parline2 = &parline2 b1arr&_i b2arr&_i b3arr&_i ;
	%end;
%end;
%if %upcase(&btype) = BMA_QUAD %then %do;
	%let parline = &parline b1 b2;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i b2arr&_i ;
		%let parline2 = &parline2 b1arr&_i b2arr&_i ;
	%end;
%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;
	%let parline = &parline b1;
	%do _i = 1 %to &g; 
		%let parline = &parline b1arr&_i ;
		%let parline2 = &parline2 b1arr&_i ;
	%end;
%end;

/* eMKF: Variance parameters (if applicable) */
%let vparline = ; 	%let _i = 0;
%if %upcase(&brndvars) = YES %then %do;
	%do _i = 1 %to &g; 
		%let vparline = &vparline varr&_i ;
	%end;
%end;

/*************************************/
/* eMKF: UDS parameters declarations */
/*************************************/
%let udsparline = ;

/* model flags updated in a separate UDS block */
%let udsparline = &udsparline parms flg %str(/uds ;);

/* etamnarr updated with the regression coefficients */
%let udsparline = &udsparline parms &parline etamnarr %str(/uds ;);

/* true states updated in a separate UDS block */
%let udsparline = &udsparline parms etaarr %str(/uds ;);

/* variances updated in a separate UDS block (when applicable) */
%if &vparline ^= %str() %then %let udsparline = &udsparline parms &vparline %str(/uds ;);

/**************************************************************************/
/* eMKF: Symbolic prior/hyperprior specifications (resolved in proc mcmc) */
/**************************************************************************/

/* eMKF: Priors for AR parameters */
%let plinetau=; %let plinepsi=; %let hplinempsi=; %let hplinespsi=;
%if %upcase(&bARmodel) = COMMON_AR %then %do; 							/* common AR parameters */
	%let plinepsi = prior psi ~ normal(&bmrho, prec=&bprho); 	
	%let plinetau = prior tau ~ uniform(&btaul, &btauu);     	
%end;
%if %upcase(&bARmodel) = INDEP_AR %then %do; 							/* group-specific AR parameters */
	%let hplinespsi = hyperprior spsi ~ uniform(0.0001,sqrt(1/&bprho)); /* Keep away from zero */
	%let hplinempsi = hyperprior mpsi ~ normal(&bmrho, prec=&bprho);
	%let plinepsi = prior &psiparline ~ normal(mpsi, sd=spsi);
	%let plinetau = prior &tauparline ~ uniform(&btaul, &btauu);   	
%end;

/* eMKF: Prior for mixture weights */
%let plinewts = ;
/*%let plinewts = prior wts ~ dirichlet(wtsshape); */

/* eMKF: Prior for latent variable flg */
%let plineflg = prior flg ~ table(wts);

/* eMKF: Prior for intercepts a1 through a&g */
%let plinea = prior &aparline ~ normal(&bmalpha, prec=&bpalpha);

/****************************************************************/
/* eMKF: Priors for regression parameters other than intercepts */
/****************************************************************/

%let plineb1=; %let plineb2=; %let plineb3=;

%if %upcase(&btype) = BMA_CUBIC %then %do;

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let plineb1 = &plineb1 if flg=4 or flg=5 or flg=6 then lpb1=lpb1+lpdfnorm(b1,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
	%let plineb1 = &plineb1 if (flg=1 or flg=2 or flg=3 or flg=7) and %str(abs(b1-mean(of b1arr1-b1arr&g))>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%let plineb1 = &plineb1 prior b1 ~ general(lpb1)%str(;);

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb1 = &plineb1 if flg=1 or flg=2 or flg=3 then lpb1=lpb1+lpdfnorm(b1arr&_i,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
		%let plineb1 = &plineb1 if flg=7 and %str(abs(b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
		%let plineb1 = &plineb1 if (flg=4 or flg=5 or flg=6) and %str(abs(b1-b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb1 = &plineb1 prior b1arr: ~ general(lpb1)%str(;);

	%let plineb2 = &plineb2 lpb2=0%str(;) ;
	%let plineb2 = &plineb2 if flg=4 or flg=5 then lpb2=lpb2+lpdfnorm(b2,&bmbeta2,sqrt(1/&bpbeta2))%str(;) ;
	%let plineb2 = &plineb2 if flg=6 and %str(abs(b2)> %sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%let plineb2 = &plineb2 if (flg=1 or flg=2 or flg=3 or flg=7) and %str(abs(b2-mean(of b2arr1-b2arr&g))>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%let plineb2 = &plineb2 prior b2 ~ general(lpb2)%str(;);

	%let plineb2 = &plineb2 lpb2=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb2 = &plineb2 if flg=1 or flg=2 then lpb2=lpb2+lpdfnorm(b2arr&_i,&bmbeta2,sqrt(1/&bpbeta2))%str(;) ;
		%let plineb2 = &plineb2 if (flg=3 or flg=7) and %str(abs(b2arr&_i)>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
		%let plineb2 = &plineb2 if (flg=4 or flg=5 or flg=6) and %str(abs(b2-b2arr&_i)>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb2 = &plineb2 prior b2arr: ~ general(lpb2)%str(;);

	%let plineb3 = &plineb3 lpb3=0%str(;) ;
	%let plineb3 = &plineb3 if flg=4 then lpb3=lpb3+lpdfnorm(b3,&bmbeta3,sqrt(1/&bpbeta3))%str(;) ;
	%let plineb3 = &plineb3 if (flg=5 or flg=6) and %str(abs(b3)>%sysevalf(1e-11)) then lpb3=lpb3-%sysevalf(1e15)%str(;) ;
	%let plineb3 = &plineb3 if (flg=1 or flg=2 or flg=3 or flg=7) and %str(abs(b3-mean(of b3arr1-b3arr&g))>%sysevalf(1e-11)) then lpb3=lpb3-%sysevalf(1e15)%str(;) ;
	%let plineb3 = &plineb3 prior b3 ~ general(lpb3)%str(;);

	%let plineb3 = &plineb3 lpb3=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb3 = &plineb3 if flg=1 then lpb3=lpb3+lpdfnorm(b3arr&_i,&bmbeta3,sqrt(1/&bpbeta3))%str(;) ;
		%let plineb3 = &plineb3 if (flg=2 or flg=3 or flg=7) and %str(abs(b3arr&_i)>%sysevalf(1e-11)) then lpb3=lpb3-%sysevalf(1e15)%str(;) ;
		%let plineb3 = &plineb3 if (flg=4 or flg=5 or flg=6) and %str(abs(b3-b3arr&_i)>%sysevalf(1e-11)) then lpb3=lpb3-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb3 = &plineb3 prior b3arr: ~ general(lpb3)%str(;);

%end;
%if %upcase(&btype) = BMA_QUAD %then %do;

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let plineb1 = &plineb1 if flg=3 or flg=4 then lpb1=lpb1+lpdfnorm(b1,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
	%let plineb1 = &plineb1 if (flg=1 or flg=2 or flg=5) and %str(abs(b1-mean(of b1arr1-b1arr&g))>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%let plineb1 = &plineb1 prior b1 ~ general(lpb1)%str(;);

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb1 = &plineb1 if flg=1 or flg=2 then lpb1=lpb1+lpdfnorm(b1arr&_i,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
		%let plineb1 = &plineb1 if flg=5 and %str(abs(b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
		%let plineb1 = &plineb1 if (flg=3 or flg=4) and %str(abs(b1-b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb1 = &plineb1 prior b1arr: ~ general(lpb1)%str(;);

	%let plineb2 = &plineb2 lpb2=0%str(;) ;
	%let plineb2 = &plineb2 if flg=3 then lpb2=lpb2+lpdfnorm(b2,&bmbeta2,sqrt(1/&bpbeta2))%str(;) ;
	%let plineb2 = &plineb2 if flg=4 and %str(abs(b2)>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%let plineb2 = &plineb2 if (flg=1 or flg=2 or flg=5) and %str(abs(b2-mean(of b2arr1-b2arr&g))>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%let plineb2 = &plineb2 prior b2 ~ general(lpb2)%str(;);

	%let plineb2 = &plineb2 lpb2=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb2 = &plineb2 if flg=1 then lpb2=lpb2+lpdfnorm(b2arr&_i,&bmbeta2,sqrt(1/&bpbeta2))%str(;) ;
		%let plineb2 = &plineb2 if (flg=2 or flg=5) and %str(abs(b2arr&_i)>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
		%let plineb2 = &plineb2 if (flg=3 or flg=4) and %str(abs(b2-b2arr&_i)>%sysevalf(1e-11)) then lpb2=lpb2-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb2 = &plineb2 prior b2arr: ~ general(lpb2)%str(;);

%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let plineb1 = &plineb1 if flg=2 then lpb1=lpb1+lpdfnorm(b1,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
	%let plineb1 = &plineb1 if (flg=1 or flg=3) and %str(abs(b1-mean(of b1arr1-b1arr&g))>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%let plineb1 = &plineb1 prior b1 ~ general(lpb1)%str(;);

	%let plineb1 = &plineb1 lpb1=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
	    %let plineb1 = &plineb1 if flg=1 then lpb1=lpb1+lpdfnorm(b1arr&_i,&bmbeta1,sqrt(1/&bpbeta1))%str(;) ;
		%let plineb1 = &plineb1 if flg=3 and %str(abs(b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
		%let plineb1 = &plineb1 if flg=2 and %str(abs(b1-b1arr&_i)>%sysevalf(1e-11)) then lpb1=lpb1-%sysevalf(1e15)%str(;) ;
	%end;
	%let plineb1 = &plineb1 prior b1arr: ~ general(lpb1)%str(;);

%end;

/* eMKF: Prior for variance parameters */
%let plinev=;
%if %upcase(&brndvars) = YES %then %let plinev = prior varr: ~ igamma(&bvshape, scale=&bvscale);

/******************************************************************************/
/* eMKF: Symbolic initialization for model parameters (resolved in proc mcmc) */
/******************************************************************************/

/* eMKF: Initial values for AR parameters */
%let initlinetau = ; %let initlinepsi = ; %let _i = 0;
%if %upcase(&bARmodel) = COMMON_AR %then %do; 	/*common AR parameters */
	%let initlinepsi = psi = &bmrho + sqrt(1/&bprho)*rand('normal');
	%let initlinetau = tau = rand('uniform', &btaul, &btauu); 
%end;
%if %upcase(&bARmodel) = INDEP_AR %then %do; /* Group-specific AR parameters */
    %let initlinepsi = &initlinepsi spsi = rand('uniform', 0.0001, sqrt(1/&bprho))%str(;) ;		
    %let initlinepsi = &initlinepsi mpsi = &bmrho + sqrt(1/&bprho)*rand('normal')%str(;) ;
	%do _i = 1 %to &g;
		%let initlinepsi = &initlinepsi psi&_i=mpsi+spsi*rand('normal')%str(;) ;
		%let initlinetau = &initlinetau tau&_i=rand('uniform',&btaul,&btauu)%str(;) ;
	%end;	
%end;

/* eMKF: Dimensionality for do loops to initialize mixture parameters */
%let _ll=0;
%if %upcase(&btype) = BMA_CUBIC  %then %let _ll = 7;
%if %upcase(&btype) = BMA_QUAD 	 %then %let _ll = 5;
%if %upcase(&btype) = BMA_LINEAR %then %let _ll = 3;
%let _ll = %eval(0+&_ll);

/* eMKF: Initial values for mixture weights */
%let initlinewts = wtssum = 0%str(;) ;
%let _l=0; 
%do _l=1 %to &_ll; 
	/*%let initlinewts = &initlinewts wts[&_l] = rand('gamma', wtsshape[&_l])%str(;) ; */
	%let initlinewts = &initlinewts wts[&_l] = wtsshape[&_l]%str(;) ; 
	%let initlinewts = &initlinewts wtssum = wtssum + wts[&_l]%str(;) ;
%end;

/* eMKF: Rescale wts to sum to one */
%let _l=0;
%do _l=1 %to &_ll;		
	%let initlinewts = &initlinewts wts[&_l] = wts[&_l]/wtssum%str(;) ;
%end;

/* eMKF: Initial values for latent variable flg */
%let initlineflg = ; %let _l=0;
%do _l=1 %to &_ll; 
	%let initlineflg = &initlineflg %str(,) wts[&_l]; 
%end;
%let initlineflg = flg = rand('table' &initlineflg);

/* eMKF: Initial/constant values for prior mean vector mbetag and precision matrix Dbetag for use with matrix operations */
%let initmbeta = call zeromatrix(Dbetag);	
%let initmbeta = &initmbeta%str(;) mbetag[1,1] = &bmalpha%str(;) Dbetag[1,1] = &bpalpha;
%if %upcase(&btype) = BMA_CUBIC %then %do;
	%let initmbeta = &initmbeta%str(;) mbetag[2,1] = &bmbeta1%str(;) Dbetag[2,2] = &bpbeta1;
	%let initmbeta = &initmbeta%str(;) mbetag[3,1] = &bmbeta2%str(;) Dbetag[3,3] = &bpbeta2;
	%let initmbeta = &initmbeta%str(;) mbetag[4,1] = &bmbeta3%str(;) Dbetag[4,4] = &bpbeta3;
%end;
%if %upcase(&btype) = BMA_QUAD %then %do;
	%let initmbeta = &initmbeta%str(;) mbetag[2,1] = &bmbeta1%str(;) Dbetag[2,2] = &bpbeta1;
	%let initmbeta = &initmbeta%str(;) mbetag[3,1] = &bmbeta2%str(;) Dbetag[3,3] = &bpbeta2;
%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;
	%let initmbeta = &initmbeta%str(;) mbetag[2,1] = &bmbeta1%str(;) Dbetag[2,2] = &bpbeta1;
%end;

/* eMKF: Initial values for intercepts */
%let initlinea=; %let _i=0;
%do _i = 1 %to &g; 
	%let initlinea = &initlinea a&_i = &bmalpha + sqrt(1/&bpalpha)*rand('normal')%str(;) ;
%end;

/*********************************************************************************/
/* eMKF: Initial values for regression parameters and applicable hyperparameters */
/*********************************************************************************/

%let initlineb1=; %let initlineb2=; %let initlineb3=; 

%if %upcase(&btype) = BMA_CUBIC %then %do;

	%let initlineb3 = &initlineb3 if flg=4 then b3=&bmbeta3+sqrt(1/&bpbeta3)*rand('normal')%str(;) ;
	%let initlineb3 = &initlineb3 if flg=5 or flg=6 then b3=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb3 = &initlineb3 if flg=1 then b3arr&_i=&bmbeta3+sqrt(1/&bpbeta3)*rand('normal')%str(;) ;
		%let initlineb3 = &initlineb3 if flg=2 or flg=3 or flg=7 then b3arr&_i=0%str(;) ;
		%let initlineb3 = &initlineb3 if flg=4 or flg=5 or flg=6 then b3arr&_i=b3%str(;) ;
	%end;
	%let initlineb3 = &initlineb3 %str(if flg=1 or flg=2 or flg=3 or flg=7 then b3=mean(of b3arr1-b3arr&g);) ;

	%let initlineb2 = &initlineb2 if flg=4 or flg=5 then b2=&bmbeta2+sqrt(1/&bpbeta2)*rand('normal')%str(;) ;
	%let initlineb2 = &initlineb2 if flg=6 then b2=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb2 = &initlineb2 if flg=1 or flg=2 then b2arr&_i=&bmbeta2+sqrt(1/&bpbeta2)*rand('normal')%str(;) ;
		%let initlineb2 = &initlineb2 if flg=3 or flg=7 then b2arr&_i=0%str(;) ;
		%let initlineb2 = &initlineb2 if flg=4 or flg=5 or flg=6 then b2arr&_i=b2%str(;) ;
	%end;
	%let initlineb2 = &initlineb2 %str(if flg=1 or flg=2 or flg=3 or flg=7 then b2= mean(of b2arr1-b2arr&g);) ;

	%let initlineb1 = &initlineb1 if flg=4 or flg=5 or flg=6 then b1=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb1 = &initlineb1 if flg=1 or flg=2 or flg=3 then b1arr&_i=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
		%let initlineb1 = &initlineb1 if flg=7 then b1arr&_i=0%str(;) ;
		%let initlineb1 = &initlineb1 if flg=4 or flg=5 or flg=6 then b1arr&_i=b1%str(;) ;
	%end;
	%let initlineb1 = &initlineb1 %str(if flg=1 or flg=2 or flg=3 or flg=7 then b1=mean(of b1arr1-b1arr&g);) ;

%end;
%if %upcase(&btype) = BMA_QUAD %then %do;

	%let initlineb2 = &initlineb2 if flg=3 then b2=&bmbeta2+sqrt(1/&bpbeta2)*rand('normal')%str(;) ;
	%let initlineb2 = &initlineb2 if flg=4 then b2=0%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb2 = &initlineb2 if flg=1 then b2arr&_i=&bmbeta2+sqrt(1/&bpbeta2)*rand('normal')%str(;) ;
		%let initlineb2 = &initlineb2 if flg=2 or flg=5 then b2arr&_i=0%str(;) ;
		%let initlineb2 = &initlineb2 if flg=3 or flg=4 then b2arr&_i=b2%str(;) ;
	%end;
	%let initlineb2 = &initlineb2 %str(if flg=1 or flg=2 or flg=5 then b2=mean(of b2arr1-b2arr&g);) ;

	%let initlineb1 = &initlineb1 if flg=3 or flg=4 then b1=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb1 = &initlineb1 if flg=1 or flg=2  then b1arr&_i=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
		%let initlineb1 = &initlineb1 if flg=5 then b1arr&_i=0%str(;) ;
		%let initlineb1 = &initlineb1 if flg=3 or flg=4 then b1arr&_i=b1%str(;) ;
	%end;
	%let initlineb1 = &initlineb1 %str(if flg=1 or flg=2 or flg=5 then b1=mean(of b1arr1-b1arr&g);) ;

%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;

	%let initlineb1 = &initlineb1 if flg=2 then b1=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
	%let _i = 0;
	%do _i = 1 %to &g; 
		%let initlineb1 = &initlineb1 if flg=1 then b1arr&_i=&bmbeta1+sqrt(1/&bpbeta1)*rand('normal')%str(;) ;
		%let initlineb1 = &initlineb1 if flg=3 then b1arr&_i=0%str(;) ;
		%let initlineb1 = &initlineb1 if flg=2 then b1arr&_i=b1%str(;) ;
	%end;
	%let initlineb1 = &initlineb1 %str(if flg=1 or flg=3 then b1=mean(of b1arr1-b1arr&g);) ;

%end;

/* eMKF: Initial values for unobserved true states predictions given regression parameters */
%let _i = 0; %let _j = 0; 
%if %upcase(&btype) = BMA_CUBIC %then %do;
  %do _i = 1 %to &g; 
  	%local initetamnarr&_i; /*eMKF: broken up into one macro variable per group instead of single combined macro variable to avoid max length error (65534) */
  	%do _j = 1 %to &n; 
	  %let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1arr&_i+X[&_j,3]*b2arr&_i+X[&_j,4]*b3arr&_i%str(;) ;
    %end;
  %end;
%end;
%if %upcase(&btype) = BMA_QUAD %then %do;
	%do _i = 1 %to &g; 
	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
	  		%let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,1]*b1arr&_i+X[&_j,2]*b2arr&_i%str(;) ;
		%end;
	%end;
%end;
%if %upcase(&btype) = BMA_LINEAR %then %do;
	%do _i = 1 %to &g; 
	    %local initetamnarr&_i;
  		%do _j = 1 %to &n; 
			%let initetamnarr&_i = &&initetamnarr&_i etamnarr[%eval((&_i-1)*&n+&_j)]=X[&_j,1]*a&_i+X[&_j,2]*b1arr&_i%str(;) ;
		%end;
	%end;
%end;

/* eMKF: Initial values for variance parameters from igamma(&bvshape, scale=&bvscale) (if applicable) */
%let initlinevarr = ; 
%if %upcase(&brndvars) = YES %then %do; 
	%let _i = 0;
	%do _i = 1 %to &g;
		%let initlinevarr = &initlinevarr varr&_i=1/rand('gamma',&bvshape,1/&bvscale)%str(;) ;
	%end;
%end;

/* eMKF: temporary dataset for building group-specific design matrix */
data _bbdata1_;
  set _bbdata_(where=(_group_ = 1));
  keep &brtm.0 &brtm.1 &brtm.2 &brtm.3;
run;

/* eMKF: dimensionality */
%let p = 1;
%if %upcase(&btype) = BMA_LINEAR %then %let p = 2;
%if %upcase(&btype) = BMA_QUAD   %then %let p = 3;
%if %upcase(&btype) = BMA_CUBIC  %then %let p = 4;
%let p = %eval(0+&p);

/* eMKF: applicable read_array statement for the design matrix */
%let rcXline = ;
%if &p = 2 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'));
%if &p = 3 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'), resolve('&brtm.2'));
%if &p = 4 %then %let rcXline = rcX = read_array('_bbdata1_', Xarr, resolve('&brtm.0'), resolve('&brtm.1'), resolve('&brtm.2'), resolve('&brtm.3'));

/* eMKF: applicable read_array statement for the effective sample sizes */
%let rcNline = ;
%if %upcase(&brndvars) = YES %then %let	rcNline = rcN = read_array('_bbdata_', Narr, '_n');

/*************************************************************************************/
/*eMKF: Applicable UDS statements - see macros gibbs_uds_compile_** for definitions  */
/*      These will be applied in the order provided here, and after any M-H samplers */
/*************************************************************************************/
%let udsline = ;
 
/*eMKF: UDS statement for model flag */
%if %upcase(&btype) = BMA_CUBIC %then
	%let udsline = &udsline uds FP_bmac(flg, wts, a, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = BMA_QUAD %then
	%let udsline = &udsline uds FP_bmaq(flg, wts, a, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;
%if %upcase(&btype) = BMA_LINEAR %then
	%let udsline = &udsline uds FP_bmal(flg, wts, a, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr)%str(;) ;

/*eMKF: UDS statement for regression coefficients  */
/*      The pseudo-parameter etamnarr is also updated in those subroutines to hold the updated regression predictions */
%if %upcase(&btype) = BMA_CUBIC %then
	%let udsline = &udsline uds CP_bmac(a, b1arr, b2arr, b3arr, b1, b2, b3, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr, flg)%str(;) ;
%if %upcase(&btype) = BMA_QUAD %then
	%let udsline = &udsline uds CP_bmaq(a, b1arr, b2arr, b1, b2, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr, flg)%str(;) ;
%if %upcase(&btype) = BMA_LINEAR %then
	%let udsline = &udsline uds CP_bmal(a, b1arr, b1, etamnarr, mbetag, Dbetag, rhoarr, nuarr, rts, X, Yarr, Sarr, flg)%str(;) ;

/*eMKF: UDS statement for true states etaarr */
%let udsline = &udsline uds EP(etaarr, etamnarr, rhoarr, nuarr, rts, Yarr, Sarr)%str(;) ;

/*eMKF: UDS statement for variances (if applicable) */
%if %upcase(&brndvars) = YES %then 
	%let udsline = &udsline uds RP(varr, vhyp, Sarr, Narr)%str(;) ;

/* eMKF: library location for pre-compiled UDS subroutines */
options cmplib = &bcmploc;

/* eMKF: Options will be the proc mcmc defaults if not specified by the user */
%let optionline=;
%if &bseed ^= %str()  %then %let optionline = &optionline seed 		= %eval(0+&bseed);;
%if &bmaxt ^= %str()  %then %let optionline = &optionline maxtune 	= %eval(0+&bmaxt);;
%if &btune ^= %str()  %then %let optionline = &optionline ntu 		= %eval(0+&btune);;
%if &bburn ^= %str()  %then %let optionline = &optionline nbi 		= %eval(0+&bburn);;
%if &biter ^= %str()  %then %let optionline = &optionline nmc 		= %eval(0+&biter);;
%if &bthin ^= %str()  %then %let optionline = &optionline thin 		= %eval(0+&bthin);;
%if &batol ^= %str()  %then %let optionline = &optionline accepttol = %sysevalf(&batol);;
%if &bttol ^= %str()  %then %let optionline = &optionline targaccept = %sysevalf(&bttol);;
%if &bprcov ^= %str() %then %let optionline = &optionline propcov 	= &bprcov;
%if &binit ^= %str()  %then %let optionline = &optionline init 		= &binit;

/* eMKF: Disable summary statistics if not requested by the user */
%if %upcase(&bprint) ^= YES %then %let optionline = &optionline stats = none;

/* eMKF: Diagnostics plots and ODS graphics enabled if requested by the user */
%if %upcase(&bplot) = YES %then %do; 
	%let optionline = &optionline plots = all;
	ods graphics on;
%end;
%else %let optionline = &optionline plots = none;

/* eMKF: Add jointmodel option (log-likelihood constructed using stored arrays) */
%let optionline = &optionline jointmodel;

/* eMKF: Monitor selected model parameters */
%let monitorline = /*wts*/ flg &parline2 etaarr &vparline;
%if %upcase(&bARmodel) = INDEP_AR %then %let monitorline = spsi mpsi &tausqparline &rhoparline &monitorline ;
%if %upcase(&bARmodel) = COMMON_AR %then %let monitorline = tausq rho &monitorline ;

/* eMKF: Empty dataset to pass to proc mcmc: data from _bbdata_ will be read directly into arrays */
data _bb_;
run;

/* eMKF: Call proc mcmc using the above customizations  */
%put Call to PROC MCMC initiated; %let _i = 0;
proc mcmc data=_bb_ outpost= &blog monitor = ( &monitorline ) &optionline;;	

	  %if %upcase(&bprint) ^=YES and %upcase(&bplot) ^=YES 	/* Disable output tables and plots as applicable */
		%then ods select none;;

	  /**********************/
	  /* Array declarations */
	  /**********************/
	  array rts[&n] (&_brtimess); 	 						/* constant array with real times */
	  array Xarr[1]						   	    /nosymbols;	/* dynamic array for predictors to read in from dataset */
	  array Yarr[1]			   	  			    /nosymbols;	/* dynamic array for _y from dataset */
	  array Sarr[1]			   	   			    /nosymbols;	/* dynamic array for _var from dataset */
	  &Narrline;;											/* dynamic array for _n from dataset (if applicable) */
	  array X[&n, &p];										/* design matrix to use in matrix multiplication */
	  array mbetag[&p, 1];									/* prior mean vector for betas (assumed common across groups) */
	  array Dbetag[&p, &p];									/* diagonal prior precision matrix for betas (assumed common across groups) */
	  %if %upcase(&bARmodel) = INDEP_AR %then %do;			/* AR-related parameters in the group-specific random effects model */
	 	  array psi[&g] psi1-psi&g;							/* group-specific psi = ln[(1-rho)/(1+rho)] */
	  	  array tau[&g] tau1-tau&g;				    		/* group-specific innovation SD tau */
	  	  array tausq[&g] tausq1-tausq&g;					/* squares of group-specific innovation SD tau */
	  	  array rho[&g] rho1-rho&g;							/* reverse-transformation for rho */
	  	  array nu[&g] nu1-nu&g;							/* innovation variance parameters under stationarity */
	  	  array dg[&g] dg1-dg&g;							/* determinants of AR variance-covariance matrices */
	  %end;
	  array rhoarr[&g]; 									/* temporary 1-dimensional array with group-specific parameters rho  */
	  array nuarr[&g]; 										/* temporary 1-dimensional array with group-specific parameters nu */
	  &wtsline;;											/* mixture weights and model flags */
	  array a[&g] a1-a&g;									/* named 1-dimensional array of group-specific intercepts */
	  &b1line;;												/* named 1-dimensional array of group-specific linear coefficients (if requested) */
 	  &b2line;;                     						/* named 1-dimensional array of group-specific quad coefficients (if requested) */
	  &b3line;; 											/* named 1-dimensional array of group-specific cubic coefficients (if requested) */
	  &etamnarrline;;										/* named 1-dimensional array etamnarr (gxn) for predictions from regression */
	  &vline;;												/* named 1-dimensional array of group-specific variance parameters (if requested) */
	  &etaarrline;;											/* named 1-dimensional array etaarr (gxn) for unobserved true states */

	  begincnst;

		  /*****************/
	  	  /* Design matrix */
		  /*****************/
	  	  &rcXline;;										/* read in dynamic array of predictors Xarr */  
		  call zeromatrix(X);								/* initialize design matrix X to all zeroes */
		  do i = 1 to &n;								    /* Xarr is a 2-dimensional array for p > 1 */
		  	  do m = 1 to &p;
				  X[i,m] = Xarr[i,m];						
			  end;
		  end;

		  /**********************/
		  /* Group sample means */
		  /**********************/
		  rcY = read_array('_bbdata_', Yarr, '_y');			/* read in 1-dimensional array of _y from dataset */

		  /**********************/
	  	  /* Sampling variances */
		  /**********************/
		  rcS = read_array('_bbdata_', Sarr, '_var');		/* read in 1-dimensional array of _var from dataset */

		  /******************************************/
		  /* Effective sample sizes (if applicable) */
		  /******************************************/
		  &rcNline;;										/* read in 1-dimensional array of _n from dataset (if applicable) */

		  /******************/
	 	  /* Initialization */
	  	  /******************/
		  call streaminit(%eval(0+&bseed));					/* set seed */

		  &initlinepsi;;									/* initialize psi = ln[(1-rho)/(1+rho)] */
		  &initlinetau;;									/* initialize innovation SD tau */
		  %if %upcase(&bARmodel) = COMMON_AR %then %do;		/* common AR parameters across groups */
			  tausq = tau**2;					 		 	/* track tau-squared */
	  		  rho = (1-exp(psi))/(1+exp(psi)); 		 		/* reverse-transformation for rho */
	  		  nu = tausq/(1-rho**2);			 		 	/* innovation variance parameter under stationarity */
			  dg = nu**&n;									/* recursive formula for determinant of Vgamma (assuming 2+ points) */
			  do i = 2 to &n;								
				  dg = dg*(1-(rho**(2*(rts[i]-rts[i-1])))); 
			  end;
			  if abs(rho) ge 1 or dg= . or dg le 0 then do; /* guard against numerical singularities */
				  rho = 0;
				  nu = tausq;
				  dg = nu**&n;
			  end;
			  do k=1 to &g;									/* temp parameter arrays (e.g., to pass to UDS subroutines) */
				  rhoarr[k] = rho;
			  	  nuarr[k] = nu;
			  end;
		  %end;
		  %if %upcase(&bARmodel) = INDEP_AR %then %do;		/* independent AR parameters across groups */
			  do k=1 to &g;
				  tausq[k] = tau[k]**2;		 	
		  		  rho[k] = (1-exp(psi[k]))/(1+exp(psi[k]));
		  		  nu[k] = tausq[k]/(1-rho[k]**2);
				  dg[k] = nu[k]**&n;
				  do i = 2 to &n;
					  dg[k] = dg[k]*(1-(rho[k]**(2*(rts[i] - rts[i-1]))));
				  end;
				  if abs(rho[k]) ge 1 or dg[k] = . or dg[k] le 0 then do;
					  rho[k] = 0;             
					  nu[k] = tausq[k]; 		 
				      dg[k] = nu[k]**&n;
				  end;
				  rhoarr[k] = rho[k];
			  	  nuarr[k] = nu[k];
			  end;
		  %end;

		  &initmbeta;;							 			/* initialize constant vector mbetag and constant matrix Dbetag */

		  &initlinewts;;									/* initialize model weights and related arrays */
		  &initlineflg;;									/* initialize model flags */

		  &initlinea;;										/* initialize intercepts */
		  &initlineb1;;										/* initialize linear coefficients and related arrays (if applicable) */
		  &initlineb2;;										/* initialize quad coefficients and related arrays (if applicable) */	
		  &initlineb3;;										/* initialize cubic coefficients and related arrays (if applicable) */

          %do _i = 1 %to &g; 
		      &&initetamnarr&_i;;						    /* initialize conditional mean for true states  */ 
		  %end;

		  do k = 1 to &g; 			  						/* initialize etaarr using Markov property of AR process */
			  etaarr[(k-1)*&n+1] = etamnarr[(k-1)*&n+1] + 
							sqrt(nuarr[k])*rand('normal'); 	/* first timepoint from stationary distribution of AR process */

		  	  do i = 2 to &n; 								/* subsequent timepoints from implied conditional distributions */
			      etaarr[(k-1)*&n+i] = etamnarr[(k-1)*&n+i] + 
							((rhoarr[k]**(rts[i] - rts[i-1]))*(etaarr[(k-1)*&n+i-1] - etamnarr[(k-1)*&n+i-1])) + 
							sqrt(nuarr[k]*(1-(rhoarr[k]**(2*(rts[i] - rts[i-1])))))*rand('normal');
		  	  end;
		  end;

		  &initlinevarr;;									/* initialize sampling variances (if applicable) */

	  endcnst;

	  /*******************/
	  /* UDS declaration */
	  /*******************/
	  &udsline;;											/* Gibbs sampling done in the order specified in udsline */
	  														/* Per SAS documentation, parameters that use M-H will be sampled first */
	  /**************************/
	  /* Parameter declarations */
	  /**************************/
	  &psiparline2;;										/* psi = ln[(1-rho)/(1+rho)] and any hyperparameters */
	  &tauparline2;;										/* innovation SD tau */
	  /* parms wts &bslice; */ 								/* Dirichlet mixture weights */
	  &udsparline;;											/* UDS parameter blocks, one for each Gibbs sampler */

	  beginnodata;

	  	  /********************/
	  	  /* Prior statements */
	  	  /********************/
		  &hplinespsi;;										/* SD hyper-prior for mean of psi = ln[(1-rho)/(1+rho)] (if applicable) */
		  &hplinempsi;;										/* Mean hyper-prior for mean of psi = ln[(1-rho)/(1+rho)] (if applicable) */
		  &plinepsi;;										/* prior for psi = ln[(1-rho)/(1+rho)] */
		  &plinetau;;										/* prior for innovation SD tau */
		  %if %upcase(&bARmodel) = COMMON_AR %then %do;		/* AR parameters in the common case */
			  tausq = tau**2;					 			/* track tau-squared */
	  		  rho = (1-exp(psi))/(1+exp(psi)); 		 		/* reverse-transformation for rho */
	  		  nu = tausq/(1-rho**2);			 			/* innovation variance parameter under stationarity */
			  dg = nu**&n;									/* recursive formula for determinant of Vgamma (assuming 2+ points) */
			  do i = 2 to &n;								
				  dg = dg*(1-(rho**(2*(rts[i]-rts[i-1]))));
			  end;
			  if abs(rho) ge 1 or dg= . or dg le 0 then do; /* guard against numerical singularities */
				  rho = 0;
				  nu = tausq;
				  dg = nu**&n;
			  end;
			  do k=1 to &g;									/* parameter arrays to pass to UDS subroutines */
				  rhoarr[k] = rho;
			  	  nuarr[k] = nu;
			  end;
		  %end;
	      %if %upcase(&bARmodel) = INDEP_AR %then %do;	 	/* Group-specific AR parameters */
			do k = 1 to &g; 
			    tausq[k] = tau[k]**2;
		  	    rho[k] = (1-exp(psi[k]))/(1+exp(psi[k]));
		  	    nu[k] = tausq[k]/(1-rho[k]**2);
			    dg[k] = nu[k]**&n;
			    do i = 2 to &n;
				    dg[k] = dg[k]*(1-(rho[k]**(2*(rts[i] - rts[i-1]))));
			    end;
			    if abs(rho[k]) ge 1 or dg[k] = . or dg[k] le 0 then do;
				    rho[k] = 0;             
				    nu[k] = tausq[k]; 		 
				    dg[k] = nu[k]**&n;
			    end;
				rhoarr[k] = rho[k];
			  	nuarr[k] = nu[k];		
			end;
		  %end;
 
		  /*&plinewts;;	*/									/* Dirichlet prior for model weights */
		  &plineflg;;										/* Discrete prior for model flag */

	  	  &plinea;; 							 			/* prior for intercepts */
	  	  &plineb1;;							 			/* conditional prior for linear coefficients given model flag (if applicable) */
	  	  &plineb2;;							 			/* conditional prior for quadratic coefficients given model flag (if applicable) */
	  	  &plineb3;;							 			/* conditional prior for cubic coefficients given model flag (if applicable) */

		  prior etamnarr ~ general(0);						/* pseudo-parameters etamnarr do no contribute to prior */

		  lpr = 0; 											/* calculation of log-prior for etaarr from univariate conditionals */
		  do k = 1 to &g;
		  	  lpr = lpr + lpdfnorm(etaarr[(k-1)*&n+1], 		/* etaarr is updated in the UDS call for the true states */
								   etamnarr[(k-1)*&n+1], 	/* etamnarr is updated in the UDS call for the regression coefficients */
								   sqrt(nuarr[k]));			/* first timepoint from stationary distribution of AR process */
		  	  do i = 2 to &n; 								/* subsequent timepoints from implied conditional distributions */
			      lpr = lpr + lpdfnorm(etaarr[(k-1)*&n+i], 
									   etamnarr[(k-1)*&n+i] + 
									    (rhoarr[k]**(rts[i] - rts[i-1]))*(etaarr[(k-1)*&n+i-1] - etamnarr[(k-1)*&n+i-1]), 
									   sqrt(nuarr[k]*(1-(rhoarr[k]**(2*(rts[i] - rts[i-1])))))); 
		  	  end;
		  end;

		  prior etaarr ~ general(lpr);						/* prior for unobserved true states */

		  &plinev;;								 			/* Inverse gamma prior for sampling variances (if applicable) */

		  /********************************/
	  	  /* Loglikelihood calculation(s) */
	  	  /********************************/
	  	  lp = 0;		  									/* log of joint distribution of sample means */
		  do k = 1 to &g*&n;
		      lp = lp + lpdfnorm(Yarr[k], etaarr[k], sqrt(Sarr[k]));
		  end;
		  %if %upcase(&brndvars) = YES %then %do;		  	/* log of joint distribution of sample variances (if applicable) */
			  do k = 1 to &g;
				  do j = 1 to &n;
					  lp = lp + lpdfgamma(Sarr[(k-1)*&n+j],
										  (Narr[(k-1)*&n+j]-1)/2,
										  (2*varr[k])/(Narr[(k-1)*&n+j]-1)); /* varr is updated in the UDS call for the variances */
			 	  end;
			  end;
		  %end;

	  endnodata;

	  /*******************/
	  /* Model statement */
	  /*******************/
	  model general(lp);

run;

/* eMKF: Re-enable output tables and plots */
%if %upcase(&bprint) ^= YES and %upcase(&bplot) ^= YES %then ods select all;;

/* eMKF: Disable ODS graphics */
%if %upcase(&bplot) = YES %then ods graphics off;;

%put Call to PROC MCMC concluded;

/*eMKF: Keep only the desired columns in the posterior log dataset */
%if %upcase(&bARmodel) = INDEP_AR %then %do;
	data &blog;
	  merge &blog(drop= etamn: Log: spsi mpsi &tauparline &psiparline
						%if %upcase(&btype) = BMA_CUBIC %then b1 b2 b3; 
  						%if %upcase(&btype) = BMA_QUAD %then b1 b2; 
						%if %upcase(&btype) = BMA_LINEAR %then b1; 
				 ) 
			&blog(keep = spsi mpsi);
	run;
%end;
%if %upcase(&bARmodel) = COMMON_AR %then %do;
	data &blog;
		set &blog(drop= etamn: Log: tau psi
						%if %upcase(&btype) = BMA_CUBIC %then b1 b2 b3; 
  						%if %upcase(&btype) = BMA_QUAD %then b1 b2; 
						%if %upcase(&btype) = BMA_LINEAR %then b1; 
				 );
	run;
%end;

/****************************************************/
/* eMKF: Reverse-transform regression coefficients  */
/****************************************************/

data _blogc2_ _tblogc2_ _tblogc_  ;
run;

%let _i = 0; 

%if %upcase(&borpoly) = YES %then %do;

	/* eMKF: order columns by group */
	data _blogc2_;
  	  retain  Iteration 
			  %do _i=1 %to &g;
				  a&_i
				  %if %upcase(&btype) = BMA_LINEAR %then b1arr&_i;
				  %if %upcase(&btype) = BMA_QUAD   %then b1arr&_i b2arr&_i;
				  %if %upcase(&btype) = BMA_CUBIC  %then b1arr&_i b2arr&_i b3arr&_i;
			  %end;
	  ;
	  set &blog(keep = Iteration a: b:);
	run;

	/* eMKF: block diagonal by group */
	%let oPPmat = ; %let _i = 0;
	%do _i=1 %to &g; 
		%if &_i = 1 %then %let oPPmat = block( oP ;
		%if &_i > 1 and &_i < &g %then %let oPPmat = &oPPmat , block ( oP ;
		%if &_i = &g and &g > 1  %then %let oPPmat = &oPPmat , oP %sysfunc(repeat( %str(%)), &g-2));
		%if &_i = &g and &g = 1  %then %let oPPmat = &oPPmat );
	%end;

	%let _i = 0;

	/* eMKF: call proc iml to perform matrix multiplication */
	proc iml;

		use _oPmat_;
		read all into oP; close _oPmat_;
		oP = oP[1:&p, 1:&p];
		oPP = &oPPmat;;

		varNames = {"Iteration"};
		%if %upcase(&btype) = BMA_LINEAR %then %do;
			%do _i = 1 %to &g;
				varNames = varNames || {"a&_i"} || {"b1arr&_i"};
			%end;
		%end;
		%if %upcase(&btype) = BMA_QUAD   %then %do; 
			%do _i = 1 %to &g;
				varNames = varNames || {"a&_i"} || {"b1arr&_i"} || {"b2arr&_i"};
			%end;
		%end;
		%if %upcase(&btype) = BMA_CUBIC  %then %do;
			%do _i = 1 %to &g;
	 			varNames = varNames || {"a&_i"} || {"b1arr&_i"} || {"b2arr&_i"} || {"b3arr&_i"};
			%end;
		%end;

		use _blogc2_;
		read all into oB;
		close _blogc2_;

		oB1 = oB[,1];
		oB = T(oB[,2:ncol(oB)]);
		oBB = oPP * oB;
		oBB = oB1 || T(oBB);

		create _tblogc2_ var varNames;
		append from oBB;
		close _tblogc2_;

	quit;

	/* eMKF: re-order columns as they were initially from PROC MCMC */
	data _tblogc_;
  	  retain  Iteration a1-a&g 
			  %if %upcase(&btype) = BMA_LINEAR  %then b1arr1-b1arr&g ; 
			  %if %upcase(&btype) = BMA_QUAD    %then b1arr1-b1arr&g b2arr1-b2arr&g ; 
			  %if %upcase(&btype) = BMA_CUBIC   %then b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g ; 
	  ;
	  set _tblogc2_;
	run;

	/* eMKF: merge into &blog */
	data &blog;
	  merge &blog(keep = Iteration flg /*wts*/)
	        _tblogc_
			&blog(drop = a: b: )
	  ;
	  by Iteration;
	run;
	  
%end;

/* eMKF: clean-up */
proc datasets nolist;
 delete _bbdata_ _bbdata1_ _bb_ _bfreqg_ _bfreqn_ _bbjunk _oXmat_ _oPmat_ _blogc2_ _tblogc2_ _tblogc_;
run ;
quit;

%mend;

data _null_;
run;




/*HTRP macro
Macro defined 03-02-2007
It allows for the estimation of parameters from the non-linear mixed effect model

Modified for eMKF in 2023 Q1-Q2 by Makram Talih.
In eMKF, this assumes data has been pre-formatted using macro reformat.

data 	: name of the dataset.
outcome : outcome of interest. 
se      : standard error of the outcome. 
time    : time variable.
by      : allows models to run for multiple strata at the same time and can be used for simulations
xtrakeep: Any variable one wants to keep in the data while runing models: weights, ... (eMKF: could be used to retain labels for multiyear data)
orpoly  : (eMKF) YES (default) for pre-transforming the design matrix using SAS IML orpol function. NO for "raw" polynomials.
          If YES, regression coefficients and their SEs will be reverse-transformed prior to macro end.
_rho_   : the value of the true rho that generated the data, if known. If not given, it will be estimated
_tausq_ : the value of the true tau-square that generated the data, if known. If not given, it will be estimated
bvalue  : Assumption about the bvalue: eMKF options are the following:
             indep_cubic	: The values of the parameters b1, b2, and b3 are computed for each group
             indep_quad		: b3=0. The values of the parameters b1 and b2 are computed for each group
             indep_linear   : (DEFAULT) b3=0 and b2=0. The value of the slope b1 is computed for each group
             common_cubic	: The values of each of the parameters b1, b2, and b3 are assumed to be the same across groups
             common_quad	: b3=0. The values of each of the parameters b1 and b2 are assumed to be the same across groups
             common_linear  : b3=0 and b2=0. The value of the slope b1 is assumed to be the same across groups
             dropped    	: A model without time trend is computed
group   : the different groups (e.g., race/ethnicity groups) variable
DF      : Non-linear model degrees of freedom. The default is set pretty high at 10000
print   : Yes will print the nlmixed results and No will not. Default is No
out     : The name of the output baseline. All the following outputs (baseline + suffix) are saved. 
          Here are the suffixes:
          	_fitstat : model fit estimates from the proc nlmixed
          	_ests    : model fit estimates formated for use in the estimation of the Kalman prediction
          	_covmat  : model fit covariance matrix
          	_pred    : Kalman prediction of the outcome of interest includes original values as well as parameters
         (e.g., for OUT=result then RESULT_PRED will be the Kalman prediction data of the outcome of interest. )
*/

%macro htrp(
             data=, 
             outcome=, 
             se=,
             group=,
             time=, 
             by=, 
             xtrakeep= ,
			 orpoly=YES,
             _rho_= , 
             _tausq_= , 
             bvalue= indep_linear , 
             DF=10000, 
             out=param , 
             print=NO
            );

%local n g p k nrep n1 n2 i iis jj lj b1line b2line b3line tline formatted dsop dscl _rtimess rtm rlag
       group_rep vmatinv vmat xmat xmat2 xmat3 xmat4 amat blmat _emkfkeep_ _emkfmu_ parList colList oPPmat;

/* eMKF: Data assumed to have been pre-formatted using macro reformat: check and reformat if not */

%let formatted = 0;

%let dsop = %sysfunc(open(&data));
%if &dsop ne 0 %then %do;
	%if %sysfunc(varnum(&dsop, inputorder)) ne 0 and %sysfunc(varnum(&dsop, &time)) ne 0 %then %let formatted = 1;
%end; 
%let dscl = %sysfunc(close(&dsop));

%let formatted = %eval(&formatted + 0);

data _sdata_;
run;

%if &formatted = 1 %then %do;
	data _sdata_;
	  set &data;
	run;
%end;
%else %do;

    %put ;
	%put Reformatting data prior to MLE-based estimation;

	%reformat(data=&data, 
		      outcome=&outcome, se=&se, outcome2=, se2=, 
			  group=&group, time=&time, by=&by, 
 			  outformat= _sdata_ );

	/*eMKF: Create copies of _group_ and _rep variables for use in proc iml matrix calculations */
	data _sdata_; 
	  set _sdata_;
	  _groupnum = _group_; 					 
	  %if &by ^= %str() %then _reps = _rep;;
	run;

	/*eMKF: Replace &group and &by macro variables by their numeric versions */
	%let group = _groupnum; 				
	%if &by ^= %str() %then %let by = _reps;;

%end;

/*eMKF: Sort by replications (if any), group, and time */
/*eMKF: Use _time as &time could be empty if data was still in format 1 when htrp was called */
proc sort data= _sdata_;
  by _rep _group_ _time ;
run;

/*eMKF: Macro variable for the number of groups */
%let g=;
data _freqg_;
run;
proc freq data=_sdata_ noprint;
  tables &group /list out=_freqg_;
run;
data _freqg_;
  set _freqg_;
  _group_ +1;
  call symput('g',_group_);
  keep _group_ &group;
run;

/*eMKF: Macro variable for the number of time points */
%let n=;
data _freqn_;
run;
proc freq data=_sdata_ noprint;
  tables _rtime /list out=_freqn_;
run;
data _freqn_;
  set _freqn_;
  _time +1;
  call symput('n',_time);
  keep _time _rtime;
run;

/*eMKF: Macro variable for the real times to use in calculations */
%let _rtimess = ;
data _freqn_;
  set _freqn_;
  retain _rts;
  if _n_= 1 then _rts = cat(_rtime);
  else _rts = catx(" ", _rts, _rtime);
  call symput('_rtimess', _rts);
  drop _rts;
run;

/* eMKF allows for irregular and fractional times points */
/* This is the variable that will be used for real time in case it is not just 1,2,3,... */
%let rtm=_rtime;

/*eMKF: macro reformat now also tracks lags between successive real time points */		
%let rlag=_rlag;

/*eMKF: Macro variable for the number of replications */
%let nrep=1;
%if &by ^=%str() %then %do;
	data _freq_;
	run;
	proc freq data=_sdata_ noprint;
	  tables &by /list out=_freq_;
	  format &by ;
	run;
	data _freq_;
	  set _freq_;
	  _rep +1;
	  call symput('nrep',_rep);
	  keep _rep &by;
	run;
%end;

/*eMKF: Set numerical values to use in code */
%let n=%eval(0+&n);
%let g=%eval(0+&g);
%let nrep=%eval(0+&nrep);

/* eMKF: Modification to set up orthogonal cubic polynomial design matrix */

data _oXmat_ _oPmat_;
run;

%if %upcase(&orpoly) = YES %then %do;

	proc iml;
	  x = { &_rtimess };	
	  x = T(x);							/* eMKF: column vector with real times */
	  oP = orpol(x, 3);					/* eMKF: orthonormal design matrix for cubic orthogonal polynomials */
	  x0 = { %cnstss(1, &n) };
	  x0 = T(x0);
	  x1 = x;
	  x2 = x#x;
	  x3 = x#x2;
	  uP = x0 || x1 || x2 || x3;		/* eMKF: raw/unstandardized design matrix */
	  oP1 = inv(T(uP)*uP)*T(uP)*oP[,1];
      oP2 = inv(T(uP)*uP)*T(uP)*oP[,2];
      oP3 = inv(T(uP)*uP)*T(uP)*oP[,3];
      oP4 = inv(T(uP)*uP)*T(uP)*oP[,4];
	  oPP = oP1 || oP2 || oP3 || oP4;	/* eMKF: right multiplication of raw uP with oPP produces orthonormal oP */
	  y = T(do(1, &n, 1));				/* eMKF: column vector of consecutive time indices */
	  yP = y || oP;
	  create _oXmat_ from yP [ colname = {"_time" "&rtm.0" "&rtm.1" "&rtm.2" "&rtm.3"} ] ;
	  append from yP; close _oXmat_;
	  create _oPmat_ from oPP [ colname = {"t0" "t1" "t2" "t3"} ] ;
	  append from oPP; close _oPmat_;
	quit;

	proc sort data=_sdata_;
	  by _time;
	run;

	data _sdata_;
	   merge _sdata_ _oXmat_;
	   by _time;
	run;

	proc sort data= _sdata_;
	  by _rep _group_ _time ;
	run;

%end;
%else %do;
	data _sdata_;
	   set _sdata_;
	   &rtm.0 = 1;
	   &rtm.1 = &rtm;
	   &rtm.2 = &rtm**2;				/* eMKF: add raw quad and cubic time terms as columns in _sdata_ dataset */
	   &rtm.3 = &rtm**3;
	run;
%end;

/* eMKF: Set up model statement symbolically for use in proc reg
 * Initial regression when no time trend is desired differs from the original MKF, here.  
 * In MKF, initial values for intercepts were based on those from linear trend model instead of dropped model.
 */
%let tline=;
%if %upcase(&bvalue) = COMMON_CUBIC  or %upcase(&bvalue) = INDEP_CUBIC  %then %let tline = &rtm.0 &rtm.1 &rtm.2 &rtm.3;
%if %upcase(&bvalue) = COMMON_QUAD   or %upcase(&bvalue) = INDEP_QUAD   %then %let tline = &rtm.0 &rtm.1 &rtm.2;
%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = INDEP_LINEAR %then %let tline = &rtm.0 &rtm.1;
%if %upcase(&bvalue) = DROPPED 											%then %let tline = &rtm.0;

/* eMKF: Modified call to proc reg to include quad and cubic terms
 * Implemented regression by _rep instead of relying only on _rep = 1 for initial values to pass to proc nlmixed
 */
proc reg data=_sdata_ outest=_beta_ noprint;
   by _rep _group_;
   model _y = &tline /noint;;
run;

/* Setup initial values (modified for eMKF to hold quad and cubic coefficients)
   Use the beta estimates from the model above. For &g groups, then:
     a1-a&g are the intercept from each of the &g groups models
     b1arr1-b1arr&g are the linear coefficients from each of the &g groups models
     b2arr1-b2arr&g are the quadratic coefficients from each of the &g groups models
     b3arr1-b3arr&g are the cubic coefficients from each of the &g groups models
*/

/*eMKF: Sort regression coefficients by _rep and _group_ */
proc sort data=_beta_;
  by _rep _group_;
run;

/*eMKF: Initial parameter values for the first replication */
%let jj = 1;
data _inits_;
	set _beta_(where=(_rep = &jj)) end=end;
	_rep = &jj;
   	array a a1-a&g;
   	array b1arr b1arr1-b1arr&g;
   	array b2arr b2arr1-b2arr&g;
   	array b3arr b3arr1-b3arr&g;
   	retain a1-a&g b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g; 
	a{_group_} = &rtm.0;
   	%if %upcase(&bvalue) ^= DROPPED %then %str( b1arr{_group_} = &rtm.1;) ;;
   	%if %upcase(&bvalue) ^= DROPPED and %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= INDEP_LINEAR %then %str( b2arr{_group_} = &rtm.2;) ;;
   	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = INDEP_CUBIC %then %str( b3arr{_group_} = &rtm.3;) ;;
   	if end then do;
   		%if &_rho_ = %str() %then logitrho = 0;;
		%if &_tausq_ = %str() %then logtau2 = log(.002);;
      	output;
   	end;
   	keep _rep 
		a1-a&g b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g
		%if &_rho_ = %str() %then logitrho;
		%if &_tausq_ = %str() %then logtau2;
    ;
run;

/*eMKF: Initial estimates from each subsequent replication */
data _initsr_;
run;
%if &nrep > 1 %then %do; 	
	%do jj = 2 %to &nrep;
		data _initsr_;
			set _beta_(where=(_rep = &jj)) end=end;
			_rep = &jj;
		   	array a a1-a&g;
		   	array b1arr b1arr1-b1arr&g;
		   	array b2arr b2arr1-b2arr&g;
		   	array b3arr b3arr1-b3arr&g;
		   	retain a1-a&g b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g; 
			a{_group_} = &rtm.0;
		   	%if %upcase(&bvalue) ^= DROPPED %then %str( b1arr{_group_} = &rtm.1;) ;;
		   	%if %upcase(&bvalue) ^= DROPPED and %upcase(&bvalue) ^= COMMON_LINEAR and 
				%upcase(&bvalue) ^= INDEP_LINEAR %then %str( b2arr{_group_} = &rtm.2;) ;;
		   	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = INDEP_CUBIC %then %str( b3arr{_group_} = &rtm.3;) ;;
		   	if end then do;
		   		%if &_rho_ = %str() %then logitrho = 0;;
				%if &_tausq_ = %str() %then logtau2 = log(.002);;
		      	output;
		   	end;
		   	keep _rep
				a1-a&g b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g
				%if &_rho_ = %str() %then logitrho;
				%if &_tausq_ = %str() %then logtau2;
		    ;
		run;
		data _inits_;
		  set _inits_ _initsr_;
		run;
		data _initsr_;
		run;
	%end;
%end;
  
/* eMKF: calculate initial estimates for common_ scenarios and remove extraneous variables */
data _inits_; 
	set _inits_;
	b1=.;
	b2=.;
	b3=.;
    %if %upcase(&bvalue)=COMMON_CUBIC or %upcase(&bvalue)=COMMON_QUAD or %upcase(&bvalue)=COMMON_LINEAR %then %str(b1 = mean(of b1arr1-b1arr&g);) ;;
    %if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD and %upcase(&bvalue) ^= INDEP_LINEAR %then %str(drop b1arr1-b1arr&g;) ;;
    %if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD %then %str(b2 = mean(of b2arr1-b2arr&g);) ;;
	%if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD %then %str(drop b2arr1-b2arr&g;) ;;
 	%if %upcase(&bvalue) = COMMON_CUBIC %then %str(b3 = mean(of b3arr1-b3arr&g);) ;;
 	%if %upcase(&bvalue) ^= INDEP_CUBIC %then %str(drop b3arr1-b3arr&g;) ;;
	%if %upcase(&bvalue) ^= COMMON_CUBIC and %upcase(&bvalue) ^= COMMON_QUAD and %upcase(&bvalue) ^= COMMON_LINEAR %then drop b1;;;
    %if %upcase(&bvalue) ^= COMMON_CUBIC and %upcase(&bvalue) ^= COMMON_QUAD %then drop b2;;;
	%if %upcase(&bvalue) ^= COMMON_CUBIC %then drop b3;;;
run;

/* eMKF: re-order columns in common trend cases so logitrho and logtau2 remain last */
%if %upcase(&bvalue)=COMMON_CUBIC or %upcase(&bvalue)=COMMON_QUAD or %upcase(&bvalue)=COMMON_LINEAR %then %do;
	data _inits_;
	  merge _inits_(drop = logitrho logtau2) _inits_(keep = _rep logitrho logtau2);
	  by _rep;
	run;
%end;

/*eMKF: Set up array declarations symbolically for use in proc nlmixed */
%let b1line=; %let b2line=; %let b3line=;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %let b1line= array b1arr(&g) b1arr1-b1arr&g ;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then %let b2line= array b2arr(&g) b2arr1-b2arr&g ;
%if %upcase(&bvalue) = INDEP_CUBIC %then %let b3line= array b3arr(&g) b3arr1-b3arr&g ;

/*eMKF: Set up model statement symbolically for use in proc nlmixed */
%let tline=;
%if %upcase(&bvalue) = INDEP_CUBIC   %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1arr[_group_]+&rtm.2*b2arr[_group_]+&rtm.3*b3arr[_group_]+gamma[_time],_se**2);
%if %upcase(&bvalue) = INDEP_QUAD    %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1arr[_group_]+&rtm.2*b2arr[_group_]+gamma[_time],_se**2);
%if %upcase(&bvalue) = INDEP_LINEAR  %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1arr[_group_]+gamma[_time],_se**2);
%if %upcase(&bvalue) = COMMON_CUBIC  %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1+&rtm.2*b2+&rtm.3*b3+gamma[_time],_se**2);
%if %upcase(&bvalue) = COMMON_QUAD   %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1+&rtm.2*b2+gamma[_time],_se**2);
%if %upcase(&bvalue) = COMMON_LINEAR %then %let tline= normal(&rtm.0*a[_group_]+&rtm.1*b1+gamma[_time],_se**2);
%if %upcase(&bvalue) = DROPPED       %then %let tline= normal(&rtm.0*a[_group_]+gamma[_time],_se**2);

/* Fit a Non-linear mixed model. Capture covariance matrix COV and covariance matrix from additional estimate ECOV */

data _ests _fitstat &out._fitstat &out._covmat &out._data;
run;

data &out._data;
  set _sdata_;
  drop &rtm.0 &rtm.1 &rtm.2 &rtm.3; /* eMKF: remove polynomial time terms from &out._data */
run;

%put Start model fitting using PROC NLMIXED; /* eMKF: substituted put statement */
%put Model is _y ~ &tline;

/*eMKF: model print handling */
%if %upcase(&print) ^= YES %then ods exclude all;; 

/*eMKF: initial counter for use in random statement */
%let i=0;

proc nlmixed data=_sdata_ DF=&DF cov ecov
	method=firo maxiter=500; 		 /* eMKF: updated method and maxiter to match those in htrp2d */
	by _rep; 						 /* eMKF: stratified by _rep (trivial if _rep is constant = 1 (i.e., no stratification) */
	array gamma(&n) gamma1-gamma&n ; /* Give the gamma values a generic name */
	array a(&g) a1-a&g; 			 /* Give the a values a generic name */
	&b1line;; 						 /* Give the bl values a generic name */
	&b2line;; 						 /* eMKF: Give the b2 values a generic name */
	&b3line;; 						 /* eMKF: Give the b3 values a generic name */ 
	parms / bydata data=_inits_; 	 /* eMKF: Setup parameters and their initial values stratified by _rep */

	/* 
	Next define tausq to be an assigned value of tausq but in the case &_tausq_ is missing 
	 then assign value exp(logtau2) a value that will be estimated.
	Same thing for _rho_ taking assigned value but in case &_rho_is missing, 
	 give it the value 2/(1+exp(-logitrho)) - 1
	*/
	%if &_tausq_ = %str() %then _tausq_ = exp(logtau2); ;
	%if &_tausq_ ^= %str() %then _tausq_ = &_tausq_; ;
	%if &_rho_ = %str() %then _rho_ = 2/(1+exp(-logitrho)) - 1;;
	%if &_rho_ ^= %str() %then _rho_ = &_rho_; ;
	nu = _tausq_/(1-(_rho_**2));

	/* Random effects gamma1 gamma2 ... gamma&n ~ normal([0,0,..,0],[nu, half triangular matrix]) */
	/*eMFK: Replaced do loop for variances with call to macro thevarcompr to allow for noninteger/unequally-spaced time points */
	random %do i = 1 %to &n; gamma&i %end;
	~ normal([%zeros(&n)], %thevarcompr(times=&_rtimess, nu=nu, vrho=_rho_)) subject=_group_;

	model _y ~ &tline;;

	%if &_rho_ = %str() %then estimate "_rho_" 2/(1+exp(-logitrho))-1;;
	%if &_tausq_ = %str() %then estimate "_tausq_" exp(logtau2);;

	ods output parameterestimates=_ests fitstatistics=_fitstat;
	ods output CovMatParmEst=&out._covmat;
run;

/*eMKF: Replace missing values in CovMatParmEst with 0s for later use */
%if %upcase(&orpoly) = YES %then %do;
	data &out._covmat;
	  set &out._covmat;
	  array NAs _numeric_;
	  do over NAs;
	      if NAs = . then NAs = 0;
	  end;
	run;
%end;

/*eMKF: Reset ODS destinations if they were turned off */
%if %upcase(&print) ^= YES %then ods exclude none;;

/* eMKF: Added put statements for log */
%put End model fitting using PROC NLMIXED;
%if %upcase(&print) ^= YES %then %put Model printout was turned off by user;;

/*eMKF: Model estimates */

data _fitstat;
   set _fitstat;
   if descr = "-2 Log Likelihood";
   _2loglike=value;
   keep _2loglike _rep;
run;

data &out._fitstat;
  merge _ests _fitstat;
  by _rep;
run;

data _llike_;
run;
data _llike_;
  set &out._fitstat;
  by _rep;
  if first._rep then output;
  keep _rep _2loglike;
run;

data &out._ests;
run;
proc transpose data=&out._fitstat out=&out._ests;
   by _rep;
   var estimate;
   id parameter;
run;

/* eMKF: Macro variables _emkfkeep_ and _emkfmu_ account for various model combinations */
%let _emkfkeep_ = ; 
%let _emkfmu_ = ;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC %then %do;
	%let _emkfkeep_ = a b1 b2 b3;
	%let _emkfmu_ = mu=a*&rtm.0+b1*&rtm.1+b2*&rtm.2+b3*&rtm.3;
%end;
%if %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD %then %do;
	%let _emkfkeep_ = a b1 b2;
	%let _emkfmu_ = mu=a*&rtm.0+b1*&rtm.1+b2*&rtm.2;
%end;
%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = COMMON_LINEAR %then %do;
	%let _emkfkeep_ = a b1;
	%let _emkfmu_ = mu=a*&rtm.0+b1*&rtm.1;
%end;
%if %upcase(&bvalue) = DROPPED %then %do;
	%let _emkfkeep_ = a; 
	%let _emkfmu_ = mu=a*&rtm.0;
%end;

/*eMKF: Re-structure estimates by _rep and _group_ */
data &out._ests;
   merge &out._ests _llike_;
   by _rep;
   array as a1-a&g;
   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %str(array b1s b1arr1-b1arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(b1 = 0;) ;;
   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then %str(array b2s b2arr1-b2arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(b2 = 0;) ;;
   %if %upcase(&bvalue) = INDEP_CUBIC %then %str(array b3s b3arr1-b3arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(b3 = 0;) ;;
   _rho_ = &_rho_ %if &_rho_ = %str() %then 2/(1+exp(-logitrho))-1  ;;;
   _tausq_ = &_tausq_ %if &_tausq_ = %str() %then exp(logtau2) ;;; 
   do _group_ = 1 to &g;
      a = as{_group_};
      %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %str(b1 = b1s{_group_};) ;;
	  %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then %str(b2 = b2s{_group_};) ;;
      %if %upcase(&bvalue) = INDEP_CUBIC %then %str(b3 = b3s{_group_};) ;;
      output;
   end;
   keep _rep _group_ &_emkfkeep_ _rho_ _tausq_ _2loglike;;
run;

/*eMKF: Merge with group labels */
proc sort data= &out._ests;
  by _group_ ;
run;
proc sort data=_freqg_;
  by _group_;
run;
data &out._ests;
  merge &out._ests _freqg_;
  by _group_;
run;

/*eMKF: Re-order columns so that intercept is always listed first */
data &out._ests; 
  merge &out._ests(drop= &_emkfkeep_) &out._ests(keep= a) %if %length(&_emkfkeep_) > 1 %then &out._ests(keep= %substr(&_emkfkeep_, 3));;; 
run;

/***************************************************************************/
/* Now let's use the Kalman technique to estimate the outcome observations */
/***************************************************************************/

data &out._pred;
run;
proc sort data=_sdata_ out=&out._pred;
  by _rep _group_ _time; 
run;

proc sort data=&out._ests;
  by _rep _group_;
run;

data &out._pred; 
  merge &out._pred &out._ests;
  by _rep _group_ ;
  &_emkfmu_;; 					/*eMKF: invoke symbolic calculation for mu */
  err=_y - mu;
run;

data _empty_;
  set &out._pred;
  if _time= 1;
  w0 = _tausq_/(1-(_rho_**2));
  gamma0 = 0 ;
  _time= 0;
  &rtm= 0;
  &rlag = .;
  keep _time &rtm &rtm.0 &rtm.1 &rtm.2 &rtm.3 &rlag &by _rep &group _group_ w0 gamma0;
run;

data &out._pred;
 set &out._pred _empty_;
 if _time = 1 then &rlag = 1;	/*eMKF: set to lag 1 instead of missing relative to time 0 */
run;

proc sort data=&out._pred;
 by _rep _group_ _time;
run;

/*eMKF: recursion formulas modified to allow lag s > 0 between time points */
data &out._pred;
  set &out._pred;
  retain wold gamold;
  delta = ((_rho_**(2*&rlag)) * wold) + _tausq_*(1 - (_rho_**(2*&rlag)))/(1 - (_rho_**2));
  if &rlag > 0 or &rlag = . then lambda = delta/(delta + (_se**2));
  else lambda = 0; 				/*eMKF: in the limiting case, recursion breaks down: set lambda = 0 instead */
  w = delta * (1-lambda);
  gamma = (lambda* err) + (1-lambda)*(_rho_**&rlag)*gamold ; /*eMKF: &rlag=1 reduces to unit-increment case */
  prediction = mu + gamma ;
  output;
  if _time = 0 then do;
	wold=w0;
	gamold=gamma0;
  end;  
  if _time ^= 0 then do; 
	wold=w;
    gamold=gamma;
  end;
run;

%let group_rep = compress(_rep || _group_);

data &out._pred;
  set &out._pred;
  group_rep= &group_rep;
  if _time ne 0;
  drop w0 gamma0 wold gamold;
run;

/*eMKF: reset lag for time 1 to missing */
data &out._pred;
  set &out._pred;
  if _time = 1 then &rlag = .; 
run;

/******************/
/* Compute the MSE*/
/******************/

data _Amat_ _Dmat_ _Vmat_ _Vgmat_ _Vemat_  _Xmat_ _junk_ _junk0_ _junk01_;
run;

proc sort data=&out._pred(keep= _rep &by _group_ &group _time &rtm 
								&rtm.0 &rtm.1 &rtm.2 &rtm.3 		/*eMKF: modification to deal with orthogonal polynomials */
								&rlag _y _se _rho_ _tausq_ err lambda prediction group_rep)
  out=_Amat_ ;
  by group_rep;
run;

data _Amat_;
  set _Amat_; 
  by group_rep;
  array hs ah1-ah&n (%zeross(&n));
  if first.group_rep then do _k=1 to &n; hs{_k}=0; end;
  if _time = 1 then ah1=lambda;
  else do;
  	 /*eMKF: updated to account for &rlag if not 1 */
     do _j=1 to _time - 1;
	    hs{_j} = hs{_j}*(1-lambda)*(_rho_**&rlag) ; 
	 end;
	 hs{_time}=lambda;
  end;
  
  drop group_rep _k _j;
  keep _rep &by _group_ &group _time &rtm &rlag ah:;
run;

proc sort data=_Amat_;
  by _rep _group_ _time  ;
run;

/* Attach row number _id to each of the A Matrix so that when computing group by group
   the appropriate Ag could be called*/
data _Amat_;
  set _Amat_;
  _id+1;
run;

data _Dmat_;
  set &out._pred(keep= _rep &by _group_ &group _time &rtm 
						&rtm.0 &rtm.1 &rtm.2 &rtm.3  /*eMKF: modification to deal with orthogonal polynomials */
						&rlag _y _rho_ _tausq_ _se);
  w0 = _tausq_/(1- (_rho_**2));
run;

/* The Matrix _Dmat_ is actually the same within group and within replication*/
proc sort data= _Dmat_ nodupkey;
  by _rep _group_ _time;
run;

/* Here the standard error is the same from one replication to another
 by different from group to group. So these matrices should just be estimated 
 differently for each group as well as replication if possible*/

data _Vmat_; 
  set _Dmat_;
  array ad ad1-ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n; 
  	 if i = _time  then ad{i} = w0 + _se**2; /* Add the variances of Y to the diagonal elements */
   	 else ad{i} = w0*(_rho_**abs(rt{i} - &rtm));
  end;
  drop i rt1-rt&n;
  keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n;
run;

/* Recreate the variance where only the Variance of gamma is estimated
 This is mostly needed for check of what the estimates are giving us
*/
data _Vgmat_;
  set _Dmat_;
  array ad ad1-ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n;
    if i = _time then ad{i} = w0;
    else ad{i} = w0*(_rho_**abs(rt{i} - &rtm));
  end;
  drop i rt1-rt&n;
  keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n;
run;

/* This is the diagonal Variance matrix of the errors */
data _Vemat_;
  set _Dmat_;
  array ad ad1-ad&n;
  do i=1 to &n;
    if i = _time then ad{i}= _se**2;
    else ad{i} = 0;
  end;
  drop i;
  keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n;
run;

/* This is just like the _Vgmat_ but with extra variables kept to be used later
 for the creation of the X matrix for example*/
data _Dmat_;
  set _Dmat_;
  array ad ad1-ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n;
    if i = _time then ad{i} = w0;
    else ad{i} = w0*(_rho_**abs(rt{i} - &rtm));
  end;
  drop i rt1-rt&n;
run;

/* This is the X matrix */
data _Xmat_;
  set _Dmat_;
  /*eMKF: modification to deal with orthogonal polynomials */
  x0 = &rtm.0; 
  x1 = &rtm.1;
  x2 = &rtm.2;
  x3 = &rtm.3;
  keep _rep &by _group_ &group x0 x1 x2 x3;
run;

data &out._H &out._PredVar &out._CovY;
run;

%let iis=;
%let vmatinv= ;
%let vmat=;
%let xmat= ;
%let xmat2= ;
%let xmat3= ; /* eMKF: xmat3 and xmat4 added to deal with quad and cubic terms */
%let xmat4= ;
%let amat= ;
%let blmat= ;

/* Now let's capture the A Matrix and turn it into a diagonal block matrix  */

data _junk_;
  set _Amat_;
  by _rep;
  if first._rep then sid=0;
  sid+1;
  repgrp=compress(_rep ||"-"|| _group_);
run;

proc sort data=_junk_;
  by repgrp sid;
run;

data _junk_;
  set _junk_;
  by repgrp;
  if first.repgrp then kp=1;
  if last.repgrp then kp=2;
  if kp=1 or kp=2;
  keep _rep &by _group_ &group _time &rtm &rlag sid kp repgrp; 	/*eMKF: also keeping &rtm and &rlag */
run;

data _junk_;
  set _junk_;
  by repgrp;
  retain minid maxid;
  array Aid(1:2) minid maxid;
  if first.repgrp then do;
    do i = 1 to 2;
      Aid[i] = .; /*initializing to missing*/
    end;
  end;
  Aid(kp) = sid; 
  if last.repgrp then output; 
  drop kp sid i;
run;

proc sort data=_junk_;
 by _rep _group_ _time;
run;

%let jj = 0; %let lj = 0;
data _junk_;
  set _junk_;
  by _rep;
  mm=0;
  if last._rep then mm=1;
  amat=compress("Z["||minid||":"||maxid||",]"); 			/* Here Z will be the standard matrix definition that can be used */
  if mm=0 then amat=compress(amat||",");
  if first._rep or mm=0 then amat=compress("block("||amat); /* eMKF: modified to allow for arbitrary number of groups */
  call symput("jj", _n_);
  if last._rep then call symput("lj", length(amat));		/* eMKF: added to capture maximal length of character variable needed */
run;

/* eMKF: modification to allow for arbitrary number of groups (function block restricted to 15) */
data _junk_;
  set _junk_;
  if _n_=&jj then mm=2;
  length amat2 $ %eval(&g+&lj+10);							/* eMKF: set large enough character length to cover block( + lj + g times ) */
  amat2 = amat;
  %if &g = 1 %then %do;
  	  if mm=1 then amat2=compress(amat||") ");
  	  if mm=2 then amat2=compress(amat||") ");
  %end;
  %if &g > 1 %then %do;
      if mm=1 then amat2=compress(amat||repeat(") ", &g-2));
      if mm=2 then amat2=compress(amat||repeat(") ", &g-2));
  %end;
  drop amat;
  rename amat2 = amat;
run;

proc sort data=_junk_;
  by _rep _group_;
run;

data _junk_;
  set _junk_;
  name=compress("&group"||_group_);
run;

data _junk0_ _junk01_;
run;

proc transpose data=_junk_(keep=name _rep &by amat) out=_junk0_;
  var amat;
  by _rep;
  id name;
run;

proc sort data=_junk_ out=_junk01_ nodupkey;
  by _group_;
run;

%let amat=;
proc sql noprint;
  select name into :amat separated by ' || ' from _junk01_;
quit;

data _junk0_;
  set _junk0_;
  i+1;
  code=compress("Z"||i);
  mcode=compress("_Z"||i);
  amat=compress(code||"="||&amat);
  keep amat _rep code mcode;
run;

%let jj=0;
%do jj=1 %to &nrep;
  %local _Z&jj ;
%end;

data _null_;
  set _junk0_;
  call symput(mcode, amat);
run;

/* eMKF: call to function block in RAND's MKF macro limited the number of blocks (groups) to 15, which produced an error if &g > 15
 * eMKF: code was revised to allow for additional blocks (groups) using nested calls to function block
 */
%let jj=0;
%do jj=1 %to &g;

	%if &jj = 1 %then %let iis = i_&n ;
	%if &jj > 1 %then %let iis = &iis // i_&n ;

	/* eMKF: Setup an (ng)x(ng) matrix to capture the full V for all groups */
	%if &jj = 1 %then %let vmat = block( V ;
	%if &jj > 1 and &jj < &g %then %let vmat = &vmat , block ( V ;
	%if &jj = &g and &g > 1  %then %let vmat = &vmat , V %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let vmat = &vmat );

	/* eMKF: V inverse */
	%if &jj = 1 %then %let vmatinv = block( invV ;
	%if &jj > 1 and &jj < &g %then %let vmatinv = &vmatinv , block ( invV ;
	%if &jj = &g and &g > 1  %then %let vmatinv = &vmatinv , invV %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let vmatinv = &vmatinv );

	/* Standard way to create the block matrices variable is B*/
	%if &jj = 1 %then %let blmat = block( B ;
	%if &jj > 1 and &jj < &g %then %let blmat = &blmat , block ( B ;
	%if &jj = &g and &g > 1  %then %let blmat = &blmat , B %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let blmat = &blmat );

	/* Look into the value of X that will be used in the different scenarios*/
	/* eMKF: modified to account for quadratic and cubic scenarios */
	%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %do; 
		/* eMKF: general form for independent trends */
		%if &jj = 1 %then %let xmat = block( X ;
		%if &jj > 1 and &jj < &g %then %let xmat = &xmat , block ( X ;
		%if &jj = &g and &g > 1  %then %let xmat = &xmat , X %sysfunc(repeat( %str(%)), &g-2));
		%if &jj = &g and &g = 1  %then %let xmat = &xmat );
  	%end;
  	%if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD  and %upcase(&bvalue) ^=INDEP_LINEAR %then %do; 
		/*eMKF: intercept always included */
		%if &jj = 1 %then %let xmat = block( X0 ;
		%if &jj > 1 and &jj < &g %then %let xmat = &xmat , block ( X0 ;
		%if &jj = &g and &g > 1  %then %let xmat = &xmat , X0 %sysfunc(repeat( %str(%)), &g-2));
		%if &jj = &g and &g = 1  %then %let xmat = &xmat );
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_LINEAR %then %do; 
		/*eMKF: linear time */
		%if &jj = 1 %then %let xmat2 = X1 ;
		%if &jj > 1 %then %let xmat2 = &xmat2 // X1 ;
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD %then %do; 
		/*eMKF: quadratic time */
		%if &jj = 1 %then %let xmat3 = X2 ;
		%if &jj > 1 %then %let xmat3 = &xmat3 // X2 ;
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC %then %do; 
		/*eMKF: cubic time */
		%if &jj = 1 %then %let xmat4 = X3 ;
		%if &jj > 1 %then %let xmat4 = &xmat4 // X3 ;
  	%end;

%end;  /* End of &jj */

/* eMKF: modified to account for quadratic and cubic scenarios */
%if %upcase(&bvalue) = COMMON_LINEAR %then  %let xmat = &xmat || ( &xmat2 );
%if %upcase(&bvalue) = COMMON_QUAD   %then  %let xmat = &xmat || ( &xmat2 ) || ( &xmat3 );
%if %upcase(&bvalue) = COMMON_CUBIC  %then  %let xmat = &xmat || ( &xmat2 ) || ( &xmat3 ) || ( &xmat4 );

%let jj=0;
%let n1=2;
%let n2= %eval(&n + 1);

/* eMKF: added dimensionality p for ease of coding of quadratic and cubic trend models in proc iml */
%let p=0;
%if %upcase(&bvalue) = INDEP_CUBIC  or %upcase(&bvalue) = COMMON_CUBIC  %then %let p = 4;
%if %upcase(&bvalue) = INDEP_QUAD   or %upcase(&bvalue) = COMMON_QUAD   %then %let p = 3;
%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = COMMON_LINEAR %then %let p = 2;
%if %upcase(&bvalue) = DROPPED 											%then %let p = 1;

proc iml;

	 %do jj=1 %to &nrep;

		/* Next for the Xmatrix */
		/* eMKF: modified to allow for quadratic and cubic trend models */
		use _Xmat_(where=(_rep=&jj) keep= _rep x:);
		read all into XX; close _Xmat_; /* eMKF: also added close statements for cleanliness */
		Z =XX[,2:(1+&p)];
		&&_Z&jj ;;
	 	Xs = Z&jj;

		/* Next for the V matrix */
		use _Vmat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VV; close _Vmat_; 
		Z =VV[,&n1:&n2];
		&&_Z&jj ;;
		Vs = Z&jj;
		invVs=inv(Vs);

		/* Next for the Ve matrix */
		use _Vemat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VVe; close _Vemat_;
		Z =VVe[,&n1:&n2];
		&&_Z&jj ;;
		Ves = Z&jj;

		/* Next for the Vg matrix */
		use _Vgmat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VVg; close _Vgmat_;
		Z =VVg[,&n1:&n2];
		&&_Z&jj ;;
		Vgs = Z&jj;

		/* Next for the A matrix */
		use _Amat_(where=(_rep=&jj)  keep= _rep ah:);
	    read all into AA; close _Amat_;
		Z =AA[,&n1:&n2];
		&&_Z&jj ;;
		As = Z&jj;

		/* Next for the data label matrix */
		use _Dmat_(where=(_rep=&jj) keep=_rep &by &group &rtm _time _y _se);  /*eMKF: also keeping &rtm */
		read all var{_y} into Y;
		read all var{_rep &by &group &rtm _time} into NM; 					/*eMKF: also keeping &rtm */
		close _Dmat_;

		/* Now do the estimations */
		i_&n = i(&n*&g);
		H = Xs * inv(t(Xs)*invVs*Xs)*t(Xs)*invVs;
		fH = H + As*(i_&n - H);
		fY = fH * Y;
		Vy = vecdiag(fH * Vs * t(fH));
		MSEy = vecdiag(  (fH - i_&n) * Vgs * t(fH - i_&n)   ) + vecdiag(fH * Ves * t(fH));
		ff= NM || fH;
		fV= NM || Vs;
		fVy=NM || fY || Vy || MSEy;

		%if &jj = 1 %then ffs= ff;;
		%if &jj = 1 %then fVs= fV;;
		%if &jj = 1 %then fVys= fVy;;
		%if &jj > 1 %then ffs= ffs // ff;;
		%if &jj > 1 %then fVs= fVs // fV;;
		%if &jj > 1 %then fVys= fVys // fVy;;

	 %end;

	 create &out._H from ffs ;
	 append from ffs; close &out._H;

	 create &out._CovY from fVs ;
	 append from fVs; close &out._CovY;
	 
	 %if &by ^=%str() %then create &out._PredVar from fVys [ colname = {"_rep" "&by" "&group" "&rtm" "_time" "Hat_y" "PredVar" "HatMSE"} ];;;
	 %if &by  =%str() %then create &out._PredVar from fVys [ colname = {"_rep" "&group" "&rtm" "_time" "Hat_y" "PredVar" "HatMSE"} ];;;
	 append from fVys; close &out._PredVar;
	 
quit; /*eMKF: ends call to proc iml with matrix calculations*/

/*eMKF: column names */
data &out._H;
  set &out._H;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4=&rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep 	       col2= &group col3=&rtm col4=_time;;;
run;

/*eMKF: column names */
data &out._CovY;
  set &out._CovY;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4=&rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep          col2= &group col3=&rtm col4=_time;;;
run;

/*eMKF: merge with predictions dataset */
data &out._pred;
  merge &out._pred &out._PredVar;
  by _rep &group _time;
  drop &rtm.0 &rtm.1 &rtm.2 &rtm.3; /* eMKF: remove polynomial time terms from &out._pred */
run;

/************************/
/* End of MSE estimation*/
/************************/

/*************************************************************************/
/* eMKF: Reverse-transform regression coefficients and covariance matrix */
/*************************************************************************/

data _tests _lests _tfits _tcmat _covmat _covmatt _tcovmat _tcovmat2 _tcovmatt _tcovmatt2 _dcovmatt ;
run;

%if %upcase(&orpoly) = YES %then %do;

	/* eMKF: reverse-transform regression coefficients */
	%let k = %eval(&p - 1); %let jj = 0; 
	proc iml;
		use _oPmat_;
		read all into oPP; close _oPmat_;
		oPP = oPP[1:&p, 1:&p];
		varNames = {"_rep" "_group_" "a"};
		%if &p > 1 %then %do;
			bNames = "b1":"b&k";	
			varNames = varNames || bNames;
		%end;
		create _tests var varNames;
		%do jj=1 %to &nrep;
			use &out._ests(where=(_rep = &jj) keep= _rep _group_ a %if &p > 1 %then b: ; ) ;
			read all into oB; close &out._ests;
			oB1 = oB[,1:2];
            oB = T(oB[,3:(2+&p)]);
            oBB = oPP * oB;
			oBB = T(oBB);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _tests;
	quit;

	/* eMKF: sort by _rep and _group_ */
	proc sort data=_tests;
  	  by _rep _group_ ;
	run;

	/* eMKF: update estimates dataset */
	data &out._ests;
  	  merge &out._ests(drop=a %if &p > 1 %then b: ;) _tests;
	  by _rep _group_;
	run;

	/* eMKF: update predictions dataset */
	data &out._pred;
  	  merge &out._pred(drop=a %if &p > 1 %then b: ; mu err delta lambda w gamma prediction group_rep Hat_y PredVar HatMSE) 
			_tests
			&out._pred(keep=_rep _group_ mu err delta lambda w gamma prediction group_rep Hat_y PredVar HatMSE) 
	  		;
	  by _rep _group_;
	run;

	/* eMKF: symbolic set up for reverse-transformation of covariance matrix */
	%let parList = ; %let colList = ; %let oPPmat = ;
	%let i = 0;
	%do i=1 %to &g;  /* eMKF: block diagonal by group */
		%if &i = 1 %then %let oPPmat = block( oP ;
		%if &i > 1 and &i < &g %then %let oPPmat = &oPPmat , block ( oP ;
		%if &i = &g and &g > 1  %then %let oPPmat = &oPPmat , oP %sysfunc(repeat( %str(%)), &g-2));
		%if &i = &g and &g = 1  %then %let oPPmat = &oPPmat );
	%end;
	%let i = 0;
	%do i=1 %to &g; /* eMKF: column and row labels */
		%let parList = &parList %bquote(")a&i%bquote(");
		%let colList = &colList a&i;
		%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList = &parList %bquote(")b1arr&i%bquote(");
			%let colList = &colList b1arr&i;
		%end;
		%if %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList = &parList %bquote(")b2arr&i%bquote(");
			%let colList = &colList b2arr&i;
		%end;
		%if %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList = &parList %bquote(")b3arr&i%bquote(");
			%let colList = &colList b3arr&i;
		%end;
	%end;
	%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList = &parList %bquote(")b1%bquote(");
		%let colList = &colList b1;
	%end;
	%if %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList = &parList %bquote(")b2%bquote(");
		%let colList = &colList b2;
	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList = &parList %bquote(")b3%bquote(");
		%let colList = &colList b3;
	%end;

	%let parList = %unquote(&parList);

	/* eMKF: square block matrix of covariances for regression coefficients by group */
	data _covmat;
	  retain _rep Row Parameter &colList;;
	  set &out._covmat(where=(Parameter in(&parList)) keep= _rep Row Parameter &colList);;
	run;

	/* eMKF: obtain new row numbers associated with modified column order */
	data _tcmat;
	  set _covmat;
      by _rep;
      if first._rep then output;
 	  drop Row Parameter;
	run;
	proc transpose data=_tcmat out=_tcmat;
      by _rep;
    run;
	data _tcmat;
	  set _tcmat;
	  nRow + 1;
	  rename _NAME_ = Parameter;
	  drop col:;
	run;

	/* eMKF: sort, merge, and re-sort using the new row numbers */
	proc sort data=_covmat;
	  by _rep Parameter;
	run;
	proc sort data=_tcmat;
	  by _rep Parameter;
	run;
	data _covmat;
	  merge _covmat _tcmat;
	  by _rep Parameter;
	run;
	proc sort data=_covmat;
	  by _rep nRow;
	run;

	/* eMKF: rectangular block matrix of covariances between regression coefficients and remaining parameters */
	data _covmatt;
	  retain _rep Row Parameter &colList;;
	  set &out._covmat(where=(Parameter not in(&parList)) keep= _rep Row Parameter &colList);;
	run;

	/* eMKF: apply reverse-transformation to both square and rectangular block matrices */
	%let i = 0; %let jj = 0;
	proc iml;

		use _oPmat_;
		read all into oP; close _oPmat_;
		oP = oP[1:&p, 1:&p];
		oPP = &oPPmat;;

		/* eMKF: re-structure block matrix in the common trend cases (where &p > 1) */
		%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
			oPP1 = oPP[do(1, &g*&p, &p), do(1, &g*&p, &p)];
			oPP2 = vecdiag(oPP[do(1, &g*&p, &p), do(2, &g*&p, &p)]);
			ToPP2 = vecdiag(oPP[do(2, &g*&p, &p), do(1, &g*&p, &p)]);
			oPP1 = oPP1 // T(ToPP2);
			oPP0 = oPP2;
			%if &p > 2 %then %do;
				oPP3 = vecdiag(oPP[do(1, &g*&p, &p), do(3, &g*&p, &p)]);
				ToPP3 = vecdiag(oPP[do(3, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP3);
				oPP0 = oPP0 || oPP3;
			%end;
			%if &p > 3 %then %do;
				oPP4 = vecdiag(oPP[do(1, &g*&p, &p), do(4, &g*&p, &p)]);
				ToPP4 = vecdiag(oPP[do(4, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP4);
				oPP0 = oPP0 || oPP4;
			%end;
			oPP0 = oPP0 // oPP[2:&p, 2:&p];
			oPP = oPP1 || oPP0;
		%end;

		varNames = {"_rep" "Row"};
		%do i=1 %to &g;
			varNames = varNames || {"a&i"};
			%if %upcase(&bvalue) = INDEP_LINEAR %then varNames = varNames || {"b1arr&i"};;
			%if %upcase(&bvalue) = INDEP_QUAD   %then varNames = varNames || {"b1arr&i"} || {"b2arr&i"};;
			%if %upcase(&bvalue) = INDEP_CUBIC  %then varNames = varNames || {"b1arr&i"} || {"b2arr&i"} || {"b3arr&i"};;
		%end;
		%if %upcase(&bvalue) = COMMON_LINEAR %then varNames = varNames || {"b1"};;
		%if %upcase(&bvalue) = COMMON_QUAD   %then varNames = varNames || {"b1" "b2"};;
		%if %upcase(&bvalue) = COMMON_CUBIC  %then varNames = varNames || {"b1" "b2" "b3"};;

		create _tcovmat var varNames;
		%do jj=1 %to &nrep;
			use _covmat(where=(_rep = &jj) drop= Parameter nRow);
			read all into oB; close _covmat;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
            oBB = oPP * oB * T(oPP);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _tcovmat;

		create _tcovmatt var varNames;
		%do jj=1 %to &nrep;
			use _covmatt(where=(_rep = &jj) drop= Parameter);
			read all into oB; close _covmatt;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
			oBB = oPP * T(oB);
			oBB = oB1 || T(oBB);
			append from oBB;
		%end;
		close _tcovmatt;

	quit;

	/* eMKF: combine both square and rectangular block matrices */
    data _tcovmat;
	  set _tcovmat _tcovmatt;
	run;

	/* eMKF: re-sort rows */
	proc sort data=_tcovmat;
  		by _rep Row ;
	run;

	/* eMKF: re-order columns as they were initially from NLMIXED */
	data _tcovmat2;
  	  retain  _rep Row a1-a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then b1arr1-b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then b1arr1-b1arr&g b2arr1-b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then b1 b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then b1 b2 b3 ; 
	  ;
	  set _tcovmat; 
	run;
	data _tcovmatt2;
  	  retain  _rep Row a1-a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then b1arr1-b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then b1arr1-b1arr&g b2arr1-b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then b1arr1-b1arr&g b2arr1-b2arr&g b3arr1-b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then b1 b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then b1 b2 b3 ; 
	  ;
	  set _tcovmatt; 
	run;
	
	/* eMKF: merge */
	data _tcovmatt2;
  	  merge &out._covmat(where=(Parameter not in(&parList)) keep = _rep Row Parameter) _tcovmatt2;
	  by _rep Row;
	run;

	/* eMKF: transpose _tcovmatt2 to add into larger matrix */
	proc transpose data=_tcovmatt2(drop=Row) out=_tcovmatt2 name = Parameter;
      by _rep;
	  id Parameter;
    run;

	/* eMKF: insert Row numbers */
	data _tcovmatt2;
	  set _tcovmatt2;
	  by _rep;
	  retain Row;
	  if first._rep then Row = 1;
	  else Row + 1;
	run;

	/* eMKF: update covariance matrix dataset */
	data &out._covmat;
  	  merge &out._covmat(keep = _rep Row Parameter) 
			_tcovmat2
 			_tcovmatt2 
			&out._covmat(where=(Parameter not in(&parList)) drop= a: %if &p > 1 %then b: ;)
			;
	  by _rep Row;
	run;

	/* eMKF: extract variances of model parameters */
	%let jj = 0; 		
	proc iml;
		create _dcovmatt var{"_rep" "Row" "Var"};
		%do jj=1 %to &nrep;
			use &out._covmat(where=(_rep = &jj) drop= Parameter);
			read all into oB; close &out._covmat;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
			oBB = vecdiag(oB);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _dcovmatt;
	quit;
	data _dcovmatt;
	  merge _dcovmatt &out._covmat(keep = _rep Row Parameter);
	  by _rep Row;
	run;

	/* eMKF: reverse-transformed regression estimates in long form */
	%if %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= COMMON_QUAD and %upcase(&bvalue) ^= COMMON_CUBIC %then %do;
		data _lests;
		  set _tests(keep = _rep a rename=(a=Est))
			  %if &p > 1 %then _tests(keep = _rep b1 rename=(b1=Est));
			  %if &p > 2 %then _tests(keep = _rep b2 rename=(b2=Est));
			  %if &p > 3 %then _tests(keep = _rep b3 rename=(b3=Est));
		  ;
		  by _rep;
		  retain Row;
		  if first._rep then Row = 1;
		  else Row + 1;
		run;
	%end;
	%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let jj = 0; 		
		data _lests;
		  set _tests(keep = _rep a rename=(a=Est))
		  	  %do jj=1 %to &nrep;
			  	  %if &p > 1 %then _tests(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep b1 rename=(b1=Est));
			   	  %if &p > 2 %then _tests(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep b2 rename=(b2=Est));
			  	  %if &p > 3 %then _tests(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep b3 rename=(b3=Est));
			  %end;
		  ;
		  by _rep;
		  retain Row;
		  if first._rep then Row = 1;
		  else Row + 1;
		run;
	%end;

	/* eMKF: merge reverse-transformed estimates and variances into fitstat dataset and update */
	data _tfits;
	  merge &out._fitstat _dcovmatt;
	  by _rep ;
	run;
	data _tfits;
	  merge _tfits _lests;
	  by _rep Row;
	run;
	data _tfits;
	  set _tfits;
	  if Est ne . then Estimate = Est;
	  if Var > 0 then StandardError = sqrt(Var);
	  else StandardError = .;
	  tValue = Estimate/StandardError;
	  Probt = (1-probt(abs(tValue), DF))*2;
	  Lower = Estimate + tinv(Alpha/2, DF)*StandardError;
	  Upper = Estimate + tinv(1-Alpha/2, DF)*StandardError;
	  drop Row Var Est;
	run;
	data &out._fitstat;
	  set _tfits;
	run;

%end;

/*eMKF: Add labels for stratification variable */
%if &by ^=%str() %then %do;
	data &out._ests;
 	  merge &out._ests _freq_;
 	  by _rep;
	run;
%end;

data &out._ests; /* eMKF: Added parameter estimates for quadratic and cubic terms */
   set &out._ests;
   label _2loglike =" -2 log-likelihood estimate"
         _rho_ ="Estimated or supplied rho of the model"
	     _tausq_="Estimated or supplied tau-square of the model"
	     _group_="Model reset group ID in case group is not ordered"
	     &group ="Numeric &&group variable"
		 /* eMKF: Added labels for by variable */
		 _rep="Model reset stratum ID in case stratification variable, if any, is not ordered"
	     &by ="&&by variable"
	     a="Parameter estimate: intercept (a)"
		 %if %upcase(&bvalue) ^=DROPPED %then b1="Parameter estimate: linear term (b1)";
   		 %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC or 
			  %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD
 				%then label b2="Parameter estimate: quadratic term (b2)" ;
   		 %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC
			   %then label b3="Parameter estimate: cubic term (b3)" ;
	;
	drop &group %if &by ^= %str() %then &by; 
    ;
run;

data &out._pred; /* eMKF: added a few useful labels */
   set &out._pred(rename=(Predvar=PredOnlyVar HatMSE=PredMSE));
   PredSE= sqrt(PredMSE);
   label _2loglike =" -2 log-likelihood estimate"
         _rho_ ="Estimated or supplied rho of the model"
	     _tausq_="Estimated or supplied tau-square of the model"
	     _group_="Model reset group ID in case group is not ordered"
	     &group ="Numeric &&group variable"
		 _rep="Model reset stratum ID in case stratification variable, if any, is not ordered"
	     &by ="&&by variable"
		 &rlag = "Elapsed real time from previous time point"
		 &rtm ="Real time used in calculations "
	     _time ="Time index variable"
	     _y  ="Original outcome"
	     _se ="Original Standard Error"
		 _avgse = "Average Standard Error across timepoints used for imputation"
		 %if &by ^= %str() %then  _avgseb = "Average Standard Error across strata used for imputation";
         impute = "Whether original Standard Error was imputed using average across timepoints"
		 %if &by ^= %str() %then imputeb = "Whether original Standard Error was imputed using average across strata";
		 inputorder = "Original ordering of the groups if it was not alphabetical"
	     prediction="Kalman estimator prediction of the outcome assuming &bvalue trend model"
	     PredMSE="Prediction variability: Mean Squared Error (MSE) assuming &bvalue trend model "
	     PredSE="Prediction standard error: Square root of MSE assuming &bvalue trend model "
		 a="Parameter estimate: intercept (a)"
	     /* eMKF: added labels for parameter estimates for quadratic and cubic terms */
   		 %if %upcase(&bvalue) ^=DROPPED %then b1="Parameter estimate: linear term (b1) assuming &bvalue trend model " ;
   		 %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC or 
			 %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD 
    			%then b2="Parameter estimate: quadratic term (b2) assuming &bvalue trend model " ;
   		 %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC 
			%then b3="Parameter estimate: cubic term (b3) assuming &bvalue trend model " ;
	;
	/* These will be deleted for now. If needed they can be useful. */
   drop mu err delta lambda w gamma Hat_y PredOnlyVar &group group_rep &rlag %if &by ^= %str() %then &by;
 ; 
run;

/*eMKF: Rename any remaining instances of numeric &group variable to _group_ + remove numeric &by variable */

data &out._predVar;
  set &out._predVar;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;

data &out._H;
  set &out._H;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;

data &out._covY;
  set &out._covY;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;

proc datasets nolist;
 delete _freqn_ _freqg_ _freq_ _sdata_ _beta_ _inits_ _initsr_ _ests _fitstat _llike_ _empty_ 
        _Amat_ _Dmat_ _Vmat_ _Vgmat_ _Vemat_  _Xmat_ _junk_ _junk0_ _junk01_ 
        _oXmat_ _oPmat_ _tests _lests _tfits _tcmat _covmat _covmatt _tcovmat _tcovmat2 _tcovmatt _tcovmatt2 _dcovmatt
        ;
run ;
quit;

%mend;

data _null_;
run;




/*HTRP 2D macro
Macro defined 08-27-2010
It allows for the estimation of parameters from the non-linear mixed effect model

Modified for eMKF in 2023 Q1-Q2 by Makram Talih.
In eMKF, this assumes data has been pre-formatted using macro reformat.

data 	: name of the dataset.
outcome : outcome of interest.
se      : standard errors of the outcome.
outcome2 : outcome2 of interest.
se2      : standard errors of outcome2.
time    : time variable.
by      : allows models to run for multiple strata at the same time and can be used for simulations
xtrakeep: Any variable one wants to keep in the data while runing models: weights, ... (eMKF: could be used to retain labels for multiyear data)
orpoly  : (eMKF) YES (default) for pre-transforming the design matrix using SAS IML orpol function. NO for "raw" polynomials.
          If YES, regression coefficients and their SEs will be reverse-transformed prior to macro end.
_rho_   : the value of the true rho that generated the data, if known. If not given, it will be estimated
_tausq_ : the value of the true tau-square that generated the data, if known. If not given, it will be estimated
bvalue  : Assumption about the bvalue: eMKF options are the following:
             indep_cubic	: The values of the parameters b1, b2, and b3 are computed for each group
             indep_quad		: b3=0. The values of the parameters b1 and b2 are computed for each group
             indep_linear   : (DEFAULT) b3=0 and b2=0. The value of the slope b1 is computed for each group
             common_cubic	: The values of each of the parameters b1, b2, and b3 are assumed to be the same across groups
             common_quad	: b3=0. The values of each of the parameters b1 and b2 are assumed to be the same across groups
             common_linear  : b3=0 and b2=0. The value of the slope b1 is assumed to be the same across groups
             dropped    	: A model without time trend is computed
group   : the different groups (e.g., race/ethnicity groups) variable
DF      : Non-linear model degrees of freedom. The default is set pretty high at 10000
print   : Yes will print the nlmixed results and No will not. Default is No
out     : The name of the output baseline. All the following outputs (baseline + suffix) are saved. 
          Here are the suffixes:
          	_fitstat : model fit estimates from the proc nlmixed
          	_ests    : model fit estimates formated for use in the estimation of the Kalman prediction
          	_covmat  : model fit covariance matrix
          	_pred    : Kalman prediction of the outcome of interest includes original values as well as parameters
         (e.g., for OUT=result then RESULT_PRED will be the Kalman prediction data of the outcome of interest. )
*/

%macro htrp2d(
             data=, 
             outcome=, 
             se=, 
             outcome2=, 
             se2=, 
             group=, 
             time=, 
             by=, 
             xtrakeep= ,
			 orpoly = YES,
             _rho_= , 
             _tausq_= , 
             _delta_= ,
             bvalue= indep_linear , 
             DF=10000, 
             out=param , 
             print=NO
            );

%local n g p k nrep n1 n2 i iis jj lj b1line b2line b3line tline jdelta _mdelta_ formatted dsop dscl _rtimess rtm rlag
       group_rep vmatinv vmat xmat xmat2 xmat3 xmat4 amat blmat _emkfkeep_ _emkfmu_ _emkfmu2_ 
	   parList1 colList1 parList2 colList2 _tcobs oPPmat;

/* eMKF: Data assumed to have been pre-formatted using macro reformat: check and reformat if not */

%let formatted = 0;

%let dsop = %sysfunc(open(&data));
%if &dsop ne 0 %then %do;
	%if %sysfunc(varnum(&dsop, inputorder)) ne 0 and %sysfunc(varnum(&dsop, &time)) ne 0 %then %let formatted = 1;
%end; 
%let dscl = %sysfunc(close(&dsop));

%let formatted = %eval(&formatted + 0);

data _sdata_ _jdata_;
run;

%if &formatted = 1 %then %do;
	data _sdata_;
	  set &data;
	run;
%end;
%else %do;

    %put ;
	%put Reformatting data prior to MLE-based estimation;

	%reformat(data=&data, 
		      outcome=&outcome, se=&se, outcome2=&outcome2, se2=&se2, 
			  group=&group, time=&time, by=&by, 
 			  outformat= _sdata_ );

	/*eMKF: Create copies of _group_ and _rep variables for use in proc iml matrix calculations */
	data _sdata_; 
	  set _sdata_;
	  _groupnum = _group_; 					 
	  %if &by ^= %str() %then _reps = _rep;;
	run;

	/*eMKF: Replace &group and &by macro variables by their numeric versions */
	%let group = _groupnum; 				
	%if &by ^= %str() %then %let by = _reps;;

%end;

/*eMKF: Sort by replications (if any), group, and time */
/*eMKF: Use _time as &time could be empty if data was still in format 1 when htrp was called */
proc sort data= _sdata_;
  by _rep _group_ _time ;
run;

/*eMKF: Macro variable for the number of groups */
%let g=;
data _freqg_;
run;
proc freq data=_sdata_ noprint;
  tables &group /list out=_freqg_;
run;
data _freqg_;
  set _freqg_;
  _group_ +1;
  call symput('g',_group_);
  keep _group_ &group;
run;

/*eMKF: Macro variable for the number of time points */
%let n=;
data _freqn_;
run;
proc freq data=_sdata_ noprint;
  tables _rtime /list out=_freqn_;
run;
data _freqn_;
  set _freqn_;
  _time +1;
  call symput('n',_time);
  keep _time _rtime;
run;

/*eMKF: Macro variable for the real times to use in calculations */
%let _rtimess = ;
data _freqn_;
  set _freqn_;
  retain _rts;
  if _n_= 1 then _rts = cat(_rtime);
  else _rts = catx(" ", _rts, _rtime);
  call symput('_rtimess', _rts);
  drop _rts;
run;

/* eMKF allows for irregular and fractional times points */
/* This is the variable that will be used for real time in case it is not just 1,2,3,... */
%let rtm=_rtime; 		

/*eMKF: macro reformat now also tracks lags between successive real time points */
%let rlag=_rlag;

/*eMKF: Macro variable for the number of replications */
%let nrep=1;
data _freq_;
run;
%if &by ^=%str() %then %do;
	proc freq data=_sdata_ noprint;
	  tables &by /list out=_freq_;
	  format &by ;
	run;
	data _freq_;
	  set _freq_;
	  _rep +1;
	  call symput('nrep',_rep);
	  keep _rep &by;
	run;
%end;

/*eMKF: Set numerical values to use in code */
%let n=%eval(0+&n);
%let g=%eval(0+&g);
%let nrep=%eval(0+&nrep);

/* Correlation between outcomes*/

%if %length(&_delta_) ^= 0  %then %let _delta_=%sysevalf(&_delta_ -1);;

data aa;
run;

%let jdelta=;
proc corr data=_sdata_ outp=aa noprint;
  var _y _y2;
run;

data _null_;
  set aa;
  id+1;
  if id=5 then call symput("jdelta",_y);
run;

/* End of Correlation estimates*/

/* eMKF: Modification to set up orthogonal cubic polynomial design matrix */

data _oXmat_ _oPmat_;
run;

%if %upcase(&orpoly) = YES %then %do;

	proc iml;
	  x = { &_rtimess };	
	  x = T(x);							/* eMKF: column vector with real times */
	  oP = orpol(x, 3);					/* eMKF: orthonormal design matrix for cubic orthogonal polynomials */
	  x0 = { %cnstss(1, &n) };
	  x0 = T(x0);
	  x1 = x;
	  x2 = x#x;
	  x3 = x#x2;
	  uP = x0 || x1 || x2 || x3;		/* eMKF: raw/unstandardized design matrix */
	  oP1 = inv(T(uP)*uP)*T(uP)*oP[,1];
      oP2 = inv(T(uP)*uP)*T(uP)*oP[,2];
      oP3 = inv(T(uP)*uP)*T(uP)*oP[,3];
      oP4 = inv(T(uP)*uP)*T(uP)*oP[,4];
	  oPP = oP1 || oP2 || oP3 || oP4;	/* eMKF: right multiplication of raw uP with oPP produces orthonormal oP */
	  y = T(do(1, &n, 1));				/* eMKF: column vector of consecutive time indices */
	  yP = y || oP;
	  create _oXmat_ from yP [ colname = {"_time" "&rtm.0" "&rtm.1" "&rtm.2" "&rtm.3"} ] ;
	  append from yP; close _oXmat_;
	  create _oPmat_ from oPP [ colname = {"t0" "t1" "t2" "t3"} ] ;
	  append from oPP; close _oPmat_;
	quit;

	proc sort data=_sdata_;
	  by _time;
	run;

	data _sdata_;
	   merge _sdata_ _oXmat_;
	   by _time;
	run;

	proc sort data= _sdata_;
	  by _rep _group_ _time ;
	run;

%end;
%else %do;
	data _sdata_;
	  set _sdata_;
	  &rtm.0 = 1;
	  &rtm.1 = &rtm;
	  &rtm.2 = &rtm**2;					/* eMKF: add raw quad and cubic time terms as columns in the _sdata_ dataset */
	  &rtm.3 = &rtm**3;
	run;
%end;

data _jdata_;
  set _sdata_(in=a rename=(_y =_oy _se=_ose)) _sdata_(in=b rename=(_y2 =_oy _se2=_ose));
  if a then orep=1;
  if b then orep=2;
  ones 	   = orep - 1;
  onetime0 = &rtm.0 * ones;
  onetime1 = &rtm.1 * ones;
  onetime2 = &rtm.2 * ones;
  onetime3 = &rtm.3 * ones;
run;

proc sort data=_jdata_;
  by _group_ orep;
run;

/* eMKF: Set up model statement symbolically for use in proc reg
 * Initial regression when no time trend is desired differs from the original MKF, here.  
 * In MKF, initial values for intercepts were based on those from linear trend model instead of dropped model.
 */
%let tline=;
%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = INDEP_CUBIC %then 
	%let tline = &rtm.0 onetime0 &rtm.1 onetime1 &rtm.2 onetime2 &rtm.3 onetime3;
%if %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = INDEP_QUAD %then 
	%let tline = &rtm.0 onetime0 &rtm.1 onetime1 &rtm.2 onetime2;
%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = INDEP_LINEAR %then 
	%let tline = &rtm.0 onetime0 &rtm.1 onetime1;
%if %upcase(&bvalue) = DROPPED %then 
	%let tline = &rtm.0 onetime0;

/* eMKF: Modified call to proc reg to include quad and cubic terms
 * Implemented regression by _rep instead of relying only on _rep = 1 for initial values to pass to proc nlmixed
 */
proc reg data=_jdata_ outest=_beta_ noprint;
	by _rep _group_ ;
	model _oy = &tline /noint;;
run;

/* Setup initial values (modified for eMKF to hold quad and cubic coefficients)
   Use the beta data estimates from the model above. For &g groups, then:
     a1-a&g are the intercept from each of the &g groups models
     b1arr1-b1arr&g are the linear coefficients from each of the &g groups models
     b2arr1-b2arr&g are the quadratic coefficients from each of the &g groups models
     b3arr1-b3arr&g are the cubic coefficients from each of the &g groups models
*/

/*eMKF: Sort regression coefficients by _rep and _group_ */
proc sort data=_beta_;
 by _rep _group_;
run;

/*eMKF: Initial parameter values for the first replication */
%let jj = 1;
data _inits_;
   set _beta_(where=(_rep = &jj)) end=end;
   _rep = &jj;
   array o1a o1a1-o1a&g;
   array o1b1arr o1b1arr1-o1b1arr&g;
   array o1b2arr o1b2arr1-o1b2arr&g;
   array o1b3arr o1b3arr1-o1b3arr&g;
   array o2a o2a1-o2a&g;
   array o2b1arr o2b1arr1-o2b1arr&g;
   array o2b2arr o2b2arr1-o2b2arr&g;
   array o2b3arr o2b3arr1-o2b3arr&g;
   retain o1a1-o1a&g o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g o1b3arr1-o1b3arr&g  
		  o2a1-o2a&g o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g o2b3arr1-o2b3arr&g; 
   o1a{_group_} = &rtm.0; o2a{_group_} = onetime0;
   %if %upcase(&bvalue) ^= DROPPED %then %str( o1b1arr{_group_} = &rtm.1; o2b1arr{_group_} = onetime1;) ;;
   %if %upcase(&bvalue) ^= DROPPED and %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= INDEP_LINEAR %then %str( o1b2arr{_group_} = &rtm.2; o2b2arr{_group_} = onetime2;) ;;
   %if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = INDEP_CUBIC %then %str( o1b3arr{_group_} = &rtm.3; o2b3arr{_group_} = onetime3;) ;;
   if end then do;
		%if %length(&_delta_) = 0 %then delta= &jdelta  - 1;;
		%if &_rho_ = %str() %then logitrho = 0; ;
		%if &_tausq_ = %str() %then logtau2 = log(.002);;
      	output;
   end;
   keep _rep
		o1a1-o1a&g o2a1-o2a&g o1b1arr1-o1b1arr&g o2b1arr1-o2b1arr&g o1b2arr1-o1b2arr&g o2b2arr1-o2b2arr&g o1b3arr1-o1b3arr&g o2b3arr1-o2b3arr&g
		%if %length(&_delta_) = 0  %then delta;
		%if &_rho_ = %str() %then logitrho;
		%if &_tausq_ = %str() %then  logtau2;
   ;
run;

/*eMKF: Initial estimates from each subsequent replication */
data _initsr_;
run;
%if &nrep > 1 %then %do;
	%do jj = 2 %to &nrep;
		data _initsr_;
		   set _beta_(where=(_rep = &jj)) end=end;
		   _rep = &jj;
		   array o1a o1a1-o1a&g;
		   array o1b1arr o1b1arr1-o1b1arr&g;
		   array o1b2arr o1b2arr1-o1b2arr&g;
		   array o1b3arr o1b3arr1-o1b3arr&g;
		   array o2a o2a1-o2a&g;
		   array o2b1arr o2b1arr1-o2b1arr&g;
		   array o2b2arr o2b2arr1-o2b2arr&g;
		   array o2b3arr o2b3arr1-o2b3arr&g;
		   retain o1a1-o1a&g o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g o1b3arr1-o1b3arr&g  
				  o2a1-o2a&g o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g o2b3arr1-o2b3arr&g; 
		   o1a{_group_} = &rtm.0; o2a{_group_} = onetime0;
		   %if %upcase(&bvalue) ^= DROPPED %then %str( o1b1arr{_group_} = &rtm.1; o2b1arr{_group_} = onetime1;) ;;
		   %if %upcase(&bvalue) ^= DROPPED and %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= INDEP_LINEAR %then %str( o1b2arr{_group_} = &rtm.2; o2b2arr{_group_} = onetime2;) ;;
		   %if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = INDEP_CUBIC %then %str( o1b3arr{_group_} = &rtm.3; o2b3arr{_group_} = onetime3;) ;;
		   if end then do;
				%if %length(&_delta_) = 0 %then delta= &jdelta  - 1;;
				%if &_rho_ = %str() %then logitrho = 0; ;
				%if &_tausq_ = %str() %then logtau2 = log(.002);;
		      	output;
		   end;
		   keep _rep
				o1a1-o1a&g o2a1-o2a&g o1b1arr1-o1b1arr&g o2b1arr1-o2b1arr&g o1b2arr1-o1b2arr&g o2b2arr1-o2b2arr&g o1b3arr1-o1b3arr&g o2b3arr1-o2b3arr&g
				%if %length(&_delta_) = 0  %then delta;
				%if &_rho_ = %str() %then logitrho;
				%if &_tausq_ = %str() %then  logtau2;
		   ;
		run;

		data _inits_;
		  set _inits_ _initsr_;
		run;

		data _initsr_;
		run;
	%end;
%end;

/*eMKF: calculate initial estimates for common_ scenarios and remove extraneous variables */
data _inits_;
  set _inits_;
  o1b1=.; 
  o1b2=.; 
  o1b3=.; 
  o2b1=.;
  o2b2=.;
  o2b3=.;
  %if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_LINEAR %then %str( o1b1 = mean(of o1b1arr1-o1b1arr&g); o2b1 = mean(of o2b1arr1-o2b1arr&g);) ;;
  %if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD and %upcase(&bvalue) ^= INDEP_LINEAR %then %str( drop o1b1arr1-o1b1arr&g o2b1arr1-o2b1arr&g;) ;;
  %if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD %then %str( o1b2 = mean(of o1b2arr1-o1b2arr&g); o2b2 = mean(of o2b2arr1-o2b2arr&g);) ;;
  %if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD %then %str( drop o1b2arr1-o1b2arr&g o2b2arr1-o2b2arr&g;) ;;
  %if %upcase(&bvalue) = COMMON_CUBIC %then %str( o1b3 = mean(of o1b3arr1-o1b3arr&g); o2b3 = mean(of o2b3arr1-o2b3arr&g);) ;;
  %if %upcase(&bvalue) ^= INDEP_CUBIC %then %str( drop o1b3arr1-o1b3arr&g o2b3arr1-o2b3arr&g;) ;;
  %if %upcase(&bvalue) ^= COMMON_CUBIC and %upcase(&bvalue) ^= COMMON_QUAD and %upcase(&bvalue) ^= COMMON_LINEAR %then drop o1b1 o2b1;;;
  %if %upcase(&bvalue) ^= COMMON_CUBIC and %upcase(&bvalue) ^= COMMON_QUAD %then drop o1b2 o2b2;;;
  %if %upcase(&bvalue) ^= COMMON_CUBIC %then drop o1b3 o2b3;;;
run;

/* eMKF: re-order columns in common trend cases so delta, logitrho, and logtau2 remain last and coefficients are grouped by outcome*/
%if %upcase(&bvalue)=COMMON_CUBIC or %upcase(&bvalue)=COMMON_QUAD or %upcase(&bvalue)=COMMON_LINEAR %then %do;
	data _inits_;
	  merge _inits_(drop = o2: delta logitrho logtau2) _inits_(keep = _rep o2:) _inits_(keep = _rep delta logitrho logtau2);
	  by _rep;
	run;
%end;

/*eMKF: Set up array declarations symbolically for use in proc nlmixed */
%let b1line=; %let b2line=; %let b3line=;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then 
	%let b1line= array o1b1arr(&g) o1b1arr1-o1b1arr&g%str(;) array o2b1arr(&g) o2b1arr1-o2b1arr&g%str(;) ;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then 
	%let b2line= array o1b2arr(&g) o1b2arr1-o1b2arr&g%str(;) array o2b2arr(&g) o2b2arr1-o2b2arr&g%str(;) ;
%if %upcase(&bvalue) = INDEP_CUBIC %then 
	%let b3line= array o1b3arr(&g) o1b3arr1-o1b3arr&g%str(;) array o2b3arr(&g) o2b3arr1-o2b3arr&g%str(;) ;

/*eMKF: Set up model statement symbolically for use in proc nlmixed */
%let _mdelta_ = delta;
%if %length(&_delta_) ^= 0 %then %let _mdelta_ = &_delta_;;
%let tline=;
%if %upcase(&bvalue) = INDEP_CUBIC %then 
	%let tline= normal(&rtm.0*o1a[_group_] + o2a[_group_]*ones + &rtm.1*o1b1arr[_group_] + onetime1*o2b1arr[_group_] 
						+ &rtm.2*o1b2arr[_group_] + onetime2*o2b2arr[_group_] + &rtm.3*o1b3arr[_group_] + onetime3*o2b3arr[_group_]
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = INDEP_QUAD %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] + &rtm.1*o1b1arr[_group_] + onetime1*o2b1arr[_group_] 
						+ &rtm.2*o1b2arr[_group_] + onetime2*o2b2arr[_group_] 
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = INDEP_LINEAR %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] + &rtm.1*o1b1arr[_group_] + onetime1*o2b1arr[_group_] 
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = COMMON_CUBIC %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] + &rtm.1*o1b1 + onetime1*o2b1 
						+ &rtm.2*o1b2 + onetime2*o2b2 + &rtm.3*o1b3 + onetime3*o2b3
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = COMMON_QUAD %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] + &rtm.1*o1b1 + onetime1*o2b1 
						+ &rtm.2*o1b2 + onetime2*o2b2
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = COMMON_LINEAR %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] + &rtm.1*o1b1 + onetime1*o2b1 
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);
%if %upcase(&bvalue) = DROPPED %then 
	%let tline= normal(&rtm.0*o1a[_group_] + onetime0*o2a[_group_] 
						+ gamma[_time] + &_mdelta_*ones*gamma[_time], _ose**2);

/* Fit a Non-linear mixed model. Capture covariance matrix COV and covariance matrix from additional estimate ECOV */

data _ests _fitstat &out._fitstat &out._covmat &out._data;
run;

data &out._data;
 set _sdata_;
run;

%put Start model fitting using PROC NLMIXED; /* eMKF: substituted put statement */
%put Model is _oy ~ &tline;

/*eMKF: model print handling */
%if %upcase(&print) ^=YES %then ods exclude all;;

/*eMKF: initical counter for use in random statement */
%let i=0;

proc nlmixed data=_jdata_ DF=&DF cov ecov 
	method=firo maxiter=500; 			/* eMKF: increased maxiter to 500 from the 200 default for dealing with cubic (tech=QUANEW default)*/
	by _rep;							/* eMKF: stratified by _rep (trivial if _rep is constant = 1 (i.e., no stratification) */
	array gamma(&n) gamma1-gamma&n ; 	/* Give the gamma values a generic name */
	array o1a(&g) o1a1-o1a&g; 			/* Give the a values a generic name */
	array o2a(&g) o2a1-o2a&g; 			/* Give the a values a generic name */
	&b1line;; 						 	/* Give the bl values a generic name */
	&b2line;; 						 	/* eMKF: Give the b2 values a generic name */
	&b3line;; 						 	/* eMKF: Give the b3 values a generic name */ 
	parms / bydata data=_inits_; 		/* eMKF: Setup parameters and their initial values stratified by _rep */
	/* 
	Next define Tausq to be an assigned value of tausq but in the case &_tausq_ is missing 
	 then assign value exp(logtau2) a value that will be estimated.
	Same thing for _rho_ taking assigned value but in case &_rho_is missing, 
	 give it the value 2/(1+exp(-logitrho)) - 1
	*/
	%if &_tausq_ = %str() %then _tausq_ = exp(logtau2); ;
	%if &_tausq_ ^= %str() %then _tausq_ = &_tausq_; ;
	%if &_rho_ = %str() %then _rho_ = 2/(1+exp(-logitrho)) - 1;;
	%if &_rho_ ^= %str() %then _rho_ = &_rho_; ;
	nu = _tausq_/(1-(_rho_**2));

	/* Random effects gamma1 gamma2 ... gamma&n ~ normal ([0,0,0,..,0],[nu, half triangular matrix]) */
	/* eMKF: Replaced call to thevarcomp with call to thevarcompr to allow noninteger/unequally-spaced time points */
	random %do i = 1 %to &n; gamma&i %end;
    ~ normal ([%zeros(&n)], %thevarcompr(times=&_rtimess, nu=nu, vrho=_rho_)) subject=_group_ out=&out._BLUP;

	model _oy ~ &tline;;

	%if &_rho_ = %str() %then estimate "_rho_" 2/(1+exp(-logitrho))-1;;
	%if &_tausq_ = %str() %then estimate "_tausq_" exp(logtau2); ;

	ods output parameterestimates=_ests fitstatistics=_fitstat;
	ods output CovMatParmEst=&out._covmat;
run;

/*eMKF: Replace missing values in CovMatParmEst with 0s for later use */
%if %upcase(&orpoly) = YES %then %do;
	data &out._covmat;
	  set &out._covmat;
	  array NAs _numeric_;
	  do over NAs;
	      if NAs = . then NAs = 0;
	  end;
	run;
%end;

/*eMKF: Reset ODS destinations if they were turned off */
%if %upcase(&print) ^=YES %then ods exclude none;;

/* eMKF: Added put statements for log */
%put End model fitting using PROC NLMIXED;
%if %upcase(&print) ^= YES %then %put Model printout was turned off by user;;

/* eMKF: Model estimates */

data _fitstat;
   set _fitstat;
   if descr = "-2 Log Likelihood";
   _2loglike=value;
   keep _2loglike _rep;
run;

data &out._fitstat;
   merge _ests _fitstat;
   by _rep;
run;

data _llike_;
run;
data _llike_;
 set &out._fitstat;
 by _rep;
 if first._rep then output;
 keep _rep _2loglike;
run;

data &out._ests;
run;
proc transpose data=&out._fitstat out=&out._ests;
   by _rep;
   var estimate;
   id parameter;
run;

/* eMKF: Macro variables _emkfkeep_, _emkfmu_, and _emkfmu2_ account for various model combinations */
%let _emkfkeep_ = ; 
%let _emkfmu_  = ; 
%let _emkfmu2_ = ;
%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC %then %do;
	%let _emkfkeep_ = o1a o2a o1b1 o2b1 o1b2 o2b2 o1b3 o2b3;
	%let _emkfmu_  = mu  = o1a*&rtm.0 + o1b1*&rtm.1 + o1b2*&rtm.2 + o1b3*&rtm.3;
	%let _emkfmu2_ = mu2 = o2a*&rtm.0 + o2b1*&rtm.1 + o2b2*&rtm.2 + o2b3*&rtm.3;
%end;
%if %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD %then %do;
	%let _emkfkeep_ = o1a o2a o1b1 o2b1 o1b2 o2b2;
	%let _emkfmu_  = mu  = o1a*&rtm.0 + o1b1*&rtm.1 + o1b2*&rtm.2;
	%let _emkfmu2_ = mu2 = o2a*&rtm.0 + o2b1*&rtm.1 + o2b2*&rtm.2;
%end;
%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = COMMON_LINEAR %then %do;
	%let _emkfkeep_ = o1a o2a o1b1 o2b1;
	%let _emkfmu_  = mu  = o1a*&rtm.0 + o1b1*&rtm.1;
	%let _emkfmu2_ = mu2 = o2a*&rtm.0 + o2b1*&rtm.1;
%end;
%if %upcase(&bvalue) = DROPPED %then %do;
	%let _emkfkeep_ = o1a o2a; 
	%let _emkfmu_  = mu  = o1a*&rtm.0; 
	%let _emkfmu2_ = mu2 = o2a*&rtm.0;
%end;

/*eMKF: Re-structure estimates by _rep and _group_ */
data &out._ests;
   merge &out._ests _llike_;
   by _rep;
   array o1as o1a1-o1a&g;
   array o2as o2a1-o2a&g;
   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %str(array o1b1s o1b1arr1-o1b1arr&g; array o2b1s o2b1arr1-o2b1arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(o1b1 = 0; o2b1 = 0;) ;;
   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then %str(array o1b2s o1b2arr1-o1b2arr&g; array o2b2s o2b2arr1-o2b2arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(o1b2 = 0; o2b2 = 0;) ;;
   %if %upcase(&bvalue) = INDEP_CUBIC %then %str(array o1b3s o1b3arr1-o1b3arr&g; array o2b3s o2b3arr1-o2b3arr&g;) ;;
   %if %upcase(&bvalue) = DROPPED %then %str(o1b3 = 0; o2b3 = 0;) ;;
   %if %length(&_delta_) ^= 0  %then delta=&_delta_;;;
   _rho_ = &_rho_ %if &_rho_ = %str() %then 2/(1+exp(-logitrho))-1  ;;;
   _tausq_ = &_tausq_ %if &_tausq_ = %str() %then exp(logtau2) ;;;
   do _group_ = 1 to &g;
      o1a = o1as{_group_};
      o2a = o2as{_group_};
      %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %str(o1b1 = o1b1s{_group_}; o2b1 = o2b1s{_group_};) ;;
	  %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD %then %str(o1b2 = o1b2s{_group_}; o2b2 = o2b2s{_group_};) ;;
      %if %upcase(&bvalue) = INDEP_CUBIC %then %str(o1b3 = o1b3s{_group_}; o2b3 = o2b3s{_group_};) ;;
      output;
      end;
   keep _rep _group_ &_emkfkeep_ _rho_ _tausq_ _2loglike delta;;
run;

/*eMKF: Reset contrasts to get outcome-specific estimates */
data &out._ests;
  set &out._ests;
  o2a = o2a + o1a;
  %if %upcase(&bvalue) ^= DROPPED %then %str(o2b1 = o2b1 + o1b1;);;
  %if %upcase(&bvalue) ^= DROPPED and %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= INDEP_LINEAR %then %str(o2b2 = o2b2 + o1b2;);;
  %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC %then %str(o2b3 = o2b3 + o1b3;);;
  delta = 1 + delta;
run;

/* eMKF: Merge with group labels */
proc sort data= &out._ests;
 by _group_ ;
run;
proc sort data=_freqg_;
 by _group_;
run;
data &out._ests;
 merge &out._ests _freqg_;
  by _group_;
run;

/*eMKF: Re-order columns so that intercept is always listed first */
data &out._ests; 
  merge &out._ests(drop= &_emkfkeep_) &out._ests(keep= o1a o2a) %if %length(&_emkfkeep_) > 7 %then &out._ests(keep= %substr(&_emkfkeep_, 9));;; 
run;

/***************************************************************************/
/* Now let's use the Kalman technique to estimate the outcome observations */
/***************************************************************************/

data &out._pred;
run;
proc sort data=_sdata_ out=&out._pred;
  by _rep _group_ _time;
run;

proc sort data=&out._ests;
  by _rep _group_;
run;

data &out._pred;

  merge &out._pred &out._ests;
  by _rep _group_ ;

  &_emkfmu_;; 		/*eMKF: invoke symbolic calculation for mu */
  err=_y - mu;

  &_emkfmu2_;; 		/*eMKF: invoke symbolic calculation for mu2 */
  err2=_y2 - mu2;

  _tausq2_= _tausq_ * (delta**2);
  _rho2_= _rho_ ;; * delta ; *This was a mistake that was fixed but need to be double checked again ;

  err3 = (err + (err2/delta) )/2;  *Add this one to use the KF on combined errors;
  _se3 = sqrt((_se**2 + (_se2/delta)**2 )/2); * The new addition;

  * These next steps are only needed if we assume the same AR(1) processes up to a scale. 
    For now, we are relaxing that assumption that both series only help estimate rho;
  *_se=_se3; *reset  standard error out assumption;
  *_se2 = delta*_se3; *reset  standard error out assumption;
  *err=err3; *reset  standard error out assumption;
  *err2=delta*err3; *reset  standard error out assumption;

  rename delta=moddelta;

run;

data _empty_;
  set &out._pred;
  if _time=1;
  w01 = _tausq_/(1-(_rho_**2));
  gamma01 = 0 ;
  w02 = _tausq2_/(1-(_rho2_**2));
  gamma02 = 0 ;
  _time=0;
  &rtm = 0;
  &rlag = .;
  keep _time &rtm &rtm.0 &rtm.1 &rtm.2 &rtm.3 &rlag &by _rep &group _group_ w01 gamma01 w02 gamma02;
run;

data &out._pred;
 set &out._pred _empty_;
 if _time = 1 then &rlag = 1;	/*eMKF: set to lag 1 instead of missing relative to time 0 */
run;

proc sort data=&out._pred;
 by _rep _group_ _time;
run;

/*eMKF recursion formulas modified to allow lag s > 0 between time points */
data &out._pred;

  set &out._pred;
  retain wold gamold wold2 gamold2;

  delta = ((_rho_**(2*&rlag)) * wold) + _tausq_*(1 - (_rho_**(2*&rlag)))/(1 - (_rho_**2));
  if &rlag > 0 or &rlag = . then lambda = delta/(delta + (_se**2));
  else lambda = 0; 	/*eMKF: in the limiting case, recursion breaks down: set lambda = 0 instead */
  w = delta * (1-lambda);
  gamma = (lambda* err) + (1-lambda)*(_rho_**&rlag)*gamold ; /*eMKF: &rlag=1 reduces to unit-increment case */
  prediction = mu + gamma ;

  delta2 = ((_rho2_**(2*&rlag)) * wold2) + _tausq2_*(1 - (_rho2_**(2*&rlag)))/(1 - (_rho2_**2));
  if &rlag > 0 or &rlag = . then lambda2 = delta2/(delta2 + (_se2**2));
  else lambda2 = 0;
  w2 = delta2 * (1-lambda2);
  gamma2 = (lambda2* err2) + (1-lambda2)*(_rho2_**&rlag)* gamold2 ;
  prediction2 = mu2 + gamma2 ;

  output;
 
  if _time =0 then do;
 	 wold=w01;
	 gamold=gamma01;
	 wold2=w02;
	 gamold2=gamma02;
  end;

  if _time ^=0 then do;
 	 wold=w;
  	 gamold=gamma;
	 wold2=w2;
 	 gamold2=gamma2;
  end;

run;

%let group_rep = compress(_rep || _group_);

data &out._pred;
  set &out._pred;
  group_rep= &group_rep;
  if _time ne 0;
  drop w01 gamma01 wold gamold w02 gamma02 wold2 gamold2 ;
run;

/* Let's re-estimate gamma as a weighted linear combination of the two gamma and use that for prediction:
 For now, let's use weight =1/2 according to the assumption of no covariance between the gammas and 
 that the variances are the same: If different variance, then we should use 
 \lambda= (1/sigma1)/(1/sigma1 +1/sigma2) and gamma= lambda gamma1 + (1- lambda) gamma2.
*/

data &out._pred;
  set &out._pred;
  gammaeff= (gamma + (gamma2 / moddelta))/2; *Take the average for now;
  predeff = mu + gammaeff ;
  predeff2 = mu2 + (moddelta*gamma2) ;
run;

/*eMKF: reset lag for time 1 to missing */
data &out._pred;
  set &out._pred;
  if _time = 1 then &rlag = .; 
run;

/******************/
/* Compute the MSE*/
/******************/

data _Amat_ _Dmat_ _Vmat_ _Vgmat_ _Vemat_  _Xmat_ _junk_ _junk0_ _junk01_;
run;

proc sort data=&out._pred(keep= _rep &by _group_ &group _time &rtm 
								&rtm.0 &rtm.1 &rtm.2 &rtm.3 		/*eMKF: modification to deal with orthogonal polynomials */
								&rlag _y _se _y2 _se2 _rho_ _rho2_  
								_tausq_ _tausq2_ err err2 lambda lambda2 prediction prediction2 group_rep)
  out=_Amat_ ;
  by group_rep;
run;

data _Amat_;
  set _Amat_; 
  by group_rep;
  array hs ah1 - ah&n (%zeross(&n) );
  array o2hs o2ah1 - o2ah&n (%zeross(&n) );
  if first.group_rep then do _k=1 to &n; hs{_k}=0; o2hs{_k}=0; end;
  if _time = 1 then do; ah1=lambda; o2ah1=lambda2; end;
  else do;
    /*eMKF: updated to account for &rlag if not 1 */
  	do _j=1 to _time -1;
	    hs{_j} 	 = hs{_j}  *(1-lambda) *(_rho_**&rlag) ;
	    o2hs{_j} = o2hs{_j}*(1-lambda2)*(_rho2_**&rlag) ;
  	end;
  	hs{_time}	= lambda;
  	o2hs{_time} = lambda2;
 end;

 drop group_rep _k _j;
 keep _rep &by _group_ &group _time &rtm &rlag ah: o2ah:;
run;

proc sort data=_Amat_;
 by _rep _group_ _time  ;
run;

/* Attach row number _id to each of the A Matrix so that when computing group by group
   the appropriate Ag could be called*/
data _Amat_;
  set _Amat_;
  _id+1;
run;

data _Dmat_;
  set &out._pred(keep= _rep &by _group_ &group _time &rtm 
						&rtm.0 &rtm.1 &rtm.2 &rtm.3  /*eMKF: modification to deal with orthogonal polynomials */
						&rlag _y _y2 _rho_ _rho2_ _tausq_ _tausq2_ _se _se2);
  w0   = _tausq_ /(1- (_rho_**2));
  o2w0 = _tausq2_/(1- (_rho2_**2));
run;

/* The Matrix _Dmat_ is actually the same within group and within replication*/
proc sort data= _Dmat_ nodupkey;
  by _rep _group_ _time;
run;

/* Here the standard error is the same from one replication to another
 by different from group to group. So these matrices should just be estimated 
 differently for each group as well as replication if possible*/

data _Vmat_;
  set _Dmat_;
  array ad ad1-ad&n;
  array o2ad o2ad1-o2ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n;
  	 if i = _time then do;  		/* Add the variances of Y to the diagonal elements */
		ad{i}  =   w0 + _se**2; 
		o2ad{i}= o2w0 + _se2**2; 
 	 end;
  	 else do; 
		ad{i}  =   w0*( _rho_**abs(rt{i} - &rtm));  
		o2ad{i}= o2w0*(_rho2_**abs(rt{i} - &rtm)); 
	 end;
  end;
  drop i rt1-rt&n;
  keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n o2ad1-o2ad&n;
run;

/* Recreate the variance where only the Variance of gamma is estimated
 This is mostly needed for check of what the estimates are giving us
*/
data _Vgmat_;
  set _Dmat_;
  array ad ad1-ad&n;
  array o2ad o2ad1-o2ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n;
  	 if i = _time then do; 
		ad{i}  =   w0 ; 
		o2ad{i}= o2w0 ; 
 	 end;
  	 else do; 
		ad{i}  =   w0*( _rho_**abs(rt{i} - &rtm));  
		o2ad{i}= o2w0*(_rho2_**abs(rt{i} - &rtm)); 
	 end;
  end;
  drop i rt1-rt&n;
  keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n o2ad1-o2ad&n;
run;

/* This is the diagonal Variance matrix of the errors */
data _Vemat_;
 set _Dmat_;
 array ad ad1-ad&n;
 array o2ad o2ad1-o2ad&n;
 do i=1 to &n;
 	if i=_time then do; 
		ad{i}= _se**2; 
		o2ad{i}= _se2**2;   
	end;
  	else do; 
		ad{i}=0; 
		o2ad{i}=0; 
	end;
 end;
 drop i;
 keep _rep &by _group_ &group _time &rtm &rlag ad1-ad&n o2ad1-o2ad&n;
run;

/* This is just like the _Vgmat_ but with extra variables kept to be used later
 for the creation of the X matrix for example*/
data _Dmat_;
  set _Dmat_;
  array ad ad1-ad&n;
  array o2ad o2ad1-o2ad&n;
  array rt rt1-rt&n (&_rtimess); 	/*eMKF: modification to deal with irregular time points */
  do i=1 to &n;
  	 if i = _time then do; 
		ad{i}  =   w0 ; 
		o2ad{i}= o2w0 ; 
 	 end;
  	 else do; 
		ad{i}  =   w0*( _rho_**abs(rt{i} - &rtm));  
		o2ad{i}= o2w0*(_rho2_**abs(rt{i} - &rtm)); 
	 end;
  end;
  drop i rt1-rt&n;
run;

/* This is the X matrix */
data _Xmat_;
  set _Dmat_;
  /*eMKF: modification to deal with orthogonal polynomials */
  x0 = &rtm.0; 
  x1 = &rtm.1;
  x2 = &rtm.2;
  x3 = &rtm.3;
  keep _rep &by _group_ &group x0 x1 x2 x3;
run;

data &out._H &out._PredVar &out._CovY;
run;

%let iis=;
%let vmatinv= ;
%let vmat=;
%let xmat= ;
%let xmat2= ;
%let xmat3= ; /* eMKF: xmat3 and xmat4 added to deal with quad and cubic trends */
%let xmat4= ;
%let amat= ;
%let blmat= ;

/* Now let's capture the A Matrix and turn it into a diagonal block matrix  */

data _junk_;
  set _Amat_;
  by _rep;
  if first._rep then sid=0;
  sid+1;
  repgrp=compress(_rep ||"-"|| _group_);
run;

proc sort data=_junk_;
  by repgrp sid;
run;

data _junk_;
  set _junk_;
  by repgrp;
  if first.repgrp then kp=1;
  if last.repgrp then kp=2;
  if kp=1 or kp=2;
  keep _rep &by _group_ &group _time &rtm &rlag sid kp repgrp; /*eMKF: also keeping &rtm and &rlag */
run;

data _junk_;
  set _junk_;
  by repgrp;
  retain minid maxid;
  array Aid(1:2) minid maxid;
  if first.repgrp then do;
    do i = 1 to 2;
      Aid[i] = .; /*initializing to missing*/
    end;
  end;
  Aid(kp) = sid; 
  if last.repgrp then output; 
  drop kp sid i;
run;

proc sort data=_junk_;
  by _rep _group_ _time;
run;

%let jj = 0; %let lj = 0;
data _junk_;
  set _junk_;
  by _rep;
  mm=0;
  if last._rep then mm=1;
  amat=compress("Z["||minid||":"||maxid||",]"); 			/* Here Z will be the standard matrix definition that can be used */
  if mm=0 then amat=compress(amat||",");
  if first._rep or mm=0 then amat=compress("block("||amat); /* eMKF: modified to allow for arbitrary number of groups */
  call symput("jj", _n_);
  if last._rep then call symput("lj", length(amat));		/* eMKF: added to capture maximal length of character variable needed */
run;

/* eMKF: modification to allow for arbitrary number of groups (function block restricted to 15) */
data _junk_;
  set _junk_;
  if _n_=&jj then mm=2;
  length amat2 $ %eval(&g+&lj+10);							/* eMKF: set large enough character length to cover block( + lj + g times ) */
  amat2 = amat;
  %if &g = 1 %then %do;
  	  if mm=1 then amat2=compress(amat||") ");
  	  if mm=2 then amat2=compress(amat||") ");
  %end;
  %if &g > 1 %then %do;
      if mm=1 then amat2=compress(amat||repeat(") ", &g-2));
      if mm=2 then amat2=compress(amat||repeat(") ", &g-2));
  %end;
  drop amat;
  rename amat2 = amat;
run;

proc sort data=_junk_;
  by _rep _group_;
run;

data _junk_;
  set _junk_;
  name=compress("&group"||_group_);
run;

data _junk0_ _junk01_;
run;

proc transpose data=_junk_(keep=name _rep &by amat) out=_junk0_;
  var amat;
  by _rep;
  id name;
run;

proc sort data=_junk_ out=_junk01_ nodupkey;
 by _group_;
run;

%let amat=;
proc sql noprint;
   select name into :amat separated by ' || ' from _junk01_;
quit;

data _junk0_;
  set _junk0_;
  i+1;
  code=compress("Z"||i);
  mcode=compress("_Z"||i);
  amat=compress(code||"="||&amat);
  keep amat _rep  code mcode;
run;

%let jj=0;
%do jj=1 %to &nrep;
	%local _Z&jj ;
%end;

data _null_;
  set _junk0_;
  call symput(mcode,amat);
run;

/* eMKF: call to function block in RAND's MKF macro limited the number of blocks (groups) to 15, which produced an error if &g > 15
 * eMKF: code was revised to allow for additional blocks (groups) using nested calls to function block
 */
%let jj=0;
%do jj=1 %to &g;

	%if &jj = 1 %then %let iis = i_&n ;
	%if &jj > 1 %then %let iis = &iis // i_&n ;

	/* eMKF: Setup an (ng)x(ng) matrix to capture the full V for all groups */
	%if &jj = 1 %then %let vmat = block( V ;
	%if &jj > 1 and &jj < &g %then %let vmat = &vmat , block ( V ;
	%if &jj = &g and &g > 1  %then %let vmat = &vmat , V %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let vmat = &vmat );

	/* eMKF: V inverse */
	%if &jj = 1 %then %let vmatinv = block( invV ;
	%if &jj > 1 and &jj < &g %then %let vmatinv = &vmatinv , block ( invV ;
	%if &jj = &g and &g > 1  %then %let vmatinv = &vmatinv , invV %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let vmatinv = &vmatinv );

	/* Standard way to create the block matrices variable is B*/
	%if &jj = 1 %then %let blmat = block( B ;
	%if &jj > 1 and &jj < &g %then %let blmat = &blmat , block ( B ;
	%if &jj = &g and &g > 1  %then %let blmat = &blmat , B %sysfunc(repeat( %str(%)), &g-2));
	%if &jj = &g and &g = 1  %then %let blmat = &blmat );

	/* Look into the value of X that will be used in the different scenarios*/
	/* eMKF: modified to account for quadratic and cubic scenarios */
	%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_LINEAR %then %do; 
		/* eMKF: general form for independent trends */
		%if &jj = 1 %then %let xmat = block( X ;
		%if &jj > 1 and &jj < &g %then %let xmat = &xmat , block ( X ;
		%if &jj = &g and &g > 1  %then %let xmat = &xmat , X %sysfunc(repeat( %str(%)), &g-2));
		%if &jj = &g and &g = 1  %then %let xmat = &xmat );
  	%end;
  	%if %upcase(&bvalue) ^= INDEP_CUBIC and %upcase(&bvalue) ^= INDEP_QUAD  and %upcase(&bvalue) ^=INDEP_LINEAR %then %do; 
		/*eMKF: intercept always included */
		%if &jj = 1 %then %let xmat = block( X0 ;
		%if &jj > 1 and &jj < &g %then %let xmat = &xmat , block ( X0 ;
		%if &jj = &g and &g > 1  %then %let xmat = &xmat , X0 %sysfunc(repeat( %str(%)), &g-2));
		%if &jj = &g and &g = 1  %then %let xmat = &xmat );
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_LINEAR %then %do; 
		/*eMKF: linear time */
		%if &jj = 1 %then %let xmat2 = X1 ;
		%if &jj > 1 %then %let xmat2 = &xmat2 // X1 ;
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC or %upcase(&bvalue) = COMMON_QUAD %then %do; 
		/*eMKF: quadratic time */
		%if &jj = 1 %then %let xmat3 = X2 ;
		%if &jj > 1 %then %let xmat3 = &xmat3 // X2 ;
  	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC %then %do; 
		/*eMKF: cubic time */
		%if &jj = 1 %then %let xmat4 = X3 ;
		%if &jj > 1 %then %let xmat4 = &xmat4 // X3 ;
  	%end;

%end;  /* End of &jj */

/* eMKF: modified to account for quadratic and cubic scenarios */
%if %upcase(&bvalue) = COMMON_LINEAR %then  %let xmat = &xmat || ( &xmat2 );
%if %upcase(&bvalue) = COMMON_QUAD   %then  %let xmat = &xmat || ( &xmat2 ) || ( &xmat3 );
%if %upcase(&bvalue) = COMMON_CUBIC  %then  %let xmat = &xmat || ( &xmat2 ) || ( &xmat3 ) || ( &xmat4 );

%let jj=0;
%let n1=2;
%let n2= %eval(&n + 1);

/* eMKF: added dimensionality p for ease of coding of quadratic and cubic trend models in proc iml */
%let p=0;
%if %upcase(&bvalue) = INDEP_CUBIC  or %upcase(&bvalue) = COMMON_CUBIC  %then %let p = 4;
%if %upcase(&bvalue) = INDEP_QUAD   or %upcase(&bvalue) = COMMON_QUAD   %then %let p = 3;
%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = COMMON_LINEAR %then %let p = 2;
%if %upcase(&bvalue) = DROPPED 											%then %let p = 1;

proc iml;

	 %do jj=1 %to &nrep;

		/* Next for the Xmatrix */
		/* eMKF: modified to allow for quadratic and cubic trend models */
		use _Xmat_(where=(_rep=&jj) keep= _rep x:);
		read all into XX; close _Xmat_; /* eMKF: added close statements for cleanliness */
		Z =XX[,2:(1+&p)];
		&&_Z&jj ;;
	 	Xs = Z&jj;

		/* Next for the V matrix */
		use _Vmat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VV; close _Vmat_;
		Z =VV[,&n1:&n2];
		&&_Z&jj ;;
		Vs = Z&jj;
		invVs=inv(Vs);

		/* The o2V matrix for outcome 2*/
		use _Vmat_(where=(_rep=&jj) keep= _rep o2ad:);
		read all into o2VV; close _Vmat_;
		Z =o2VV[,&n1:&n2];_;
		&&_Z&jj ;;
		o2Vs = Z&jj;
		invo2Vs=inv(o2Vs);

		/* Next for the Ve matrix */
		use _Vemat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VVe; close _Vemat_;
		Z =VVe[,&n1:&n2];
		&&_Z&jj ;;
		Ves = Z&jj;

		/* The Ve matrix for outcome 2*/
		use _Vemat_(where=(_rep=&jj) keep= _rep o2ad:);
		read all into o2VVe; close _Vemat_;
		Z =o2VVe[,&n1:&n2];
		&&_Z&jj ;;
		o2Ves = Z&jj;

		/* Next for the Vg matrix */
		use _Vgmat_(where=(_rep=&jj) keep= _rep ad:);
		read all into VVg; close _Vgmat_;
		Z =VVg[,&n1:&n2];
		&&_Z&jj ;;
		Vgs = Z&jj;

		/* The Vg matrix for outcome 2*/
		use _Vgmat_(where=(_rep=&jj) keep= _rep o2ad:);
		read all into o2VVg; close _Vgmat_;
		Z =o2VVg[,&n1:&n2];
		&&_Z&jj ;;
		o2Vgs = Z&jj;

		/* Next for the A matrix */
		use _Amat_(where=(_rep=&jj)  keep= _rep ah:);
		read all into AA; close _Amat_;
		Z =AA[,&n1:&n2];
		&&_Z&jj ;;
		As = Z&jj;

		/* The A matrix for outcome 2*/
		use _Amat_(where=(_rep=&jj)  keep= _rep o2ah:);
		read all into o2AA; close _Amat_;
		Z =o2AA[,&n1:&n2];
		&&_Z&jj ;;
		o2As = Z&jj;

		/* Next for the data label matrix */
		use _Dmat_(where=(_rep=&jj) keep=_rep &by &group &rtm _time _y _se _y2 _se2); 	/*eMKF: also keeping &rtm */
		read all var{_y} into Y;
		read all var{_y2} into o2Y;
		read all var{_rep &by &group &rtm _time} into NM;							/*eMKF: also keeping &rtm */
		close _Dmat_;

		/* Now do the estimations */
		i_&n = i(&n*&g);
		H = Xs * inv(t(Xs)*invVs*Xs)*t(Xs)*invVs;
		fH = H + As*(i_&n - H);
		fY = fH * Y;
		Vy = vecdiag(fH * Vs * t(fH));
		MSEy = vecdiag(  (fH - i_&n) * Vgs * t(fH - i_&n)   ) + vecdiag(fH * Ves * t(fH));

		/* Then for outcome 2 */
		o2H = Xs * inv(t(Xs)*invo2Vs*Xs)*t(Xs)*invo2Vs;
		fo2H = o2H + o2As*(i_&n - o2H);
		fo2Y = fo2H * o2Y;
		Vo2y = vecdiag(fo2H * o2Vs * t(fo2H));
		MSEo2y = vecdiag(  (fo2H - i_&n) * o2Vgs * t(fo2H - i_&n)   ) + vecdiag(fo2H * o2Ves * t(fo2H));

		ff= NM || fH;
		o2ff= NM || fo2H;
		fV= NM || Vs;
		o2fV= NM || o2Vs;

		fVy=NM || fY || Vy || MSEy || fo2Y || Vo2y || MSEo2y;

		%if &jj = 1 %then ffs= ff;;
		%if &jj = 1 %then o2ffs= o2ff;;
		%if &jj = 1 %then fVs= fV;;
		%if &jj = 1 %then o2fVs= o2fV;;
		%if &jj = 1 %then fVys= fVy;;
		%if &jj > 1 %then ffs= ffs // ff;;
		%if &jj > 1 %then o2ffs= o2ffs // o2ff;;
		%if &jj > 1 %then fVs= fVs // fV;;
		%if &jj > 1 %then o2fVs= o2fVs // o2fV;;
		%if &jj > 1 %then fVys= fVys // fVy;;

	 %end;

	 create &out._H from ffs ;
	 append from ffs; close &out._H;
	 create &out._o2H from o2ffs ;
	 append from o2ffs; close &out._o2H;

	 create &out._CovY from fVs ;
	 append from fVs; close &out._CovY;
	 create &out._o2CovY from o2fVs ;
	 append from o2fVs; close &out._o2CovY;

	 %if &by ^=%str() %then create &out._PredVar from fVys [ colname = {"_rep" "&by" "&group" "&rtm" "_time" "Hat_y" "PredVar" "HatMSE" "Hat_y2" "PredVar2" "HatMSE2"} ];;;
	 %if &by  =%str() %then create &out._PredVar from fVys [ colname = {"_rep" "&group" "&rtm" "_time" "Hat_y" "PredVar" "HatMSE" "Hat_y2" "PredVar2" "HatMSE2"} ];;;
	 append from fVys; close &out._PredVar;

quit; /*eMKF: ends call to proc iml with matrix calculations */

/*eMKF: column names */
data &out._H;
  set &out._H;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4= &rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep 		   col2= &group col3= &rtm col4=_time;;;
run;
data &out._o2H;
  set &out._o2H;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4= &rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep 		   col2= &group col3= &rtm col4=_time;;;
run;

/*eMKF: column names */
data &out._CovY;
  set &out._CovY;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4= &rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep 		   col2= &group col3= &rtm col4=_time;;;
run;
data &out._o2CovY;
  set &out._o2CovY;
  %if &by ^=%str() %then rename col1=_rep col2=&by col3= &group col4= &rtm col5=_time;;;
  %if &by  =%str() %then rename col1=_rep 		   col2= &group col3= &rtm col4=_time;;;
run;

/*eMKF: merge with predictions dataset */
data &out._pred;
  merge &out._pred &out._PredVar;
  by _rep &group _time;
  drop &rtm.0 &rtm.1 &rtm.2 &rtm.3; /* eMKF: remove polynomial time terms from &out._pred */
run;

/*eMKF: BLUP estimation results */
data &out._BLUP;
  set &out._BLUP;
  _l2 = length(effect)-5;
  _time = 1*substr(effect,6,_l2);
  rename estimate=gamma_blup stderrpred=gamma_se;
  drop _l2;
run;

/*eMKF: sort and merge with predictions dataset */

proc sort data=&out._BLUP;
  by _rep _group_ _time;
run;

proc sort data=&out._pred;
  by _rep _group_ _time;
run;

data &out._pred; 
  merge &out._pred &out._BLUP(keep=_rep _group_ _time gamma_blup gamma_se);
  by  _rep _group_ _time;
  pred_blup = mu + gamma_blup ;
  pred_blup2 = mu2 + (moddelta*gamma_blup) ;
  label gamma_blup = " BLUP estimate of gamma from NLMIXED"
        gamma_se   = " Standard error of the BLUP estimate of gamma from NLMIXED"
        pred_blup  = "Outcome &outcome prediction using the BLUP estimate of gamma from NLMIXED"
        pred_blup2 = "Outcome &outcome2 prediction using the BLUP estimate of gamma from NLMIXED"
  ;
run;

/************************/
/* End of MSE estimation*/
/************************/

/*************************************************************************/
/* eMKF: Reverse-transform regression coefficients and covariance matrix */
/*************************************************************************/

data _tests1 _tests2 _lests _tfits _tcmat1 _tcmat2 _covmat1 _covmat2 _covmat12 _covmat21 
	 _covmat _covmatt _tcovmat _tcovmat2 _tcovmatt _tcovmatt2 _dcovmatt ;
run;

%if %upcase(&orpoly) = YES %then %do;

	/* eMKF: reverse-transform regression coefficients */
	%let k = %eval(&p - 1); %let jj = 0; 
	proc iml;

		use _oPmat_;
		read all into oPP; close _oPmat_;
		oPP = oPP[1:&p, 1:&p];

		/* Outcome 1 */
		varNames = {"_rep" "_group_" "o1a" };
		%if &p > 1 %then %do;
			o1bNames = "o1b1":"o1b&k";	
			varNames = varNames || o1bNames;
		%end;
		create _tests1 var varNames;
		%do jj=1 %to &nrep;
			use &out._ests(where=(_rep = &jj) keep= _rep _group_ o1a %if &p > 1 %then o1b: ; ) ;
			read all into oB; close &out._ests;
			oB1 = oB[,1:2];
            oB = T(oB[,3:(2+&p)]);
            oBB = oPP * oB;
			oBB = T(oBB);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _tests1;

		/* Outcome 2 */
		varNames = {"_rep" "_group_" "o2a" };
		%if &p > 1 %then %do;
			o1bNames = "o2b1":"o2b&k";	
			varNames = varNames || o1bNames;
		%end;
		create _tests2 var varNames;
		%do jj=1 %to &nrep;
			use &out._ests(where=(_rep = &jj) keep= _rep _group_ o2a %if &p > 1 %then o2b: ; ) ;
			read all into oB; close &out._ests;
			oB1 = oB[,1:2];
            oB = T(oB[,3:(2+&p)]);
            oBB = oPP * oB;
			oBB = T(oBB);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _tests2;

	quit;

	/* eMKF: sort by _rep and _group_ */
	proc sort data=_tests1;
  	  by _rep _group_ ;
	run;
	proc sort data=_tests2;
  	  by _rep _group_ ;
	run;

	/* eMKF: update estimates dataset */
	data &out._ests;
  	  merge &out._ests(drop=o1a o2a %if &p > 1 %then o1b: o2b: ;) _tests1 _tests2;
	  by _rep _group_;
	run;

	/* eMKF: update predictions dataset */
	data &out._pred;
  	  merge &out._pred(drop=o1a o2a %if &p > 1 %then o1b: o2b: ; mu err mu2 err2 _tausq2_ _rho2_ err3 _se3 
							delta lambda w gamma prediction delta2 lambda2 w2 gamma2 prediction2 
							group_rep gammaeff predeff predeff2 Hat_y PredVar HatMSE Hat_y2 PredVar2 HatMSE2 
							gamma_blup gamma_se pred_blup pred_blup2) 
			_tests1 _tests2
			&out._pred(keep=_rep _group_ mu err mu2 err2 _tausq2_ _rho2_ err3 _se3 
							delta lambda w gamma prediction delta2 lambda2 w2 gamma2 prediction2 
							group_rep gammaeff predeff predeff2 Hat_y PredVar HatMSE Hat_y2 PredVar2 HatMSE2 
							gamma_blup gamma_se pred_blup pred_blup2) 
	  		;
	  by _rep _group_;
	run;

	/* eMKF: symbolic set up for reverse-transformation of covariance matrix */
	%let parList1 = ; %let colList1 = ; %let parList2 = ; %let colList2 = ; %let oPPmat = ;
	%let i = 0;
	%do i=1 %to &g;  /* eMKF: block diagonal by group */
		%if &i = 1 %then %let oPPmat = block( oP ;
		%if &i > 1 and &i < &g %then %let oPPmat = &oPPmat , block ( oP ;
		%if &i = &g and &g > 1  %then %let oPPmat = &oPPmat , oP %sysfunc(repeat( %str(%)), &g-2));
		%if &i = &g and &g = 1  %then %let oPPmat = &oPPmat );
	%end;
	%let i = 0;
	%do i=1 %to &g; /* eMKF: column and row labels */
		%let parList1 = &parList1 %bquote(")o1a&i%bquote(");
		%let parList2 = &parList2 %bquote(")o2a&i%bquote(");
		%let colList1 = &colList1 o1a&i;
		%let colList2 = &colList2 o2a&i;
		%if %upcase(&bvalue) = INDEP_LINEAR or %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList1 = &parList1 %bquote(")o1b1arr&i%bquote(");
			%let parList2 = &parList2 %bquote(")o2b1arr&i%bquote(");
			%let colList1 = &colList1 o1b1arr&i;
			%let colList2 = &colList2 o2b1arr&i;
		%end;
		%if %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList1 = &parList1 %bquote(")o1b2arr&i%bquote(");
			%let parList2 = &parList2 %bquote(")o2b2arr&i%bquote(");
			%let colList1 = &colList1 o1b2arr&i;
			%let colList2 = &colList2 o2b2arr&i;
		%end;
		%if %upcase(&bvalue) = INDEP_CUBIC %then %do;
			%let parList1 = &parList1 %bquote(")o1b3arr&i%bquote(");
			%let parList2 = &parList2 %bquote(")o2b3arr&i%bquote(");
			%let colList1 = &colList1 o1b3arr&i;
			%let colList2 = &colList2 o2b3arr&i;
		%end;
	%end;
	%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList1 = &parList1 %bquote(")o1b1%bquote(");
		%let parList2 = &parList2 %bquote(")o2b1%bquote(");
		%let colList1 = &colList1 o1b1;
		%let colList2 = &colList2 o2b1;
	%end;
	%if %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList1 = &parList1 %bquote(")o1b2%bquote(");
		%let parList2 = &parList2 %bquote(")o2b2%bquote(");
		%let colList1 = &colList1 o1b2;
		%let colList2 = &colList2 o2b2;
	%end;
	%if %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let parList1 = &parList1 %bquote(")o1b3%bquote(");
		%let parList2 = &parList2 %bquote(")o2b3%bquote(");
		%let colList1 = &colList1 o1b3;
		%let colList2 = &colList2 o2b3;
	%end;

	%let parList1 = %unquote(&parList1);
	%let parList2 = %unquote(&parList2);

	/* eMKF: square block matrix of covariances for regression coefficients by group */
	data _covmat1; 	/* Outcome 1 */
	  retain _rep Row Parameter &colList1;;
	  set &out._covmat(where=(Parameter in(&parList1)) keep= _rep Row Parameter &colList1);;
	run;
	data _covmat2; /* Outcome 2 */
	  retain _rep Row Parameter &colList2;;
	  set &out._covmat(where=(Parameter in(&parList2)) keep= _rep Row Parameter &colList2);;
	run;
	data _covmat12; /* Outcome 1 by 2 */
	  retain _rep Row Parameter &colList2;;
	  set &out._covmat(where=(Parameter in(&parList1)) keep= _rep Row Parameter &colList2);;
	run;
	data _covmat21; /* Outcome 2 by 1 */
	  retain _rep Row Parameter &colList1;;
	  set &out._covmat(where=(Parameter in(&parList2)) keep= _rep Row Parameter &colList1);;
	run;

	/* eMKF: obtain new row numbers associated with modified column order */
	data _tcmat1; /* Outcome 1 */
	  set _covmat1;
      by _rep;
      if first._rep then output;
 	  drop Row Parameter;
	run;
	proc transpose data=_tcmat1 out=_tcmat1;
      by _rep;
    run;
	data _tcmat1;
	  set _tcmat1;
	  nRow + 1;
	  rename _NAME_ = Parameter;
	  drop col:;
	run;
	data _tcmat2; /* Outcome 2 */
	  set _covmat2;
      by _rep;
      if first._rep then output;
 	  drop Row Parameter;
	run;
	proc transpose data=_tcmat2 out=_tcmat2;
      by _rep;
    run;
	data _tcmat2;
	  set _tcmat2;
	  nRow + 1;
	  rename _NAME_ = Parameter;
	  drop col:;
	run;

	/* eMKF: Determine offset to add to row numbers in _tcmat2*/
	%let _tcobs = 0;
	data _null_;
	  if 0 then set _tcmat1 nobs=n;
	  call symput("_tcobs", n);
	  stop;
	run;
	%let _tcobs = %eval(0 + &_tcobs);
	data _tcmat2;
	  set _tcmat2;
	  nRow = nRow + &_tcobs;
	run;

	/* eMKF: sort, merge, and re-sort using the new row numbers */
	proc sort data=_covmat1; /* Outcome 1 */
	  by _rep Parameter;
	run;
	proc sort data=_covmat12;
	  by _rep Parameter;
	run;
	proc sort data=_tcmat1;
	  by _rep Parameter;
	run;
	data _covmat1;
	  merge _covmat1 _covmat12 _tcmat1;
	  by _rep Parameter;
	run;
	proc sort data=_covmat1;
	  by _rep nRow;
	run;
	proc sort data=_covmat2; /* Outcome 2 */
	  by _rep Parameter;
	run;
	proc sort data=_covmat21;
	  by _rep Parameter;
	run;
	proc sort data=_tcmat2;
	  by _rep Parameter;
	run;
	data _covmat2;
	  merge _covmat21 _covmat2 _tcmat2;
	  by _rep Parameter;
	run;
	proc sort data=_covmat2;
	  by _rep nRow;
	run;

	/*eMKF: set into single square block matrix */
	data _covmat;
	  set _covmat1 _covmat2;
	run;
	data _covmat;
	  merge _covmat(drop=nRow) _covmat(keep=nRow);
	run;

	/* eMKF: rectangular block matrix of covariances between regression coefficients and remaining parameters */
	data _covmatt;
	  retain _rep Row Parameter &colList1 &colList2;;
	  set &out._covmat(where=(Parameter not in(&parList1 &parlist2)) keep= _rep Row Parameter &colList1 &colList2);;
	run;

	/* eMKF: apply reverse-transformation to both square and rectangular block matrices */
	%let i = 0; %let jj = 0;
	proc iml;

		use _oPmat_;
		read all into oP; close _oPmat_;
		oP = oP[1:&p, 1:&p];
		oPP = &oPPmat;;

		/* eMKF: re-structure block matrix in the common trend cases (where &p > 1) */
		%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
			oPP1 = oPP[do(1, &g*&p, &p), do(1, &g*&p, &p)];
			oPP2 = vecdiag(oPP[do(1, &g*&p, &p), do(2, &g*&p, &p)]);
			ToPP2 = vecdiag(oPP[do(2, &g*&p, &p), do(1, &g*&p, &p)]);
			oPP1 = oPP1 // T(ToPP2);
			oPP0 = oPP2;
			%if &p > 2 %then %do;
				oPP3 = vecdiag(oPP[do(1, &g*&p, &p), do(3, &g*&p, &p)]);
				ToPP3 = vecdiag(oPP[do(3, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP3);
				oPP0 = oPP0 || oPP3;
			%end;
			%if &p > 3 %then %do;
				oPP4 = vecdiag(oPP[do(1, &g*&p, &p), do(4, &g*&p, &p)]);
				ToPP4 = vecdiag(oPP[do(4, &g*&p, &p), do(1, &g*&p, &p)]);
				oPP1 = oPP1 // T(ToPP4);
				oPP0 = oPP0 || oPP4;
			%end;
			oPP0 = oPP0 // oPP[2:&p, 2:&p];
			oPP = oPP1 || oPP0;
		%end;

		oPP = block(oPP, oPP); /* one block per outcome */

		varNames = {"_rep" "Row"};
		%do i=1 %to &g; /* Outcome 1 */
			varNames = varNames || {"o1a&i"};
			%if %upcase(&bvalue) = INDEP_LINEAR %then varNames = varNames || {"o1b1arr&i"};;
			%if %upcase(&bvalue) = INDEP_QUAD   %then varNames = varNames || {"o1b1arr&i"} || {"o1b2arr&i"};;
			%if %upcase(&bvalue) = INDEP_CUBIC  %then varNames = varNames || {"o1b1arr&i"} || {"o1b2arr&i"} || {"o1b3arr&i"};;
		%end;
		%if %upcase(&bvalue) = COMMON_LINEAR %then varNames = varNames || {"o1b1"};;
		%if %upcase(&bvalue) = COMMON_QUAD   %then varNames = varNames || {"o1b1" "o1b2"};;
		%if %upcase(&bvalue) = COMMON_CUBIC  %then varNames = varNames || {"o1b1" "o1b2" "o1b3"};;
		%do i=1 %to &g; /* Outcome 2 */
			varNames = varNames || {"o2a&i"};
			%if %upcase(&bvalue) = INDEP_LINEAR %then varNames = varNames || {"o2b1arr&i"};;
			%if %upcase(&bvalue) = INDEP_QUAD   %then varNames = varNames || {"o2b1arr&i"} || {"o2b2arr&i"};;
			%if %upcase(&bvalue) = INDEP_CUBIC  %then varNames = varNames || {"o2b1arr&i"} || {"o2b2arr&i"} || {"o2b3arr&i"};;
		%end;
		%if %upcase(&bvalue) = COMMON_LINEAR %then varNames = varNames || {"o2b1"};;
		%if %upcase(&bvalue) = COMMON_QUAD   %then varNames = varNames || {"o2b1" "o2b2"};;
		%if %upcase(&bvalue) = COMMON_CUBIC  %then varNames = varNames || {"o2b1" "o2b2" "o2b3"};;

		create _tcovmat var varNames;
		%do jj=1 %to &nrep;
			use _covmat(where=(_rep = &jj) drop= Parameter nRow);
			read all into oB; close _covmat;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
            oBB = oPP * oB * T(oPP);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _tcovmat;

		create _tcovmatt var varNames;
		%do jj=1 %to &nrep;
			use _covmatt(where=(_rep = &jj) drop= Parameter);
			read all into oB; close _covmatt;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
			oBB = oPP * T(oB);
			oBB = oB1 || T(oBB);
			append from oBB;
		%end;
		close _tcovmatt;

	quit;

	/* eMKF: combine both square and rectangular block matrices */
    data _tcovmat;
	  set _tcovmat _tcovmatt;
	run;

	/* eMKF: re-sort rows */
	proc sort data=_tcovmat;
  		by _rep Row ;
	run;

	/* eMKF: re-order columns as they were initially from NLMIXED */
	data _tcovmat2;
  	  retain  _rep Row 
			  o1a1-o1a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then o1b1arr1-o1b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g o1b3arr1-o1b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then o1b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then o1b1 o1b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then o1b1 o1b2 o1b3 ; 
			  o2a1-o2a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then o2b1arr1-o2b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g o2b3arr1-o2b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then o2b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then o2b1 o2b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then o2b1 o2b2 o2b3 ; 
	  ;
	  set _tcovmat; 
	run;
	data _tcovmatt2;
  	  retain  _rep Row 
			  o1a1-o1a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then o1b1arr1-o1b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then o1b1arr1-o1b1arr&g o1b2arr1-o1b2arr&g o1b3arr1-o1b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then o1b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then o1b1 o1b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then o1b1 o1b2 o1b3 ; 
			  o2a1-o2a&g 
			  %if %upcase(&bvalue) = INDEP_LINEAR  %then o2b1arr1-o2b1arr&g ; 
			  %if %upcase(&bvalue) = INDEP_QUAD    %then o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g ; 
			  %if %upcase(&bvalue) = INDEP_CUBIC   %then o2b1arr1-o2b1arr&g o2b2arr1-o2b2arr&g o2b3arr1-o2b3arr&g ; 
			  %if %upcase(&bvalue) = COMMON_LINEAR %then o2b1 ; 
			  %if %upcase(&bvalue) = COMMON_QUAD   %then o2b1 o2b2 ; 
			  %if %upcase(&bvalue) = COMMON_CUBIC  %then o2b1 o2b2 o2b3 ; 
	  ;
	  set _tcovmatt; 
	run;
	
	/* eMKF: merge */
	data _tcovmatt2;
  	  merge &out._covmat(where=(Parameter not in(&parList1 &parList2)) keep = _rep Row Parameter) _tcovmatt2;
	  by _rep Row;
	run;

	/* eMKF: transpose _tcovmatt2 to add into larger matrix */
	proc transpose data=_tcovmatt2(drop=Row) out=_tcovmatt2 name = Parameter;
      by _rep;
	  id Parameter;
    run;

	/* eMKF: insert Row numbers */
	data _tcovmatt2;
	  set _tcovmatt2;
	  by _rep;
	  retain Row;
	  if first._rep then Row = 1;
	  else Row + 1;
	run;

	/* eMKF: update covariance matrix dataset */
	data &out._covmat;
  	  merge &out._covmat(keep = _rep Row Parameter) 
			_tcovmat2
 			_tcovmatt2 
			&out._covmat(where=(Parameter not in(&parList1 &parList2)) drop= o1a: o2a: %if &p > 1 %then o1b: o2b: ;)
			;
	  by _rep Row;
	run;

	/* eMKF: extract variances of model parameters */
	%let jj = 0; 		
	proc iml;
		create _dcovmatt var{"_rep" "Row" "Var"};
		%do jj=1 %to &nrep;
			use &out._covmat(where=(_rep = &jj) drop= Parameter);
			read all into oB; close &out._covmat;
			oB1 = oB[,1:2];
			oB = oB[,3:ncol(oB)];
			oBB = vecdiag(oB);
			oBB = oB1 || oBB;
			append from oBB;
		%end;
		close _dcovmatt;
	quit;
	data _dcovmatt;
	  merge _dcovmatt &out._covmat(keep = _rep Row Parameter);
	  by _rep Row;
	run;

	/* eMKF: reset regression coefficients for o2 as differences relative to o1 */
	data _tests2;
  	  merge _tests2 _tests1;
	  by _rep _group_;
  	  o2a = o2a - o1a;
  	  %if &p > 1 %then %str(o2b1 = o2b1 - o1b1;);;
  	  %if &p > 2 %then %str(o2b2 = o2b2 - o1b2;);;
  	  %if &p > 3 %then %str(o2b3 = o2b3 - o1b3;);;
	run;

	/* eMKF: reverse-transformed regression estimates in long form */
	%if %upcase(&bvalue) ^= COMMON_LINEAR and %upcase(&bvalue) ^= COMMON_QUAD and %upcase(&bvalue) ^= COMMON_CUBIC %then %do;
		data _lests;
		  set _tests1(keep = _rep o1a rename=(o1a=Est))
			  %if &p > 1 %then _tests1(keep = _rep o1b1 rename=(o1b1=Est));
			  %if &p > 2 %then _tests1(keep = _rep o1b2 rename=(o1b2=Est));
			  %if &p > 3 %then _tests1(keep = _rep o1b3 rename=(o1b3=Est));
			  _tests2(keep = _rep o2a rename=(o2a=Est))
			  %if &p > 1 %then _tests2(keep = _rep o2b1 rename=(o2b1=Est));
			  %if &p > 2 %then _tests2(keep = _rep o2b2 rename=(o2b2=Est));
			  %if &p > 3 %then _tests2(keep = _rep o2b3 rename=(o2b3=Est));
		  ;
		  by _rep;
		  retain Row;
		  if first._rep then Row = 1;
		  else Row + 1;
		run;
	%end;
	%if %upcase(&bvalue) = COMMON_LINEAR or %upcase(&bvalue) = COMMON_QUAD or %upcase(&bvalue) = COMMON_CUBIC %then %do;
		%let jj = 0; 		
		data _lests;
		  set _tests1(keep = _rep o1a rename=(o1a=Est))
		  	  %do jj=1 %to &nrep;
			  	  %if &p > 1 %then _tests1(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o1b1 rename=(o1b1=Est));
			   	  %if &p > 2 %then _tests1(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o1b2 rename=(o1b2=Est));
			  	  %if &p > 3 %then _tests1(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o1b3 rename=(o1b3=Est));
			  %end;
			  _tests2(keep = _rep o2a rename=(o2a=Est))
		  	  %do jj=1 %to &nrep;
			  	  %if &p > 1 %then _tests2(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o2b1 rename=(o2b1=Est));
			   	  %if &p > 2 %then _tests2(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o2b2 rename=(o2b2=Est));
			  	  %if &p > 3 %then _tests2(where=(_rep=&jj) firstobs=1 obs=1 keep = _rep o2b3 rename=(o2b3=Est));
			  %end;
		  ;
		  by _rep;
		  retain Row;
		  if first._rep then Row = 1;
		  else Row + 1;
		run;
	%end;

	/* eMKF: merge reverse-transformed estimates and variances into fitstat dataset and update */
	data _tfits;
	  merge &out._fitstat _dcovmatt;
	  by _rep ;
	run;
	data _tfits;
	  merge _tfits _lests;
	  by _rep Row;
	run;
	data _tfits;
	  set _tfits;
	  if Est ne . then Estimate = Est;
	  if Var > 0 then StandardError = sqrt(Var);
	  else StandardError = .;
	  tValue = Estimate/StandardError;
	  Probt = (1-probt(abs(tValue), DF))*2;
	  Lower = Estimate + tinv(Alpha/2, DF)*StandardError;
	  Upper = Estimate + tinv(1-Alpha/2, DF)*StandardError;
	  drop Row Var Est;
	run;
	data &out._fitstat;
	  set _tfits;
	run;

%end;


/*eMKF: Add labels for stratification variable */
%if &by ^=%str() %then %do;
	data &out._ests;
	  merge &out._ests _freq_;
	  by _rep;
	run;
%end;

data &out._ests; /* eMKF: added parameter estimates for quadratic and cubic terms */
  set &out._ests;
  label _2loglike =" -2 log-likelihood estimate"
        _rho_ ="Estimated or supplied rho of the model"
	    _tausq_="Estimated or supplied tau-square of the model"
        delta ="Multiplicative factor of the latent AR(1) process for &outcome2"
	    _group_="Model reset group ID just in case group is not ordered"
	    &group ="Numeric &&group variable"
		/*eMKF: Added labels for by variable */
		_rep="Model reset stratum ID just in case stratification variable, if any, is not ordered"
		&by="&&by variable"
	    o1a="Parameter estimate (a) for the outcome &outcome"
        o2a="Parameter estimate (a) for the outcome &outcome2"
		%if %upcase(&bvalue) ^=DROPPED
			%then o1b1="Parameter estimate (b1) for the outcome &outcome assuming &bvalue trend model "
            	  o2b1="Parameter estimate (b1) for the outcome &outcome2 assuming &bvalue trend model " ;
  		%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC or 
			 %upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD 
				%then o1b2="Parameter estimate (b2) for the outcome &outcome assuming &bvalue trend model "
                	  o2b2="Parameter estimate (b2) for the outcome &outcome2 assuming &bvalue trend model " ;
  		%if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC 
			%then o1b3="Parameter estimate (b3) for the outcome &outcome assuming &bvalue trend model "
                  o2b3="Parameter estimate (b3) for the outcome &outcome2 assuming &bvalue trend model " ;
	;
	drop &group %if &by ^= %str() %then &by;
	;
run;

data &out._pred;  /*eMKF: added a few useful labels */
	 set &out._pred(rename=(Predvar=PredOnlyVar HatMSE=PredMSE Predvar2=PredOnlyVar2 HatMSE2=PredMSE2));
	 PredSE= sqrt(PredMSE);
	 PredSE2= sqrt(PredMSE2);
	 label _2loglike =" -2 log-likelihood estimate"
	       _rho_ ="Estimated or supplied rho of the model for the outcome &outcome"
		   _tausq_="Estimated or supplied tau-square of the model for the outcome &outcome"
	       _rho2_  ="Estimated or supplied rho of the model for the outcome &outcome2"
		   _tausq2_="Estimated or supplied tau-square of the model for the outcome &outcome2"
		   _group_="Model reset group ID just in case group is not ordered"
		   &group ="Numeric &&group variable"
		   _rep="Model reset stratum ID just in case stratification variable, if any, is not ordered"
		   &by="&&by variable"
		   &rlag = "Elapsed real time from previous time point"
		   &rtm = "Real time used in calculations"
		   _time ="Time "
		   _y  ="Original outcome for the outcome &outcome"
		   _se ="Original Standard Error for the outcome &outcome"
		   _avgse ="Average Standard Error across timepoints for the outcome &outcome used for imputation"
	       %if &by ^= %str() %then _avgseb ="Average Standard Error across strata for the outcome &outcome used for imputation";
           _y2  ="Original outcome for the outcome &outcome2"
		   _se2 ="Original Standard Error for the outcome &outcome2"
		   _avgse2 ="Average Standard Error across timepoints for the outcome &outcome2 used for imputation"
		   %if &by ^= %str() %then _avgse2b ="Average Standard Error across strata for the outcome &outcome2 used for imputation";
		   impute= "Whether original Standard Error for outcome &outcome or &outcome2 was imputed using average across timepoints"
		   %if &by ^= %str() %then imputeb= "Whether original Standard Error for outcome &outcome or &outcome2 was imputed using average across strata";
		   inputorder="Original orderng of the groups if it was not alphabetical"
		   prediction="Kalman estimator prediction of the outcome &outcome assuming &bvalue trend model"
		   PredMSE="Prediction variability: Mean Squared Error (MSE) for the outcome &outcome assuming &bvalue trend model "
		   PredSE="Prediction standard error: Square root of MSE for the outcome &outcome assuming &bvalue trend model "
		   predeff ="Kalman estimator prediction of the outcome &outcome assuming &bvalue trend model and using a combined (averaged) gamma"
		   prediction2="Kalman estimator prediction of the outcome &outcome2 assuming &bvalue trend model"
		   PredMSE2="Prediction variability: Mean Squared Error (MSE) for the outcome &outcome2 assuming &bvalue trend model "
		   PredSE2="Prediction standard error: Square root of MSE for the outcome &outcome2 assuming &bvalue trend model "
	       predeff2="Kalman estimator prediction of the outcome &outcome2 assuming &bvalue trend model and using a combined (averaged) gamma"
		   o1a="Parameter estimate (a) for the outcome &outcome"
	       o2a="Parameter estimate (a) for the outcome &outcome2"
	 	   %if %upcase(&bvalue) ^=DROPPED 
		   	   %then o1b1="Parameter estimate (b1) for the outcome &outcome assuming &bvalue trend model "
	                 o2b1="Parameter estimate (b1) for the outcome &outcome2 assuming &bvalue trend model " ;
	   	   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC or 
				%upcase(&bvalue) = INDEP_QUAD or %upcase(&bvalue) = COMMON_QUAD 
					%then o1b2="Parameter estimate (b2) for the outcome &outcome assuming &bvalue trend model "
	                 	  o2b2="Parameter estimate (b2) for the outcome &outcome2 assuming &bvalue trend model " ;
		   %if %upcase(&bvalue) = INDEP_CUBIC or %upcase(&bvalue) = COMMON_CUBIC
				%then o1b3="Parameter estimate (b3) for the outcome &outcome assuming &bvalue trend model "
	                  o2b3="Parameter estimate (b3) for the outcome &outcome2 assuming &bvalue trend model " ;
	 ;
	 /* These will be deleted for now. If needed they can be useful */
	 drop mu mu2 err err2 delta delta2 lambda lambda2 w w2 Hat_y PredOnlyVar Hat_y2 PredOnlyVar2 &group group_rep &rlag
  	      %if &by ^= %str() %then &by;
	; 
run;

/*eMKF: Rename any remaining instances of numeric &group variable to _group_ + remove numeric &by variable */

data &out._predVar;
  set &out._predVar;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;

data &out._H;
  set &out._H;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;
data &out._o2H;
  set &out._o2H;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;

data &out._covY;
  set &out._covY;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;
data &out._o2covY;
  set &out._o2covY;
  rename &group = _group_;
  %if &by ^= %str() %then drop &by;;
run;


proc datasets nolist;
 delete _freqn_ _freqg_ _freq_ _sdata_ _jdata_ _beta_ _inits_ _initsr_ _ests _fitstat _llike_ _empty_ 
        _Amat_ _Dmat_ _Vmat_ _Vgmat_ _Vemat_ _Xmat_ _junk_ _junk0_ _junk01_ aa _oXmat_ _oPmat_
        _tests1 _tests2 _lests _tfits _tcmat1 _tcmat2 _covmat1 _covmat2 _covmat12 _covmat21 
	 	_covmat _covmatt _tcovmat _tcovmat2 _tcovmatt _tcovmatt2 _dcovmatt 
        ;
run ;
quit; 

%mend;

data _null_;
run;


/* 
This macro allows us to reformat the data in a way that only one 
outcome is defined with time and group variables
*/

/* eMKF: Modified to track real time for use in calculations and calculate lags between succesive time points
 *		 Modified to allow for effective sample sizes to be specified, and to allow for imputing SEs across strata if needed
 * 		 Added random variance option
 */

%macro reformat(
			 data=, 
             outcome=, 
             se=,
             neff=,    			
			 outcome2=, 
             se2=, 
			 neff2=,
             group=, 
             time=, 
             by= ,
			 randomVars = NO, 	
			 outformat= 
             );

%local wt1 wt2 wn wi wa wb;

/* eMKF: inputorder variable definition */
data &outformat _freqg_ _freq_ _freqn_;
run;
data &outformat;
 set &data;
 inputorder +1;
run;

/* eMKF: _group_ variable definition */
proc freq data=&data noprint;
 tables &group /list out=_freqg_;
run;
data _freqg_;
 set _freqg_;
 _group_ +1;
 keep _group_ &group;
run;

/* eMKF: sort and merge */
proc sort data=&outformat;
 by &group;
run;
proc sort data=_freqg_;
 by &group;
run;
data &outformat;
 merge &outformat _freqg_;
 by &group;
run;

/* eMKF: _rep variable definition */

%if &by =%str() %then %do; /* eMKF: constant _rep */
	data &outformat;
	 set &outformat;
	 _rep=1;
	run;
%end;
%else %do; 	/* eMKF: incremented _rep */
	proc freq data=&outformat noprint;
	 tables &by /list out=_freq_;
	run;
	data _freq_;
	 set _freq_;
	 _rep +1;
	 keep _rep &by;
	run;
	proc sort data=&outformat; /* eMKF: sort and merge */
	 by &by;
	run;
	data &outformat;
	 merge &outformat _freq_;
	 by &by;
	run;
%end;

%let wn=; %let wt1=1; %let wt2=;
%let wi=0; %let wa=; %let wb=;

/* eMKF: If format 2, then find number of time points and define generic variables for outcome(s) and SE(s) */
%if %scan(&outcome,2) = %str() %then %do;

	/* eMKF: _time variable definition */
	proc freq data=&outformat noprint;
	 tables &time /list out=_freqn_;
	run;

	/*eMKF: modified to also calculate real time assuming &time is numeric: this allows handling unequally-spaced time points */
	data _freqn_;
	  set _freqn_;
	  retain base _rtimeold;
	  if _n_ = 1 then do;
		base=&time;
		_rtimeold=base;
	  end;
	  _rlag = &time - _rtimeold;
	  _rtimeold = &time;
	  _rtime = &time - base + 1;
	  _time + 1;
	  call symput("wt2", _time);	
	  keep _rlag _rtime _time &time;
	run;

	data _freqn_; /*eMKF: Set first lag to missing instead of 0*/
	  set _freqn_;
	  if _n_ = 1 then _rlag = .;
	run;

	/* eMKF: number of time points */
	%let wt2 = %eval(0 + &wt2 );

	/* eMKF: sort prior to merging */
	proc sort data=&outformat;
	 by &time;
	run;

	/* eMKF: merge and define outcomes(s) and SE(s) generic variables */
	data &outformat;
	 merge &outformat _freqn_;
	 by &time;
	 _y =&outcome;
	 _se=&se;
	 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _y2 =&outcome2;;
	 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _se2 =&se2;;
	run;

%end;
%else %do; /* Data is in format 1: need to reshape it into format 2 */

	%if %scan(&outcome,2,'-') =%str() %then %let wt2=%_counts_(&outcome);
	%if %scan(&outcome,2,'-') ^=%str() %then %do;
		%let wn=%length(%scan(&outcome,2,'-')); /*eMKF: number of time points */
		%do wi=1 %to &wn;
			%let wa=%scan(&outcome,1,'-');
			%let wb=%scan(&outcome,2,'-');
			%if %substr(&wa,1,&wi) ^= %substr(&wb,1,&wi) and &wt2=%str() %then %do ;
	           %let wt1 = %eval(0 + %substr(&wa,&wi, 1 + &wn - &wi) );
			   %let wt2 = %eval(0 + %substr(&wb,&wi, 1 + &wn - &wi) );
			%end;
		%end;
 	%end;

	/* eMKF: Define _time _y _se _y2 _se2 columns */
	data &outformat;
	 set &outformat;
	 array Ays(&wt1.:&wt2) &outcome;
	 array Ases(&wt1.:&wt2) &se;
	 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then array Ays2(&wt1.:&wt2) &outcome2;;
	 %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then array Ases2(&wt1.:&wt2) &se2;;
	 do _time = &wt1 to &wt2;
	    _y = Ays(_time);
		_se = Ases(_time);
		%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _y2 = Ays2(_time);;
		%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _se2 = Ases2(_time);;
		output;
	 end;
	 _rtime = _time; /*eMKF: Added for consistency but will just be assumed equally-spaced if in format 1 */
	 _rlag = 1;
	 keep _y _se &group _group_ _rep &by _time _rtime _rlag inputorder 
		  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _y2 _se2 ;
	 ;
	run;

%end; /* eMKF: end reshape format 1 into format 2 */

/*eMKF: mean imputations for missing or zero SEs */

data _means_ _means2_ _meansb_ _means2b_;
run;

data &outformat;
  set &outformat;
  if _se le 0 then _se=.;
  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _se2 le 0 then _se2=.;;;
run;

proc sort data= &outformat;
  by _rep _group_ _time;
run;

/* eMKF: Averages across timepoints */
proc means data=&outformat noprint;
  var _se;
  where _se ne .;
  by _rep _group_;
  output out= _means_(drop= _type_ _freq_) mean=_avgse;
run;
%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	proc means data=&outformat noprint;
	 var _se2;
	 where _se2 ne .;
	 by _rep _group_;
	 output out= _means2_(drop= _type_ _freq_) mean=_avgse2;
	run;
%end;

/* eMKF: Averages across strata */
%if &by ^= %str() %then %do; 
  proc sort data= &outformat;
    by _group_ _time _rep;
  run;
  proc means data=&outformat noprint;
    var _se;
    where _se ne .;
    by _group_ _time;
    output out= _meansb_(drop= _type_ _freq_) mean=_avgseb;
  run;
  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
    proc means data=&outformat noprint;
	  var _se2;
	  where _se2 ne .;
	  by _group_ _time;
	  output out= _means2b_(drop= _type_ _freq_) mean=_avgse2b;
	run;
  %end;
  proc sort data= &outformat;
    by _rep _group_ _time;
  run;
%end;

data &outformat;
  merge &outformat _means_ %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _means2_ ;;
  by _rep _group_;
  impute=0;
  if _se=. and _avgse > 0 then do;
	impute=1;
	_se = _avgse;
  end;
  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	if _se2=. and _avgse2 > 0 then do;
		impute=1;
		_se2 = _avgse2;
	end;
  %end;
run;

/*eMKF: (new in eMKF) use average SE across strata if still have missing SEs */
%if &by ^= %str() %then %do; 
  proc sort data= &outformat;
    by _group_ _time _rep;
  run;
  data &outformat;
    merge &outformat _meansb_ %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _means2b_ ;;
    by _group_ _time;
    imputeb=0;
    if _se=. and _avgseb > 0 then do;
	  imputeb=1;
	  _se = _avgseb;
    end;
    %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	  if _se2=. and _avgse2b > 0 then do;
		  imputeb=1;
		  _se2 = _avgse2b;
	  end;
    %end;
  run;
  proc sort data= &outformat;
    by _rep _group_ _time;
  run;
%end;

/*eMKF: Added effective sample size imputations */

data _means_ _means2_ _meansb_ _means2b_;
run;

%if %upcase(&randomVars) = YES %then %do;

	%if &neff = %str() %then %do;
		%put ERROR: (Effective) sample sizes neff must be specified to fit random sampling variances;
		proc iml;
			print "  Error Note:";
			print "  (Effective) sample sizes neff must be specified to fit random sampling variances. ";
		quit;
		%return;
	%end;

	data &outformat;
	  set &outformat;
	  _n = &neff;
	  if _n le 0 then _n=.;
	  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
	  	_n2 = &neff2;
		if _n2 le 0 then _n2=.;
	  %end;
	run;

	/* eMKF: Averages across timepoints */
	proc means data=&outformat noprint;
	  var _n;
	  where _n ne .;
	  by _rep _group_;
	  output out= _means_(drop= _type_ _freq_) mean=_avgn;
	run;
	%if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
		proc means data=&outformat noprint;
		  var _n2;
		  where _n2 ne .;
		  by _rep _group_;
		  output out= _means2_(drop= _type_ _freq_) mean=_avgn2;
		run;
	%end;

	/* eMKF: Averages across strata */
    %if &by ^= %str() %then %do; 
      proc sort data= &outformat;
        by _group_ _time _rep;
      run;
      proc means data=&outformat noprint;
        var _n;
        where _n ne .;
        by _group_ _time;
        output out= _meansb_(drop= _type_ _freq_) mean=_avgnb;
      run;
      %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then %do;
        proc means data=&outformat noprint;
	      var _n2;
	      where _n2 ne .;
	      by _group_ _time;
	      output out= _means2b_(drop= _type_ _freq_) mean=_avgn2b;
	    run;
      %end;
      proc sort data= &outformat;
        by _rep _group_ _time;
      run;
    %end;

	/*eMKF: first pass at imputation using average over timepoints */
	data &outformat;
	  merge &outformat _means_ %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _means2_ ;;
	  by _rep _group_;
	  if _n=. and _avgn > 0 then _n = _avgn;
	  %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _n2=. and _avgn2 > 0 then _n2 = _avgn2;;;
	run;

	/*eMKF: second pass at imputation using average across strata */
	%if &by ^= %str() %then %do; 
      proc sort data= &outformat;
        by _group_ _time _rep;
      run;
	  data &outformat;
	    merge &outformat _meansb_ %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then _means2b_ ;;
	    by _group_ _time;
	    if _n=. and _avgnb > 0 then _n = _avgnb;
	    %if %scan(&outcome2,1) ^=%str() and %scan(&se2,1) ^=%str() %then if _n2=. and _avgn2b > 0 then _n2 = _avgn2b;;;
	  run;
      proc sort data= &outformat;
        by _rep _group_ _time;
      run;
	%end;

%end;

data _means_;
run;
data _means_;
  set &outformat;
  %if &by = %str() %then if impute=1;;
  %if &by ^= %str() %then if impute=1 or imputeb=1;;
run;

data &outformat;
  merge &outformat(drop=impute %if &by ^= %str() %then imputeb;) 
        _means_(keep= _rep _group_ impute %if &by ^= %str() %then imputeb;);
  by _rep _group_;
  if impute=. then impute=0;
  %if &by ^= %str() %then if imputeb=. then imputeb=0;;
run;

data _means_;
  set &outformat;
  if _time= &wt2;
  order=inputorder;
  keep order _rep _group_ _time;
run;
proc sort data=_means_;
  by order;
run;
data _means_;
  set _means_;
  inputorder +1;
  drop order;
run;
proc sort data= _means_;
  by _rep _group_;
run;

proc sort data= &outformat;
  by _rep _group_;
run;
data &outformat;
  merge &outformat(drop=inputorder) _means_(keep= _rep _group_ inputorder);
  by _rep _group_;
run;
proc sort data= &outformat;
  by _rep inputorder;
run;

proc datasets nolist;
  delete _means_ _means2_ _meansb_ _means2b_ _freqg_ _freqn_ _freq_ ;
run ;
quit; 

%mend;

data _null_;
run;




/************************************************************************************/
/*                                 Utility macros                                   */
/************************************************************************************/
data _null_;
run;

/* 
Macro defined 03-02-2007
This macro %_counts_ allows to count number of words in a string.
The first argument is the string
split allows to specify the character that splits between the words
If split is empty it will just estimate the number of different words.

This will be used to count the number of variables in a list of variables.
This is just a utility macro
*/

%macro _counts_(arg1, text=&arg1, notes=, split=%str( ));

%local i;
%let i=0;

%do %while(%length(%nrbquote(%scan(&text, &i+1, &split))));
	%let i=%eval(&i+1);
%end;

&i

%if %length(&notes) %then %put NOTE: _COUNT is returning the value: &i..;
%mend _counts_;

/* 
Macro defined 03-02-2007
These 2 %zeros and %zeross macros just prints a number of zeros.
They respectively print them with and without comma.
 This is used in the non-linear mixed model to define parameter mean values.
*/

%macro zeros(t);
0
%do j = 1 %to %eval(&t-1);
 ,0
%end;
%mend; 

%macro zeross(t);
0
%do j = 1 %to %eval(&t-1);
  0
%end;
%mend; 

/* eMKF: Just like %zeross, this macro prints the specified constant a number of times */
%macro cnstss(v, t);
%sysevalf(&v)
%do j = 1 %to %eval(&t-1);
  %sysevalf(&v)
%end;
%mend; 

/* eMKF: This is the autoregressive covariance matrix from RAND used for random effects in proc nlmixed */
/* Compare to covariance matrix used in the built-in MVNAR distribution in proc mcmc */
%macro thevarcomp(t=, nu=, vrho= , delta= );
	%local ii_ jj_ kk_ ll_  _nn _tvar1 _tvar2 rr_ rr2_ nu_;
	%let _tvar1=;
	%let _tvar2=;
	%let ii_=0; %let jj_=0; %let kk_=0; %let ll_=0; /*eMKF: added to initialize values */
	%let _nn=0;
	%let _nn= %_counts_(&delta);
	%let _nn= %eval(&_nn + 1);

	[
	%do ii_ = 1 %to &_nn;
		%if &ii_=1 %then %do;
	  		%do jj_ = 1 %to &t;
   				%if &jj_ > 1 %then %do;
    				%do kk_ = 1 %to %eval(&jj_ -1);
     					%if %eval(&jj_ -&kk_ ) ne 1 %then ,&vrho.**%eval(&jj_ -&kk_ )*&nu ;
     					%if %eval(&jj_ -&kk_ ) =  1 %then ,&vrho.*&nu ;
    				%end; 
   				%end; 
   				%if &jj_ > 1 %then , ;
   				&nu
  			%end; 
 		%end;
	 	%if &ii_ > 1 %then %do;
	    	%let _tvar1=%scan( &delta, %eval(&ii_-1));
	     	%let nu_=;
	  		%do jj_ = 1 %to &t;
	  			%do ll_ =1 %to &ii_;
	   				%if &ll_ > 1 %then  %let _tvar2=%scan(&delta, %eval(&ll_-1));
	   				%if &ll_ < &ii_ %then %do;
	   					%let nu_=;
	   					%if &ll_ = 1 %then %let nu_ = &_tvar1.*&nu ;
	   					%if &ll_ > 1 %then %let nu_ = &_tvar2.*&_tvar1.*&nu ;
	    				%do kk_ = 1 %to &t;
	    					%let rr_= %eval(&jj_ -&kk_ );
	    					%let rr2_= %eval(&kk_ -&jj_ );
	      					%if &rr_ =  0  %then %do; ,&nu_ %end;
	      					%if &rr_ < 0 and &rr2_ ne 1 %then %do; ,&vrho.**&rr2_.*&nu_ %end;
	      					%if &rr_ > 0 and &rr_ ne 1 %then %do; ,&vrho.**&rr_.*&nu_ %end; 
	      					%if &rr_ < 0 and &rr2_ = 1 %then %do; ,&vrho.*&nu_ %end;
	      					%if &rr_ > 0 and &rr_  = 1 %then %do; ,&vrho.*&nu_ %end; 
	    				%end; 
	    			%end; 
	      			%if &ll_ = &ii_ %then %do;
	      				%let nu_=;
	      				%let nu_ = &_tvar1.**2*&nu ;
	   					%if &jj_ > 1 %then %do;
	    					%do kk_ = 1 %to %eval(&jj_ -1);
	     						%if %eval(&jj_ -&kk_ ) ne 1 %then ,&vrho.**%eval(&jj_ -&kk_ )*&nu_ ;
	     						%if %eval(&jj_ -&kk_ ) =  1 %then ,&vrho.*&nu_ ;
	    					%end; 
	   					%end; 
	   					,&nu_
	    			%end; 
	   			%end; 
  			%end; 
 		%end; 
	%end;
	]

%mend;

/* eMKF: This is a modified/simplified version of original thevarcomp macro to allow for unequally spaced or noninteger times */
/* eMKF: Also corrected how delta = 0 was handled, which should have been equivalent to delta = %str() */

%macro thevarcompr(times= , nu=, vrho= );
	%local jj_ kk_ tj_ tk_ t;
	%let jj_ = 0; %let kk_ = 0; 
	%let tj_ = 0; %let tk_ = 0; 
	%let t=0;
	%let t = %_counts_(&times);
	[
  	%do jj_ = 1 %to &t;
	    %let tj_ = %scan(&times, &jj_, character=" ");
   		%if &jj_ > 1 %then %do;
    		%do kk_ = 1 %to %eval(&jj_ - 1);
			   	%let tk_ = %scan(&times, &kk_, character=" ");
     			, &vrho.**%sysevalf(&tj_ - &tk_)*&nu 
    		%end; 
   		%end; 
   		%if &jj_ > 1 %then , ;
   		&nu
  	%end; 
	]
%mend;

data _null_;
run;





/****************************************************************************
 * eMKF
 * Wrapper macros to compile applicable Gibbs samplers for use with proc mcmc
 * Makram Talih, Ph.D. (NCHS contractor)
 ****************************************************************************/
data _null_;
run;

/**************************************************/
/* eMKF: Gibbs sampler for true state predictions */
/**************************************************/
%macro gibbs_uds_compile_EP(g=, n=, loc=);

%local uloc;
%let uloc = &loc..uds;
 
proc fcmp outlib=&uloc; 			

	subroutine EP(etaarr[*],			/* 1-dimensional array (length gn) of updated values of true states */
				  etamnarr[*],			/* 1-dimensional array (length gn) of predictions from regression */
				  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
				  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
				  rts[*],				/* 1-dimensional array (length n) of real times */
				  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
				  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
				  );

	outargs etaarr;						/* arguments that are updated after execution */

	array Yvec[&n, 1]					/nosym; /* vector (nx1) for use in calculations */
	array Xbeta[&n, 1]					/nosym;	/* holds matrix multiplication */
	array Wgamma[&n, &n]   				/nosym;	/* Vgamma^{-1} */
	array Vg[&n, &n]  					/nosym;	/* Vgamma + sampling variances */
	array Wg[&n, &n]   					/nosym;	/* (Vgamma + sampling variances)^{-1} */
	array petas[&n, 1]					/nosym;	/* contribution to posterior mean from random effects prior */
	array yetas[&n, 1]					/nosym;	/* contribution to posterior mean from data*/
	array etas[&n, 1]					/nosym;	/* sampled vector (nx1) of group-specific etas */
	array EE[&n, &n] 					/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
	array EI[&n, &n] 					/nosym;	/* inverse of EE */

	do k = 1 to &g;
		do i = 1 to &n;					/* predictions from regression for group k */
			Xbeta[i,1] = etamnarr[(k-1)*&n+i];	
		end;
	    call zeromatrix(Wgamma);												/* Wgamma = Vgamma^{-1} is tridiagonal */
		Wgamma[1,1] = (1/((1-(rhoarr[k]**(2*(rts[2]-rts[1]))))*nuarr[k]));		/* First main diagonal entry */
		do j = 2 to &n-1; 														/* Intermediate main diagonal entries */
			Wgamma[j,j] = ((1-(rhoarr[k]**(2*(rts[j+1]-rts[j-1]))))
				 	        /((1-(rhoarr[k]**(2*(rts[j+1]-rts[j]))))
				  		  	  *(1-(rhoarr[k]**(2*(rts[j]-rts[j-1]))))*nuarr[k])); 
  		end; 
		Wgamma[&n,&n]= (1/((1-(rhoarr[k]**(2*(rts[&n]-rts[&n-1]))))*nuarr[k])); /* Last main diagonal entry */
		do i = 1 to &n-1;														/* Entries in second and third diagonals */
			Wgamma[i,i+1] = -(rhoarr[k]**(rts[i+1]-rts[i]))
						 	/((1-(rhoarr[k]**(2*(rts[i+1]-rts[i]))))*nuarr[k]); 
			Wgamma[i+1,i] = Wgamma[i,i+1]; 
  		end; 
		call mult(Wgamma, Xbeta, petas);	/* contribution to posterior mean from random effects prior */
		call zeromatrix(Vg); 				/* initialize Vg to all zeroes */
		do j = 1 to &n;							
			Yvec[j,1] = Yarr[(k-1)*&n + j]; /* populate nx1 data vector Yvec */
			Vg[j,j] = 1/Sarr[(k-1)*&n + j]; /* populate diagonal entries of Vg using sampling precisions */
		end;
		call mult(Vg, Yvec, yetas);			/* contribution to posterior mean from data */
		call addmatrix(yetas, petas,petas); /* sum of prior and data contributions */
		call addmatrix(Wgamma, Vg, Wg);		/* precision matrix Wg = sum of Wgamma and diagonal of sampling precisions */
		do i = 1 to &n;
			etas[i,1] = rand('normal');		/* sample from univariate standard normal */
		end;
		call chol(Wg, EE);					/* calculate Cholesky decomposition for precision matrix (returns lower triangular) */
		call inv(EE, EI);					/* inverse of lower triangular matrix from Cholesky decomposition */
		call mult(EI, petas, petas);		/* re-scale petas (part 1) */
		call transpose(EI, EI);				/* transpose */
		call mult(EI, petas, petas);		/* re-scale petas (part 2) */
		call mult(EI, etas, etas);			/* re-scale etas */
		call addmatrix(petas, etas, etas);	/* re-center */
		do j = 1 to &n;
			etaarr[(k-1)*&n+j] = etas[j,1];	/* output argument etaarr is an array of group by time state predictions */
		end;
	end;

	endsub;
run;
quit;

%mend;

data _null_;
run;

/*****************************************************/
/* eMKF: Gibbs sampler for random sampling variances */
/*****************************************************/
%macro gibbs_uds_compile_RP(g=, n=, loc=);

%local uloc;
%let uloc = &loc..uds;
 
proc fcmp outlib=&uloc; 			

	subroutine RP(varr[*],				/* 1-dimensional array (length g) of updated values of true variances */
				  vhyp[*],				/* shape and scale for inverse gamma prior */
				  Sarr[*],				/* 1-dimensional array (length gn) for _var from dataset */
				  Narr[*]				/* 1-dimensional array for _n from dataset */
				  );

	outargs varr;						/* arguments that are updated after execution */

	do k = 1 to &g;	
		pvshape = vhyp[1];				/* shape parameter for inverse gamma prior */
		pvscale = vhyp[2];				/* scale parameter for inverse gamma prior */
		do j = 1 to &n;
			pvshape = pvshape + ((Narr[(k-1)*&n + j] - 1) / 2);
			pvscale = pvscale + ((Narr[(k-1)*&n + j] - 1) * Sarr[(k-1)*&n + j] / 2);
		end;
		varr[k] = 1/rand('gamma', pvshape, 1/pvscale);
	end;

	endsub;
run;
quit;

%mend;

data _null_;
run;

/*************************************************************************************/
/* eMKF: Gibbs samplers for mean hyper-parameters in the fully Bayesian trend models */
/*************************************************************************************/
%macro gibbs_uds_compile_MP(uvar=, g=, loc=);

/* eMKF: return if no applicable model is indicated */
%if %upcase(&uvar) ^= FULL_LINEAR and %upcase(&uvar) ^= FULL_QUAD  and %upcase(&uvar) ^= FULL_CUBIC  %then %do;
		%put ERROR: No Gibbs sampler for mean hyper-parameters was found for the specified model &uvar: Please check!;
		proc iml;
		  print "  Error Note:";
		  print "  No Gibbs sampler for mean hyper-parameters was found for the specified model. Please check. ";
	    quit;
	%return;
%end;

%local uloc;
%let uloc = &loc..uds;

%if %upcase(&uvar) = FULL_LINEAR %then %do;
	proc fcmp outlib=&uloc; 			

		subroutine MP_bfl(mb1,					/* updated value of prior mean mb1 */
						  mbetag[*,*], 			/* updated prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* updated diagonal matrix (p x p) of prior precisions for regression coefficients */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  mb1hyp[*],			/* hyper prior mean [1] and precision [2] for mb1 */
						  sb1					/* updated value of prior SD sb1 from M-H sampler */
						  );

		outargs mb1, mbetag, Dbetag;			/* arguments that are updated after execution */

		/*************************/
		/* Update prior mean mb1 */
		/*************************/
		b1mn = 0;
		do k = 1 to &g;
			b1mn = b1mn + b1arr[k];
		end;
		mb1mn = ((b1mn/mb1hyp[2]) + (mb1hyp[1]*(sb1**2)))/((&g/mb1hyp[2]) + (sb1**2));
		mb1sd = 1/sqrt(mb1hyp[2] + (&g/(sb1**2)));
		mb1 = mb1mn + mb1sd*rand('normal');

		mbetag[2,1] = mb1; 
		Dbetag[2,2] = 1/sb1**2;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = FULL_QUAD %then %do;
	proc fcmp outlib=&uloc; 			

		subroutine MP_bfq(mb1,					/* updated value of prior mean mb1 */
						  mb2,					/* updated value of prior mean mb2 */
						  mbetag[*,*], 			/* updated prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* updated diagonal matrix (p x p) of prior precisions for regression coefficients */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						  mb1hyp[*],			/* hyper prior mean [1] and precision [2] for mb1 */
						  mb2hyp[*],			/* hyper prior mean [1] and precision [2] for mb2 */
						  sb1,					/* updated value of prior SD sb1 from M-H sampler */
						  sb2					/* updated value of prior SD sb2 from M-H sampler */
						  );

		outargs mb1, mb2, mbetag, Dbetag;		/* arguments that are updated after execution */

		/*************************/
		/* Update prior mean mb1 */
		/*************************/
		b1mn = 0;
		do k = 1 to &g;
			b1mn = b1mn + b1arr[k];
		end;
		mb1mn = ((b1mn/mb1hyp[2]) + (mb1hyp[1]*(sb1**2)))/((&g/mb1hyp[2]) + (sb1**2));
		mb1sd = 1/sqrt(mb1hyp[2] + (&g/(sb1**2)));
		mb1 = mb1mn + mb1sd*rand('normal');

		mbetag[2,1] = mb1; 
		Dbetag[2,2] = 1/sb1**2;

		/*************************/
		/* Update prior mean mb2 */
		/*************************/
		b2mn = 0;
		do k = 1 to &g;
			b2mn = b2mn + b2arr[k];
		end;
		mb2mn = ((b2mn/mb2hyp[2]) + (mb2hyp[1]*(sb2**2)))/((&g/mb2hyp[2]) + (sb2**2));
		mb2sd = 1/sqrt(mb2hyp[2] + (&g/(sb2**2)));
		mb2 = mb2mn + mb2sd*rand('normal');

		mbetag[3,1] = mb2; 
		Dbetag[3,3] = 1/sb2**2;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = FULL_CUBIC %then %do;
	proc fcmp outlib=&uloc; 			

		subroutine MP_bfc(mb1,					/* updated value of prior mean mb1 */
						  mb2,					/* updated value of prior mean mb2 */
						  mb3,					/* updated value of prior mean mb3 */
						  mbetag[*,*], 			/* updated prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* updated diagonal matrix (p x p) of prior precisions for regression coefficients */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						  b3arr[*], 			/* 1-dimensional array (length g) of updated values of cubic coefficients by group */
						  mb1hyp[*],			/* hyper prior mean [1] and precision [2] for mb1 */
						  mb2hyp[*],			/* hyper prior mean [1] and precision [2] for mb2 */
						  mb3hyp[*],			/* hyper prior mean [1] and precision [2] for mb3 */
						  sb1,					/* updated value of prior SD sb1 from M-H sampler */
						  sb2,					/* updated value of prior SD sb2 from M-H sampler */
						  sb3					/* updated value of prior SD sb3 from M-H sampler */
						  );

		outargs mb1, mb2, mb3, mbetag, Dbetag;	/* arguments that are updated after execution */

		/*************************/
		/* Update prior mean mb1 */
		/*************************/
		b1mn = 0;
		do k = 1 to &g;
			b1mn = b1mn + b1arr[k];
		end;
		mb1mn = ((b1mn/mb1hyp[2]) + (mb1hyp[1]*(sb1**2)))/((&g/mb1hyp[2]) + (sb1**2));
		mb1sd = 1/sqrt(mb1hyp[2] + (&g/(sb1**2)));
		mb1 = mb1mn + mb1sd*rand('normal');

		mbetag[2,1] = mb1; 
		Dbetag[2,2] = 1/sb1**2;

		/*************************/
		/* Update prior mean mb2 */
		/*************************/
		b2mn = 0;
		do k = 1 to &g;
			b2mn = b2mn + b2arr[k];
		end;
		mb2mn = ((b2mn/mb2hyp[2]) + (mb2hyp[1]*(sb2**2)))/((&g/mb2hyp[2]) + (sb2**2));
		mb2sd = 1/sqrt(mb2hyp[2] + (&g/(sb2**2)));
		mb2 = mb2mn + mb2sd*rand('normal');

		mbetag[3,1] = mb2; 
		Dbetag[3,3] = 1/sb2**2;

		/*************************/
		/* Update prior mean mb3 */
		/*************************/
		b3mn = 0;
		do k = 1 to &g;
			b3mn = b3mn + b3arr[k];
		end;
		mb3mn = ((b3mn/mb3hyp[2]) + (mb3hyp[1]*(sb3**2)))/((&g/mb3hyp[2]) + (sb3**2));
		mb3sd = 1/sqrt(mb3hyp[2] + (&g/(sb3**2)));
		mb3 = mb3mn + mb3sd*rand('normal');

		mbetag[4,1] = mb3; 
		Dbetag[4,4] = 1/sb3**2;

		endsub;
	run;
	quit;

%end;

%mend;

data _null_;
run;

/**********************************************************************************/
/* eMKF: Gibbs samplers for regression coefficients in the supported trend models */
/**********************************************************************************/
%macro gibbs_uds_compile_CP(uvar=, g=, n=, loc=);

/* eMKF: return if no applicable model is indicated */
%if %upcase(&uvar) ^= DROPPED and 
	%upcase(&uvar) ^= FULL_LINEAR and %upcase(&uvar) ^= INDEP_LINEAR and %upcase(&uvar) ^= COMMON_LINEAR and
	%upcase(&uvar) ^= FULL_QUAD   and %upcase(&uvar) ^= INDEP_QUAD   and %upcase(&uvar) ^= COMMON_QUAD and 
	%upcase(&uvar) ^= FULL_CUBIC  and %upcase(&uvar) ^= INDEP_CUBIC  and %upcase(&uvar) ^= COMMON_CUBIC and
    %upcase(&uvar) ^= BMA_CUBIC   and %upcase(&uvar) ^= BMA_QUAD     and %upcase(&uvar) ^= BMA_LINEAR
	%then %do;
		%put ERROR: No Gibbs sampler was found for the specified model &uvar: Please check!;
		proc iml;
		  print "  Error Note:";
		  print "  No Gibbs sampler was found for the specified model. Please check. ";
	    quit;
	%return;
%end;

%local p uloc;

/* eMKF: dimensionality (needed for UDS set up) */
%let p = 0;
%if %upcase(&uvar) = DROPPED %then %let p = 1;
%if %upcase(&uvar) = BMA_LINEAR  or %upcase(&uvar) = FULL_LINEAR or %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = COMMON_LINEAR %then %let p = 2;
%if %upcase(&uvar) = BMA_QUAD    or %upcase(&uvar) = FULL_QUAD   or %upcase(&uvar) = INDEP_QUAD   or %upcase(&uvar) = COMMON_QUAD   %then %let p = 3;
%if %upcase(&uvar) = BMA_CUBIC   or %upcase(&uvar) = FULL_CUBIC  or %upcase(&uvar) = INDEP_CUBIC  or %upcase(&uvar) = COMMON_CUBIC  %then %let p = 4;
%let p = %eval(0+&p);

%let uloc = &loc..uds;

%if %upcase(&uvar) = DROPPED %then %do;
	/***********************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the dropped (no trend) model */
	/***********************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine  CP_b0(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  mbetag[*,*], 			/* prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  X[*,*], 				/* design matrix (n x p) using real times */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, etamnarr;					/* arguments that are updated after execution */

		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array Xt[&p, &n]   						/nosym;	/* transpose of design matrix */
		array XtW[&p, &n]					  	/nosym; /* matrix multiplication of Xt and Wg */
		array XtWX[&p, &p] 						/nosym; /* precision matrix of WLS regression estimators */
		array DXtWX[&p, &p]					 	/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array prbeta[&p, 1] 	    			/nosym;	/* vector (p x 1) of regression estimates from prior */
		array pbeta[&p, 1] 	    			   	/nosym;	/* vector (p x 1) of regression estimates from pooled posterior */
		array ybeta[&p, 1] 	       				/nosym;	/* vector (p x 1) of regression estimates from WLS */
		array beta[&p, 1] 	       				/nosym;	/* sampled vector (p x 1) of regression coefficients */
		array CC[&p, &p]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array CI[&p, &p] 					  	/nosym;	/* inverse of CC */
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		call transpose(X, Xt);					/* transpose X */
		call mult(Dbetag, mbetag, prbeta);		/* contribution to posterior mean from prior */
		do k = 1 to &g;							/* cycle through each group independently */
		    do i = 1 to &n;						
			  Yvec[i,1]= Yarr[(k-1)*&n + i];    /* populate nx1 data vector Yvec */
			  Vg[i,i]=nuarr[k]+Sarr[(k-1)*&n+i];/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;					/* off-diagonal elements are those of AR matrix Vgamma */
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);					/* Wg = Vg^{-1} */
			call mult(Xt, Wg, XtW);				/* multiply Xt and Wg */
			call mult(XtW, X, XtWX);			/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(Dbetag,XtWX, DXtWX); /* posterior precision matrix for beta is Dbetag + XtWX */
			call mult(XtW, Yvec, ybeta);		/* contribution to posterior mean from WLS */
			call addmatrix(prbeta,ybeta,pbeta); /* sum of prior and WLS contributions */
			do m = 1 to &p;
				beta[m,1] = rand('normal');		/* sample from univariate standard normal */
			end;
			call chol(DXtWX, CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(CC, CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 1) */
			call transpose(CI, CI);				/* transpose */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 2) */
			call mult(CI, beta, beta);			/* re-scale beta */
			call addmatrix(pbeta, beta, beta);	/* re-center */
			call mult(X, beta, Xbeta);			/* updated vector Xb */
			do i = 1 to &n; 							
		      etamnarr[(k-1)*&n+i] = Xbeta[i,1];/* updated predictions from regression */
			end;
												/* p = 1 here */
			a[k] = beta[1,1];					/* 1st output argument is 1-d array of group-specific intercepts */
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = INDEP_LINEAR or %upcase(&uvar) = FULL_LINEAR %then %do;
	/********************************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the group-specific linear trend model */
	/********************************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_bgl(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  mbetag[*,*], 			/* prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  X[*,*], 				/* design matrix (n x p) using real times */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1arr, etamnarr;				/* arguments that are updated after execution */

		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array Xt[&p, &n]   						/nosym;	/* transpose of design matrix */
		array XtW[&p, &n]					  	/nosym; /* matrix multiplication of Xt and Wg */
		array XtWX[&p, &p] 						/nosym; /* precision matrix of WLS regression estimators */
		array DXtWX[&p, &p]					 	/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array prbeta[&p, 1] 	    			/nosym;	/* vector (p x 1) of regression estimates from prior */
		array pbeta[&p, 1] 	    			   	/nosym;	/* vector (p x 1) of regression estimates from pooled posterior */
		array ybeta[&p, 1] 	       				/nosym;	/* vector (p x 1) of regression estimates from WLS */
		array beta[&p, 1] 	       				/nosym;	/* sampled vector (p x 1) of regression coefficients */
		array CC[&p, &p]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array CI[&p, &p] 					  	/nosym;	/* inverse of CC */
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		call transpose(X, Xt);					/* transpose X */
		call mult(Dbetag, mbetag, prbeta);		/* contribution to posterior mean from prior */
		do k = 1 to &g;							/* cycle through each group independently */
		    do i = 1 to &n;
			  Yvec[i,1]= Yarr[(k-1)*&n + i];    /* populate nx1 data vector Yvec */
			  Vg[i,i]=nuarr[k]+Sarr[(k-1)*&n+i];/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;					/* off-diagonal elements are those of AR matrix Vgamma */
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);					/* Wg = Vg^{-1} */
			call mult(Xt, Wg, XtW);				/* multiply Xt and Wg */
			call mult(XtW, X, XtWX);			/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(Dbetag,XtWX, DXtWX); /* posterior precision matrix for beta is Dbetag + XtWX */
			call mult(XtW, Yvec, ybeta);		/* contribution to posterior mean from WLS */
			call addmatrix(prbeta,ybeta,pbeta); /* sum of prior and WLS contributions */
			do m = 1 to &p;
				beta[m,1] = rand('normal');		/* sample from univariate standard normal */
			end;
			call chol(DXtWX, CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(CC, CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 1) */
			call transpose(CI, CI);				/* transpose */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 2) */
			call mult(CI, beta, beta);			/* re-scale beta */
			call addmatrix(pbeta, beta, beta);	/* re-center */
			call mult(X, beta, Xbeta);			/* updated vector Xb */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i]= Xbeta[i,1]; /* updated predictions from regression */
			end;
												/* p = 2 here */
			a[k] = beta[1,1];					/* 1st output argument is 1-d array of group-specific intercepts */
			b1arr[k] = beta[2,1];				/* 2nd output argument is 1-d array of group-specific linear coefficients */
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = INDEP_QUAD or %upcase(&uvar) = FULL_QUAD %then %do;
	/******************************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the group-specific quad trend model */
	/******************************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_bgq(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  mbetag[*,*], 			/* prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  X[*,*], 				/* design matrix (n x p) using real times */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1arr, b2arr, etamnarr;		/* arguments that are updated after execution */

		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array Xt[&p, &n]   						/nosym;	/* transpose of design matrix */
		array XtW[&p, &n]					  	/nosym; /* matrix multiplication of Xt and Wg */
		array XtWX[&p, &p] 						/nosym; /* precision matrix of WLS regression estimators */
		array DXtWX[&p, &p]					 	/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array prbeta[&p, 1] 	    			/nosym;	/* vector (p x 1) of regression estimates from prior */
		array pbeta[&p, 1] 	    			   	/nosym;	/* vector (p x 1) of regression estimates from pooled posterior */
		array ybeta[&p, 1] 	       				/nosym;	/* vector (p x 1) of regression estimates from WLS */
		array beta[&p, 1] 	       				/nosym;	/* sampled vector (p x 1) of regression coefficients */
		array CC[&p, &p]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array CI[&p, &p] 					  	/nosym;	/* inverse of CC */
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		call transpose(X, Xt);					/* transpose X */
		call mult(Dbetag, mbetag, prbeta);		/* contribution to posterior mean from prior */
		do k = 1 to &g;							/* cycle through each group independently */
		    do i = 1 to &n;	
			  Yvec[i,1]= Yarr[(k-1)*&n + i];    /* populate nx1 data vector Yvec */
			  Vg[i,i]=nuarr[k]+Sarr[(k-1)*&n+i];/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;					/* off-diagonal elements are those of AR matrix Vgamma */
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);					/* Wg = Vg^{-1} */
			call mult(Xt, Wg, XtW);				/* multiply Xt and Wg */
			call mult(XtW, X, XtWX);			/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(Dbetag,XtWX, DXtWX); /* posterior precision matrix for beta is Dbetag + XtWX */
			call mult(XtW, Yvec, ybeta);		/* contribution to posterior mean from WLS */
			call addmatrix(prbeta,ybeta,pbeta); /* sum of prior and WLS contributions */
			do m = 1 to &p;
				beta[m,1] = rand('normal');		/* sample from univariate standard normal */
			end;
			call chol(DXtWX, CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(CC, CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 1) */
			call transpose(CI, CI);				/* transpose */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 2) */
			call mult(CI, beta, beta);			/* re-scale beta */
			call addmatrix(pbeta, beta, beta);	/* re-center */
			call mult(X, beta, Xbeta);			/* updated vector Xb */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i]= Xbeta[i,1]; /* updated predictions from regression */
			end;
												/* p = 3 here */
			a[k] = beta[1,1];					/* 1st output argument is 1-d array of group-specific intercepts */
			b1arr[k] = beta[2,1];				/* 2nd output argument is 1-d array of group-specific linear coefficients */
			b2arr[k] = beta[3,1];				/* 3rd output argument is 1-d array of group-specific quad coefficients */
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = INDEP_CUBIC or %upcase(&uvar) = FULL_CUBIC %then %do;
	/*******************************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the group-specific cubic trend model */
	/*******************************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_bgc(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						  b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						  b3arr[*], 			/* 1-dimensional array (length g) of updated values of cubic coefficients by group */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  mbetag[*,*], 			/* prior mean vector (p x 1) for regression coefficients */
						  Dbetag[*,*], 			/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  X[*,*], 				/* design matrix (n x p) using real times */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1arr,b2arr,b3arr, etamnarr;	/* arguments that are updated after execution */

		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array Xt[&p, &n]   						/nosym;	/* transpose of design matrix */
		array XtW[&p, &n]					  	/nosym; /* matrix multiplication of Xt and Wg */
		array XtWX[&p, &p] 						/nosym; /* precision matrix of WLS regression estimators */
		array DXtWX[&p, &p]					 	/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array prbeta[&p, 1] 	    			/nosym;	/* vector (p x 1) of regression estimates from prior */
		array pbeta[&p, 1] 	    			   	/nosym;	/* vector (p x 1) of regression estimates from pooled posterior */
		array ybeta[&p, 1] 	       				/nosym;	/* vector (p x 1) of regression estimates from WLS */
		array beta[&p, 1] 	       				/nosym;	/* sampled vector (p x 1) of regression coefficients */
		array CC[&p, &p]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array CI[&p, &p] 					  	/nosym;	/* inverse of CC */
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		call transpose(X, Xt);					/* transpose X */
		call mult(Dbetag, mbetag, prbeta);		/* contribution to posterior mean from prior */
		do k = 1 to &g;							/* cycle through each group independently */
		    do i = 1 to &n;
			  Yvec[i,1]= Yarr[(k-1)*&n + i];    /* populate nx1 data vector Yvec */
			  Vg[i,i]=nuarr[k]+Sarr[(k-1)*&n+i];/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;					/* off-diagonal elements are those of AR matrix Vgamma */
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);					/* Wg = Vg^{-1} */
			call mult(Xt, Wg, XtW);				/* multiply Xt and Wg */
			call mult(XtW, X, XtWX);			/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(Dbetag,XtWX, DXtWX); /* posterior precision matrix for beta is Dbetag + XtWX */
			call mult(XtW, Yvec, ybeta);		/* contribution to posterior mean from WLS */
			call addmatrix(prbeta,ybeta,pbeta); /* sum of prior and WLS contributions */
			do m = 1 to &p;
				beta[m,1] = rand('normal');		/* sample from univariate standard normal */
			end;
			call chol(DXtWX, CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(CC, CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 1) */
			call transpose(CI, CI);				/* transpose */
			call mult(CI, pbeta, pbeta);		/* re-scale pbeta (part 2) */
			call mult(CI, beta, beta);			/* re-scale beta */
			call addmatrix(pbeta, beta, beta);	/* re-center */
			call mult(X, beta, Xbeta);			/* updated vector Xb */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i]= Xbeta[i,1]; /* updated predictions from regression */
			end;
												/* p = 4 here */
			a[k] = beta[1,1];					/* 1st output argument is 1-d array of group-specific intercepts */
			b1arr[k] = beta[2,1];				/* 2nd output argument is 1-d array of group-specific linear coefficients */
			b2arr[k] = beta[3,1];				/* 3rd output argument is 1-d array of group-specific quad coefficients */
			b3arr[k] = beta[4,1];				/* 4th output argument is 1-d array of group-specific cubic coefficients */
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = COMMON_LINEAR %then %do;
	/******************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the common linear model */
	/******************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_b1l(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1, 					/* updated value of common linear coefficient */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  ambetag[*,*], 		/* prior mean vector (1 x 1) for intercepts */
						  bmbetag[*,*], 		/* prior mean vector ((p-1) x 1) for remaining coefficients */
						  aDbetag[*,*], 		/* diagonal matrix (1 x 1) of prior precisions for intercepts */
						  bDbetag[*,*], 		/* diagonal matrix ((p-1) x (p-1)) of prior precisions for remaining coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  aX[*,*], 				/* conformal design submatrix (n x 1) */
						  bX[*,*], 				/* conformal design submatrix (n x (p-1)) */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1, etamnarr;				/* arguments that are updated after execution */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array sumbXtWX[%eval(&p-1),%eval(&p-1)] /nosym;	/* cumulative sum of group-specific precision matrices */
		array sumbzbeta[%eval(&p-1), 1] 	    /nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array bXt[%eval(&p-1), &n] 				/nosym; /* transpose of design matrix (excl. intercept) */
		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array bXtW[%eval(&p-1), &n]				/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array bXtWX[%eval(&p-1), %eval(&p-1)] 	/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array bDXtWX[%eval(&p-1), %eval(&p-1)]  /nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array bprbeta[%eval(&p-1), 1]           /nosym; /* vector of regression estimates (excl. intercept) from prior */
		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array bpbeta[%eval(&p-1), 1]            /nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array bzbeta[%eval(&p-1), 1]	 		/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array bbeta[%eval(&p-1), 1]		 	    /nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array bCC[%eval(&p-1), %eval(&p-1)]     /nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array bCI[%eval(&p-1), %eval(&p-1)]     /nosym; /* inverse of CC */

		/********************************/
		/* Update common coefficient(s) */
		/********************************/
		call zeromatrix(sumbXtWX);						/* initialize cumulative sums to all zeroes */
		call zeromatrix(sumbzbeta);
		call transpose(bX, bXt);						/* transpose bX */
		call mult(bDbetag, bmbetag, bprbeta);			/* contribution to posterior mean from prior */
		do k = 1 to &g;									/* cycle through each group independently */
			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(bXt, Wg, bXtW);					/* multiply bXt and Wg */
			call mult(bXtW, bX, bXtWX);					/* calculate bXtWX */
			call addmatrix(sumbXtWX, bXtWX, sumbXtWX);	/* cumulative matrix sum */
			call mult(bXtW, Zvec, bzbeta);			 	/* contributions to posterior mean from WLS */
			call addmatrix(sumbzbeta,bzbeta,sumbzbeta);	/* cumulative matrix sum */
		end;
		call addmatrix(bDbetag, sumbXtWX, bDXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(bprbeta, sumbzbeta, bpbeta);		/* sum of prior and the cumulative WLS contributions */
		do m = 2 to &p;
			bbeta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
		end;
		call chol(bDXtWX, bCC);							/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
		call inv(bCC, bCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 1) */
		call transpose(bCI, bCI);						/* transpose */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 2) */
		call mult(bCI, bbeta, bbeta);					/* re-scale beta */
		call addmatrix(bpbeta, bbeta, bbeta);			/* re-center */
		call mult(bX, bbeta, bXbeta);					/* updated vector bX (used below) */
														/* p=2 here */
		b1 = bbeta[1,1];								/* output argument b1 is updated value of common linear coefficient */

		/************************************************/
		/* Update intercepts and regression predictions */
		/************************************************/
		call transpose(aX, aXt);						/* transpose aX */
		call mult(aDbetag, ambetag, aprbeta);			/* contribution to posterior mean from prior */
		do k = 1 to &g;									/* cycle through each group independently */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - bX */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(aXt, Wg, aXtW);					/* multiply aXt and Wg */
			call mult(aXtW, aX, aXtWX);					/* calculate aXtWX, the precision matrix from WLS */
			call addmatrix(aDbetag,aXtWX,aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
			call mult(aXtW, Zvec, azbeta);				/* contribution to posterior mean from WLS */
			call addmatrix(aprbeta, azbeta, apbeta); 	/* sum of prior and WLS contributions */
			abeta[1,1] = rand('normal');				/* sample intercept from univariate normal */
			call chol(aDXtWX, aCC);						/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(aCC, aCI);							/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 1) */
			call transpose(aCI, aCI);					/* transpose */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 2) */
			call mult(aCI, abeta, abeta);				/* re-scale beta */
			call addmatrix(apbeta, abeta, abeta);		/* re-center */
			a[k] = abeta[1,1];							/* output argument a is 1-dimensional array of group-specific intercepts */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
			end;
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = COMMON_QUAD %then %do;
	/**********************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the common quad trend model */
	/**********************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_b1q(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1, 					/* updated value of common linear coefficient */
						  b2, 					/* updated value of common quad coefficient */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  ambetag[*,*], 		/* prior mean vector (1 x 1) for intercepts */
						  bmbetag[*,*], 		/* prior mean vector ((p-1) x 1) for remaining coefficients */
						  aDbetag[*,*], 		/* diagonal matrix (1 x 1) of prior precisions for intercepts */
						  bDbetag[*,*], 		/* diagonal matrix ((p-1) x (p-1)) of prior precisions for remaining coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  aX[*,*], 				/* conformal design submatrix (n x 1) */
						  bX[*,*], 				/* conformal design submatrix (n x (p-1)) */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1, b2, etamnarr;			/* arguments that are updated after execution */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array sumbXtWX[%eval(&p-1),%eval(&p-1)] /nosym;	/* cumulative sum of group-specific precision matrices */
		array sumbzbeta[%eval(&p-1), 1] 	    /nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array bXt[%eval(&p-1), &n] 				/nosym; /* transpose of design matrix (excl. intercept) */
		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array bXtW[%eval(&p-1), &n]				/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array bXtWX[%eval(&p-1), %eval(&p-1)] 	/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array bDXtWX[%eval(&p-1), %eval(&p-1)]  /nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array bprbeta[%eval(&p-1), 1]           /nosym; /* vector of regression estimates (excl. intercept) from prior */
		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array bpbeta[%eval(&p-1), 1]            /nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array bzbeta[%eval(&p-1), 1]	 		/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array bbeta[%eval(&p-1), 1]		 	    /nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array bCC[%eval(&p-1), %eval(&p-1)]     /nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array bCI[%eval(&p-1), %eval(&p-1)]     /nosym; /* inverse of CC */

		/********************************/
		/* Update common coefficient(s) */
		/********************************/
		call zeromatrix(sumbXtWX);						/* initialize cumulative sums to all zeroes */
		call zeromatrix(sumbzbeta);
		call transpose(bX, bXt);						/* transpose bX */
		call mult(bDbetag, bmbetag, bprbeta);			/* contribution to posterior mean from prior */	
		do k = 1 to &g;									/* cycle through each group independently */
			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(bXt, Wg, bXtW);					/* multiply bXt and Wg */
			call mult(bXtW, bX, bXtWX);					/* calculate bXtWX */
			call addmatrix(sumbXtWX, bXtWX, sumbXtWX);	/* cumulative matrix sum */
			call mult(bXtW, Zvec, bzbeta);			 	/* contributions to posterior mean from WLS */
			call addmatrix(sumbzbeta,bzbeta,sumbzbeta);	/* cumulative matrix sum */
		end;
		call addmatrix(bDbetag, sumbXtWX, bDXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(bprbeta, sumbzbeta, bpbeta);		/* sum of prior and the cumulative WLS contributions */
		do m = 2 to &p;
			bbeta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
		end;
		call chol(bDXtWX, bCC);							/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
		call inv(bCC, bCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 1) */
		call transpose(bCI, bCI);						/* transpose */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 2) */
		call mult(bCI, bbeta, bbeta);					/* re-scale beta */
		call addmatrix(bpbeta, bbeta, bbeta);			/* re-center */
		call mult(bX, bbeta, bXbeta);					/* updated vector bX (used below) */
														/* p=3 here */
		b1 = bbeta[1,1];								/* output argument b1 is updated value of common linear coefficient */
		b2 = bbeta[2,1];								/* output argument b2 is updated value of common quad coefficient */

		/************************************************/
		/* Update intercepts and regression predictions */
		/************************************************/
		call transpose(aX, aXt);						/* transpose aX */
		call mult(aDbetag, ambetag, aprbeta);			/* contribution to posterior mean from prior */
		do k = 1 to &g;									/* cycle through each group independently */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - bX */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(aXt, Wg, aXtW);					/* multiply aXt and Wg */
			call mult(aXtW, aX, aXtWX);					/* calculate aXtWX, the precision matrix from WLS */
			call addmatrix(aDbetag,aXtWX,aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
			call mult(aXtW, Zvec, azbeta);				/* contribution to posterior mean from WLS */
			call addmatrix(aprbeta, azbeta, apbeta); 	/* sum of prior and WLS contributions */
			abeta[1,1] = rand('normal');				/* sample intercept from univariate normal */
			call chol(aDXtWX, aCC);						/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(aCC, aCI);							/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 1) */
			call transpose(aCI, aCI);					/* transpose */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 2) */
			call mult(aCI, abeta, abeta);				/* re-scale beta */
			call addmatrix(apbeta, abeta, abeta);		/* re-center */
			a[k] = abeta[1,1];							/* output argument a is 1-dimensional array of group-specific intercepts */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
			end;
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = COMMON_CUBIC %then %do;
	/***********************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the common cubic trend model */
	/***********************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_b1c(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						  b1, 					/* updated value of common linear coefficient */
						  b2, 					/* updated value of common quad coefficient */
						  b3, 					/* updated value of common cubic coefficient */
						  etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						  ambetag[*,*], 		/* prior mean vector (1 x 1) for intercepts */
						  bmbetag[*,*], 		/* prior mean vector ((p-1) x 1) for remaining coefficients */
						  aDbetag[*,*], 		/* diagonal matrix (1 x 1) of prior precisions for intercepts */
						  bDbetag[*,*], 		/* diagonal matrix ((p-1) x (p-1)) of prior precisions for remaining coefficients */
						  rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						  nuarr[*],				/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						  rts[*],				/* 1-dimensional array (length n) of real times */
						  aX[*,*], 				/* conformal design submatrix (n x 1) */
						  bX[*,*], 				/* conformal design submatrix (n x (p-1)) */
						  Yarr[*], 				/* 1-dimensional array (length gn) for _y from dataset */
						  Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						  );

		outargs a, b1, b2, b3, etamnarr;		/* arguments that are updated after execution */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */
		array sumbXtWX[%eval(&p-1),%eval(&p-1)] /nosym;	/* cumulative sum of group-specific precision matrices */
		array sumbzbeta[%eval(&p-1), 1] 	    /nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array bXt[%eval(&p-1), &n] 				/nosym; /* transpose of design matrix (excl. intercept) */
		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array bXtW[%eval(&p-1), &n]				/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array bXtWX[%eval(&p-1), %eval(&p-1)] 	/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array bDXtWX[%eval(&p-1), %eval(&p-1)]  /nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array bprbeta[%eval(&p-1), 1]           /nosym; /* vector of regression estimates (excl. intercept) from prior */
		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array bpbeta[%eval(&p-1), 1]            /nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array bzbeta[%eval(&p-1), 1]	 		/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array bbeta[%eval(&p-1), 1]		 	    /nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array bCC[%eval(&p-1), %eval(&p-1)]     /nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array bCI[%eval(&p-1), %eval(&p-1)]     /nosym; /* inverse of CC */

		/********************************/
		/* Update common coefficient(s) */
		/********************************/
		call zeromatrix(sumbXtWX);						/* initialize cumulative sums to all zeroes */
		call zeromatrix(sumbzbeta);
		call transpose(bX, bXt);						/* transpose bX */
		call mult(bDbetag, bmbetag, bprbeta);			/* contribution to posterior mean from prior */		
		do k = 1 to &g;									/* cycle through each group independently */
			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(bXt, Wg, bXtW);					/* multiply bXt and Wg */
			call mult(bXtW, bX, bXtWX);					/* calculate bXtWX */
			call addmatrix(sumbXtWX, bXtWX, sumbXtWX);	/* cumulative matrix sum */
			call mult(bXtW, Zvec, bzbeta);			 	/* contributions to posterior mean from WLS */
			call addmatrix(sumbzbeta,bzbeta,sumbzbeta);	/* cumulative matrix sum */
		end;
		call addmatrix(bDbetag, sumbXtWX, bDXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(bprbeta, sumbzbeta, bpbeta);		/* sum of prior and the cumulative WLS contributions */
		do m = 2 to &p;
			bbeta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
		end;
		call chol(bDXtWX, bCC);							/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
		call inv(bCC, bCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 1) */
		call transpose(bCI, bCI);						/* transpose */
		call mult(bCI, bpbeta, bpbeta);					/* re-scale pbeta (part 2) */
		call mult(bCI, bbeta, bbeta);					/* re-scale beta */
		call addmatrix(bpbeta, bbeta, bbeta);			/* re-center */
		call mult(bX, bbeta, bXbeta);					/* updated vector bX (used below) */
														/* p=4 here */
		b1 = bbeta[1,1];								/* output argument b1 is updated value of common linear coefficient */
		b2 = bbeta[2,1];								/* output argument b2 is updated value of common quad coefficient */
		b3 = bbeta[3,1];								/* output argument b3 is updated value of common cubic coefficient */

		/************************************************/
		/* Update intercepts and regression predictions */
		/************************************************/
		call transpose(aX, aXt);						/* transpose aX */
		call mult(aDbetag, ambetag, aprbeta);			/* contribution to posterior mean from prior */
		do k = 1 to &g;									/* cycle through each group independently */
		    do i = 1 to &n;						
			  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];/* populate nx1 data vector Zvec = Yvec - bX */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */
			call mult(aXt, Wg, aXtW);					/* multiply aXt and Wg */
			call mult(aXtW, aX, aXtWX);					/* calculate aXtWX, the precision matrix from WLS */
			call addmatrix(aDbetag,aXtWX,aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
			call mult(aXtW, Zvec, azbeta);				/* contribution to posterior mean from WLS */
			call addmatrix(aprbeta, azbeta, apbeta); 	/* sum of prior and WLS contributions */
			abeta[1,1] = rand('normal');				/* sample intercept from univariate normal */
			call chol(aDXtWX, aCC);						/* Cholesky decomposition for precision matrix (returns lower triangular) */
			call inv(aCC, aCI);							/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 1) */
			call transpose(aCI, aCI);					/* transpose */
			call mult(aCI, apbeta, apbeta);				/* re-scale pbeta (part 2) */
			call mult(aCI, abeta, abeta);				/* re-scale beta */
			call addmatrix(apbeta, abeta, abeta);		/* re-center */
			a[k] = abeta[1,1];							/* output argument a is 1-dimensional array of group-specific intercepts */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
			do i = 1 to &n; 								
		      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
			end;
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = BMA_LINEAR %then %do;
	/*********************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the BMA linear trend model */
	/*********************************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine CP_bmal(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						   b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						   b1, 					/* updated value of common linear coefficient */
						   etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*],				/* 1-dimensional array (length gn) for _var from dataset */
						   flg					/* model flag (1...3 in the linear case) */
						   );

		outargs a, b1arr, b1, etamnarr;			/* arguments that are updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		/*************************************************************************************/
		/* Array structures for indep trend models in the full dimensional linear BMA: p = 2 */
		/*************************************************************************************/
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array q2X[&n, 2]						/nosym; /* 2-column version of the design matrix X */
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */

		array q2mbetag[2, 1]					/nosym; /* 2-dimensional version of mbetag */
		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */

		array q2Dbetag[2, 2]					/nosym; /* 2-dimensional version of Dbetag */
		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */

		array q2Xt[2, &n]   					/nosym;	/* transpose of design matrix */
		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */

		array q2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators */
		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */

		array q2DXtWX[2, 2]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q2prbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from prior */
		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */

		array q2pbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from pooled posterior */
		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */

		array q2ybeta[2, 1] 	       			/nosym;	/* vector (2 x 1) of regression estimates from WLS */
		array q1ybeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */

		array q2beta[2, 1] 	       				/nosym;	/* sampled vector (2 x 1) of regression coefficients */
		array q1beta[1, 1] 	       				/nosym;	/* sampled vector (1 x 1) of regression coefficients */

		array q2CC[2, 2]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q1CC[1, 1]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */

		array q2CI[2, 2] 					  	/nosym;	/* inverse of CC */
		array q1CI[1, 1] 					  	/nosym;	/* inverse of CC */

		/***************************************************************************************/
		/* Array structures for common trend models in the full dimensional linear BMA: p = 2  */
		/***************************************************************************************/
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array ambetag[1, 1]						/nosym; /* prior mean vector (1 x 1) for intercepts */
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for remaining coefficients */

		array aDbetag[1, 1] 					/nosym; /* diagonal matrix (1 x 1) of prior precisions for intercepts */
		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for remaining coefficients */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */

		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array b1beta[1, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */

		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array b1CC[1, 1]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */

		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array b1CI[1, 1]     					/nosym; /* inverse of CC */

		/****************************************************************/
		/* Populate needed array structures depending on the model flag */	
		/****************************************************************/
		if flg = 1 then do;								/* indep linear: q = 2 */
			do i = 1 to &n;
				do m = 1 to 2;
					q2X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q2Dbetag);
			do m = 1 to 2;
				q2mbetag[m, 1] = mbetag[m, 1];
				q2Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q2X, q2Xt);					/* transpose qX */
			call mult(q2Dbetag, q2mbetag, q2prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 3 then do;								/* dropped: q = 1 */
			do i = 1 to &n;
				do m = 1 to 1;
					q1X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q1Dbetag);
			do m = 1 to 1;
				q1mbetag[m, 1] = mbetag[m, 1];
				q1Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q1X, q1Xt);					/* transpose qX */
			call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 2 then do;								/* common linear: q = 2 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 2;
					b1X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b1Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 2;
				b1mbetag[m-1, 1]   = mbetag[m, 1];
			    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b1X, b1Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb1zbeta);	
		end;

		/*******************************/
		/* Group-specific trend models */
		/*******************************/
		if flg = 1 or flg = 3 then do;
			/**********************************/
			/* Update regression coefficients */
			/**********************************/
			do k = 1 to &g;									/* cycle through each group independently */
			    do i = 1 to &n;								
					Yvec[i,1] = Yarr[(k-1)*&n + i]; 		/* populate nx1 data vector Yvec */
				    Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
				    do j = i+1 to &n;
					    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
					    Vg[j,i] = Vg[i,j];
				    end; 
			  	end; 
				call inv(Vg, Wg);							/* Wg = Vg^{-1} */
				if flg = 1 then do;							/* group-specific linear trend model (q=2)*/
					call mult(q2Xt, Wg, q2XtW);				/* multiply Xt and Wg */
					call mult(q2XtW, q2X, q2XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q2Dbetag,q2XtWX,q2DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q2XtW, Yvec, q2ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q2prbeta,q2ybeta,q2pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 2;
						q2beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q2DXtWX, q2CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q2CC, q2CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q2CI, q2CI);				/* transpose */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 2) */
					call mult(q2CI, q2beta, q2beta);		/* re-scale beta */
					call addmatrix(q2pbeta,q2beta,q2beta);	/* re-center */
					call mult(q2X, q2beta, Xbeta);			/* updated vector Xb */
					a[k] = q2beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q2beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
				end;
				if flg = 3 then do;							/* group-specific intercept-only model (q=1)*/
					call mult(q1Xt, Wg, q1XtW);				/* multiply Xt and Wg */
					call mult(q1XtW, q1X, q1XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q1Dbetag,q1XtWX,q1DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q1XtW, Yvec, q1ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q1prbeta,q1ybeta,q1pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 1;
						q1beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q1DXtWX, q1CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q1CC, q1CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q1CI, q1CI);				/* transpose */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 2) */
					call mult(q1CI, q1beta, q1beta);		/* re-scale beta */
					call addmatrix(q1pbeta,q1beta, q1beta);	/* re-center */
					call mult(q1X, q1beta, Xbeta);			/* updated vector Xb */
					a[k] = q1beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = 0;							/* 2nd output argument is 1-d array of group-specific linear coefficients */
				end;
				do i = 1 to &n; 								
			     	etamnarr[(k-1)*&n+i]= Xbeta[i,1]; 		/* updated predictions from regression */
				end;
			end;
			tmpb1 = 0;
			do k = 1 to &g;
				tmpb1 = tmpb1 + b1arr[k]; 		 			/* common coefficient(s) set to average of group-specific coefficients */
			end;
			b1 = tmpb1/&g;
		end;

		/**********************/
		/* Common trend model */
		/**********************/
		if flg = 2 then do;
			/********************************/
			/* Update common coefficient(s) */
			/********************************/
			do k = 1 to &g;										/* cycle through each group independently */
				abeta[1,1] = a[k];								/* group-specific abeta vector */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - a */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
			    call mult(b1Xt, Wg, b1XtW);						/* multiply bXt and Wg */
			    call mult(b1XtW, b1X, b1XtWX);					/* calculate bXtWX */
			    call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
			    call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
			    call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);	/* cumulative matrix sum */
			end;												/* end cycle through groups */
			call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
			call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);		/* sum of prior and the cumulative WLS contributions */
			do m = 2 to 2;
				b1beta[m-1,1] = rand('normal');					/* sample from univariate standard normal(s) */
			end;
			call chol(b1DXtWX, b1CC);							/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
			call inv(b1CC, b1CI);								/* inverse of lower triangular matrix from Cholesky decomposition */
			call mult(b1CI, b1pbeta, b1pbeta);					/* re-scale pbeta (part 1) */
			call transpose(b1CI, b1CI);							/* transpose */
			call mult(b1CI, b1pbeta, b1pbeta);					/* re-scale pbeta (part 2) */
			call mult(b1CI, b1beta, b1beta);					/* re-scale beta */
			call addmatrix(b1pbeta, b1beta, b1beta);			/* re-center */
			call mult(b1X, b1beta, bXbeta);						/* updated vector bX */
			b1 = b1beta[1,1];									/* output argument b1 is updated value of common linear coefficient */
			do k = 1 to &g;										/* arrays of group-specific coefficient also updated to reflect common value */
				b1arr[k] = b1;
			end;
			/************************************************/
			/* Update intercepts and regression predictions */
			/************************************************/
			do k = 1 to &g;
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - bX */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
				call mult(aXt, Wg, aXtW);						/* multiply aXt and Wg */
				call mult(aXtW, aX, aXtWX);						/* calculate aXtWX, the precision matrix from WLS */
				call addmatrix(aDbetag, aXtWX, aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
				call mult(aXtW, Zvec, azbeta);					/* contribution to posterior mean from WLS */
				call addmatrix(aprbeta, azbeta, apbeta); 		/* sum of prior and WLS contributions */
				abeta[1,1] = rand('normal');					/* sample intercept from univariate normal */
				call chol(aDXtWX, aCC);							/* Cholesky decomposition for precision matrix (returns lower triangular) */
				call inv(aCC, aCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 1) */
				call transpose(aCI, aCI);						/* transpose */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 2) */
				call mult(aCI, abeta, abeta);					/* re-scale beta */
				call addmatrix(apbeta, abeta, abeta);			/* re-center */
				a[k] = abeta[1,1];								/* output argument a is 1-dimensional array of group-specific intercepts */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
				do i = 1 to &n; 								
			      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
				end;
			end;
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = BMA_QUAD %then %do;
	/*******************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the BMA quad trend model */
	/*******************************************************************************/
	proc fcmp outlib=&uloc;		

		subroutine CP_bmaq(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						   b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						   b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						   b1, 					/* updated value of common linear coefficient */
						   b2, 					/* updated value of common quad coefficient */
						   etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*],				/* 1-dimensional array (length gn) for _var from dataset */
						   flg					/* model flag (1...5 in the quad case) */
						   );

		outargs a, b1arr, b2arr, b1, b2, etamnarr;		/* arguments that are updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		/***********************************************************************************/
		/* Array structures for indep trend models in the full dimensional quad BMA: p = 3 */
		/***********************************************************************************/
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array q3X[&n, 3]						/nosym; /* 3-column version of the design matrix X */
		array q2X[&n, 2]						/nosym; /* 2-column version of the design matrix X */
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */

		array q3mbetag[3, 1]					/nosym; /* 3-dimensional version of mbetag */
		array q2mbetag[2, 1]					/nosym; /* 2-dimensional version of mbetag */
		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */

		array q3Dbetag[3, 3]					/nosym; /* 3-dimensional version of Dbetag */
		array q2Dbetag[2, 2]					/nosym; /* 2-dimensional version of Dbetag */
		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */

		array q3Xt[3, &n]   					/nosym;	/* transpose of design matrix */
		array q2Xt[2, &n]   					/nosym;	/* transpose of design matrix */
		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */

		array q3XtW[3, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q3XtWX[3, 3] 						/nosym; /* precision matrix of WLS regression estimators */
		array q2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators */
		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */

		array q3DXtWX[3, 3]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q2DXtWX[2, 2]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q3prbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from prior */
		array q2prbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from prior */
		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */

		array q3pbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from pooled posterior */
		array q2pbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from pooled posterior */
		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */

		array q3ybeta[3, 1] 	       			/nosym;	/* vector (3 x 1) of regression estimates from WLS */
		array q2ybeta[2, 1] 	       			/nosym;	/* vector (2 x 1) of regression estimates from WLS */
		array q1ybeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */

		array q3beta[3, 1] 	       				/nosym;	/* sampled vector (3 x 1) of regression coefficients */
		array q2beta[2, 1] 	       				/nosym;	/* sampled vector (2 x 1) of regression coefficients */
		array q1beta[1, 1] 	       				/nosym;	/* sampled vector (1 x 1) of regression coefficients */

		array q3CC[3, 3]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q2CC[2, 2]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q1CC[1, 1]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */

		array q3CI[3, 3] 					  	/nosym;	/* inverse of CC */
		array q2CI[2, 2] 					  	/nosym;	/* inverse of CC */
		array q1CI[1, 1] 					  	/nosym;	/* inverse of CC */

		/*************************************************************************************/
		/* Array structures for common trend models in the full dimensional quad BMA: p = 3  */
		/*************************************************************************************/
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array ambetag[1, 1]						/nosym; /* prior mean vector (1 x 1) for intercepts */
		array b2mbetag[2, 1] 					/nosym;	/* prior mean vector (2 x 1) for remaining coefficients */
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for remaining coefficients */

		array aDbetag[1, 1] 					/nosym; /* diagonal matrix (1 x 1) of prior precisions for intercepts */
		array b2Dbetag[2, 2] 					/nosym;	/* diagonal matrix (2 x 2) of prior precisions for remaining coefficients */
		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for remaining coefficients */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array b2X[&n, 2]						/nosym; /* 2-dimensional conformal design submatrix X */
		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */

		array sumb2XtWX[2, 2] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb2zbeta[2, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array b2Xt[2, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array b2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array b2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array b2DXtWX[2, 2]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array b2prbeta[2, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array b2pbeta[2, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array b2zbeta[2, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array b2beta[2, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array b1beta[1, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */

		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array b2CC[2, 2]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array b1CC[1, 1]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */

		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array b2CI[2, 2]     					/nosym; /* inverse of CC */
		array b1CI[1, 1]     					/nosym; /* inverse of CC */

		/****************************************************************/
		/* Populate needed array structures depending on the model flag */	
		/****************************************************************/
		if flg = 1 then do;								/* indep quad: q = 3 */
			do i = 1 to &n;
				do m = 1 to 3;
					q3X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q3Dbetag);
			do m = 1 to 3;
				q3mbetag[m, 1] = mbetag[m, 1];
				q3Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q3X, q3Xt);					/* transpose qX */
			call mult(q3Dbetag, q3mbetag, q3prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 2 then do;								/* indep linear: q = 2 */
			do i = 1 to &n;
				do m = 1 to 2;
					q2X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q2Dbetag);
			do m = 1 to 2;
				q2mbetag[m, 1] = mbetag[m, 1];
				q2Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q2X, q2Xt);					/* transpose qX */
			call mult(q2Dbetag, q2mbetag, q2prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 5 then do;								/* dropped: q = 1 */
			do i = 1 to &n;
				do m = 1 to 1;
					q1X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q1Dbetag);
			do m = 1 to 1;
				q1mbetag[m, 1] = mbetag[m, 1];
				q1Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q1X, q1Xt);					/* transpose qX */
			call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 3 then do;								/* common quad: q = 3 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 3;
					b2X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b2Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 3;
				b2mbetag[m-1, 1]   = mbetag[m, 1];
			    b2Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b2X, b2Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b2Dbetag, b2mbetag, b2prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb2XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb2zbeta);	
		end;
		if flg = 4 then do;								/* common linear: q = 2 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 2;
					b1X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b1Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 2;
				b1mbetag[m-1, 1]   = mbetag[m, 1];
			    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b1X, b1Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb1zbeta);	
		end;

		/*******************************/
		/* Group-specific trend models */
		/*******************************/
		if flg = 1 or flg = 2 or flg = 5 then do;
			/**********************************/
			/* Update regression coefficients */
			/**********************************/
			do k = 1 to &g;									/* cycle through each group independently */
			    do i = 1 to &n;								
					Yvec[i,1] = Yarr[(k-1)*&n + i]; 		/* populate nx1 data vector Yvec */
				    Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
				    do j = i+1 to &n;
					    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
					    Vg[j,i] = Vg[i,j];
				    end; 
			  	end; 
				call inv(Vg, Wg);							/* Wg = Vg^{-1} */
				if flg = 1 then do;							/* group-specific quad trend model (q=3)*/
					call mult(q3Xt, Wg, q3XtW);				/* multiply Xt and Wg */
					call mult(q3XtW, q3X, q3XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q3Dbetag,q3XtWX,q3DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q3XtW, Yvec, q3ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q3prbeta,q3ybeta,q3pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 3;
						q3beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q3DXtWX, q3CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q3CC, q3CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q3CI, q3pbeta, q3pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q3CI, q3CI);				/* transpose */
					call mult(q3CI, q3pbeta, q3pbeta);		/* re-scale pbeta (part 2) */
					call mult(q3CI, q3beta, q3beta);		/* re-scale beta */
					call addmatrix(q3pbeta,q3beta,q3beta);	/* re-center */
					call mult(q3X, q3beta, Xbeta);			/* updated vector Xb */
					a[k] = q3beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q3beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = q3beta[3,1];					/* 3rd output argument is 1-d array of group-specific quad coefficients */
				end;
				if flg = 2 then do;							/* group-specific linear trend model (q=2)*/
					call mult(q2Xt, Wg, q2XtW);				/* multiply Xt and Wg */
					call mult(q2XtW, q2X, q2XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q2Dbetag,q2XtWX,q2DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q2XtW, Yvec, q2ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q2prbeta,q2ybeta,q2pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 2;
						q2beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q2DXtWX, q2CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q2CC, q2CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q2CI, q2CI);				/* transpose */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 2) */
					call mult(q2CI, q2beta, q2beta);		/* re-scale beta */
					call addmatrix(q2pbeta,q2beta,q2beta);	/* re-center */
					call mult(q2X, q2beta, Xbeta);			/* updated vector Xb */
					a[k] = q2beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q2beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = 0;							/* 3rd output argument is 1-d array of group-specific quad coefficients */
				end;
				if flg = 5 then do;							/* group-specific intercept-only model (q=1)*/
					call mult(q1Xt, Wg, q1XtW);				/* multiply Xt and Wg */
					call mult(q1XtW, q1X, q1XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q1Dbetag,q1XtWX,q1DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q1XtW, Yvec, q1ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q1prbeta,q1ybeta,q1pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 1;
						q1beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q1DXtWX, q1CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q1CC, q1CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q1CI, q1CI);				/* transpose */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 2) */
					call mult(q1CI, q1beta, q1beta);		/* re-scale beta */
					call addmatrix(q1pbeta,q1beta,q1beta);	/* re-center */
					call mult(q1X, q1beta, Xbeta);			/* updated vector Xb */
					a[k] = q1beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = 0;							/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = 0;							/* 3rd output argument is 1-d array of group-specific quad coefficients */
				end;
				do i = 1 to &n; 								
			     	etamnarr[(k-1)*&n+i]= Xbeta[i,1]; 		/* updated predictions from regression */
				end;
			end;
			tmpb1 = 0;
			tmpb2 = 0;
			do k = 1 to &g;
				tmpb1 = tmpb1 + b1arr[k]; 		 			/* common coefficient(s) set to average of group-specific coefficients */
				tmpb2 = tmpb2 + b2arr[k];
			end;
			b1 = tmpb1/&g;
			b2 = tmpb2/&g;
		end;

		/***********************/
		/* Common trend models */
		/***********************/
		if flg = 3 or flg = 4 then do;
			/********************************/
			/* Update common coefficient(s) */
			/********************************/
			do k = 1 to &g;										/* cycle through each group independently */
				abeta[1,1] = a[k];								/* group-specific abeta vector */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - a */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
				if flg = 3 then do;								/* common quad trend model */
				  call mult(b2Xt, Wg, b2XtW);					/* multiply bXt and Wg */
				  call mult(b2XtW, b2X, b2XtWX);				/* calculate bXtWX */
				  call addmatrix(sumb2XtWX, b2XtWX, sumb2XtWX);	/* cumulative matrix sum */
				  call mult(b2XtW, Zvec, b2zbeta);			 	/* contributions to posterior mean from WLS */
				  call addmatrix(sumb2zbeta,b2zbeta,sumb2zbeta);/* cumulative matrix sum */
				end;
				if flg = 4 then do;								/* common linear trend model */
				  call mult(b1Xt, Wg, b1XtW);					/* multiply bXt and Wg */
				  call mult(b1XtW, b1X, b1XtWX);				/* calculate bXtWX */
				  call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
				  call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
				  call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);/* cumulative matrix sum */
				end;
			end;												/* end cycle through groups */
			if flg = 3 then do;									/* common quad trend model */
				call addmatrix(b2Dbetag, sumb2XtWX, b2DXtWX); 	/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
				call addmatrix(b2prbeta, sumb2zbeta, b2pbeta);	/* sum of prior and the cumulative WLS contributions */
				do m = 2 to 3;
					b2beta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
				end;
				call chol(b2DXtWX, b2CC);						/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
				call inv(b2CC, b2CI);							/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(b2CI, b2pbeta, b2pbeta);				/* re-scale pbeta (part 1) */
				call transpose(b2CI, b2CI);						/* transpose */
				call mult(b2CI, b2pbeta, b2pbeta);				/* re-scale pbeta (part 2) */
				call mult(b2CI, b2beta, b2beta);				/* re-scale beta */
				call addmatrix(b2pbeta, b2beta, b2beta);		/* re-center */
				call mult(b2X, b2beta, bXbeta);					/* updated vector bX */
				b1 = b2beta[1,1];								/* output argument b1 is updated value of common linear coefficient */
				b2 = b2beta[2,1];								/* output argument b2 is updated value of common quad coefficient */
			end;
			if flg = 4 then do;									/* common linear trend model */
				call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 	/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
				call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);	/* sum of prior and the cumulative WLS contributions */
				do m = 2 to 2;
					b1beta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
				end;
				call chol(b1DXtWX, b1CC);						/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
				call inv(b1CC, b1CI);							/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(b1CI, b1pbeta, b1pbeta);				/* re-scale pbeta (part 1) */
				call transpose(b1CI, b1CI);						/* transpose */
				call mult(b1CI, b1pbeta, b1pbeta);				/* re-scale pbeta (part 2) */
				call mult(b1CI, b1beta, b1beta);				/* re-scale beta */
				call addmatrix(b1pbeta, b1beta, b1beta);		/* re-center */
				call mult(b1X, b1beta, bXbeta);					/* updated vector bX */
				b1 = b1beta[1,1];								/* output argument b1 is updated value of common linear coefficient */
				b2 = 0;											/* output argument b2 is updated value of common quad coefficient */
			end;
			do k = 1 to &g;										/* arrays of group-specific coefficient also updated to reflect common value */
				b1arr[k] = b1;
				b2arr[k] = b2;
			end;
			/************************************************/
			/* Update intercepts and regression predictions */
			/************************************************/
			do k = 1 to &g;
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - bX */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
				call mult(aXt, Wg, aXtW);						/* multiply aXt and Wg */
				call mult(aXtW, aX, aXtWX);						/* calculate aXtWX, the precision matrix from WLS */
				call addmatrix(aDbetag, aXtWX, aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
				call mult(aXtW, Zvec, azbeta);					/* contribution to posterior mean from WLS */
				call addmatrix(aprbeta, azbeta, apbeta); 		/* sum of prior and WLS contributions */
				abeta[1,1] = rand('normal');					/* sample intercept from univariate normal */
				call chol(aDXtWX, aCC);							/* Cholesky decomposition for precision matrix (returns lower triangular) */
				call inv(aCC, aCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 1) */
				call transpose(aCI, aCI);						/* transpose */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 2) */
				call mult(aCI, abeta, abeta);					/* re-scale beta */
				call addmatrix(apbeta, abeta, abeta);			/* re-center */
				a[k] = abeta[1,1];								/* output argument a is 1-dimensional array of group-specific intercepts */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
				do i = 1 to &n; 								
			      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
				end;
			end;
		end;

		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = BMA_CUBIC %then %do;
	/********************************************************************************/
	/* eMKF: Gibbs sampler for regression coefficients in the BMA cubic trend model */
	/********************************************************************************/
	proc fcmp outlib=&uloc;		

		subroutine CP_bmac(a[*], 				/* 1-dimensional array (length g) of updated values of intercepts by group */
						   b1arr[*], 			/* 1-dimensional array (length g) of updated values of linear coefficients by group */
						   b2arr[*], 			/* 1-dimensional array (length g) of updated values of quad coefficients by group */
						   b3arr[*], 			/* 1-dimensional array (length g) of updated values of cubic coefficients by group */
						   b1, 					/* updated value of common linear coefficient */
						   b2, 					/* updated value of common quad coefficient */
						   b3, 					/* updated value of common cubic coefficient */
						   etamnarr[*],			/* 1-dimensional array (length gn) of updated values of regression predictions */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*],				/* 1-dimensional array (length gn) for _var from dataset */
						   flg					/* model flag (1...7 in the cubic case) */
						   );

		outargs a,b1arr,b2arr,b3arr,b1,b2,b3,etamnarr;	/* arguments that are updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array Yvec[&n, 1]						/nosym; /* vector (nx1) for use in calculations */
		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		/***********************************************************************************/
		/* Array structures for indep trend models in the full dimensional cubic BMA: p = 4*/
		/***********************************************************************************/
		array Xbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array q4X[&n, 4]						/nosym; /* 4-column version of the design matrix X */
		array q3X[&n, 3]						/nosym; /* 3-column version of the design matrix X */
		array q2X[&n, 2]						/nosym; /* 2-column version of the design matrix X */
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */

		array q4mbetag[4, 1]					/nosym; /* 4-dimensional version of mbetag */
		array q3mbetag[3, 1]					/nosym; /* 3-dimensional version of mbetag */
		array q2mbetag[2, 1]					/nosym; /* 2-dimensional version of mbetag */
		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */

		array q4Dbetag[4, 4]					/nosym; /* 4-dimensional version of Dbetag */
		array q3Dbetag[3, 3]					/nosym; /* 3-dimensional version of Dbetag */
		array q2Dbetag[2, 2]					/nosym; /* 2-dimensional version of Dbetag */
		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */

		array q4Xt[4, &n]   					/nosym;	/* transpose of design matrix */
		array q3Xt[3, &n]   					/nosym;	/* transpose of design matrix */
		array q2Xt[2, &n]   					/nosym;	/* transpose of design matrix */
		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */

		array q4XtW[4, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q3XtW[3, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q4XtWX[4, 4] 						/nosym; /* precision matrix of WLS regression estimators */
		array q3XtWX[3, 3] 						/nosym; /* precision matrix of WLS regression estimators */
		array q2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators */
		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */

		array q4DXtWX[4, 4]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q3DXtWX[3, 3]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q2DXtWX[2, 2]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q4prbeta[4, 1] 	    			/nosym;	/* vector (4 x 1) of regression estimates from prior */
		array q3prbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from prior */
		array q2prbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from prior */
		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */

		array q4pbeta[4, 1] 	    			/nosym;	/* vector (4 x 1) of regression estimates from pooled posterior */
		array q3pbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from pooled posterior */
		array q2pbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from pooled posterior */
		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */

		array q4ybeta[4, 1] 	       			/nosym;	/* vector (4 x 1) of regression estimates from WLS */
		array q3ybeta[3, 1] 	       			/nosym;	/* vector (3 x 1) of regression estimates from WLS */
		array q2ybeta[2, 1] 	       			/nosym;	/* vector (2 x 1) of regression estimates from WLS */
		array q1ybeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */

		array q4beta[4, 1] 	       				/nosym;	/* sampled vector (4 x 1) of regression coefficients */
		array q3beta[3, 1] 	       				/nosym;	/* sampled vector (3 x 1) of regression coefficients */
		array q2beta[2, 1] 	       				/nosym;	/* sampled vector (2 x 1) of regression coefficients */
		array q1beta[1, 1] 	       				/nosym;	/* sampled vector (1 x 1) of regression coefficients */

		array q4CC[4, 4]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q3CC[3, 3]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q2CC[2, 2]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array q1CC[1, 1]   						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */

		array q4CI[4, 4] 					  	/nosym;	/* inverse of CC */
		array q3CI[3, 3] 					  	/nosym;	/* inverse of CC */
		array q2CI[2, 2] 					  	/nosym;	/* inverse of CC */
		array q1CI[1, 1] 					  	/nosym;	/* inverse of CC */

		/*************************************************************************************/
		/* Array structures for common trend models in the full dimensional cubic BMA: p = 4 */
		/*************************************************************************************/
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array bXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */

		array ambetag[1, 1]						/nosym; /* prior mean vector (1 x 1) for intercepts */
		array b3mbetag[3, 1] 					/nosym;	/* prior mean vector (3 x 1) for remaining coefficients */
		array b2mbetag[2, 1] 					/nosym;	/* prior mean vector (2 x 1) for remaining coefficients */
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for remaining coefficients */

		array aDbetag[1, 1] 					/nosym; /* diagonal matrix (1 x 1) of prior precisions for intercepts */
		array b3Dbetag[3, 3] 					/nosym;	/* diagonal matrix (3 x 3) of prior precisions for remaining coefficients */
		array b2Dbetag[2, 2] 					/nosym;	/* diagonal matrix (2 x 2) of prior precisions for remaining coefficients */
		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for remaining coefficients */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array b3X[&n, 3]						/nosym; /* 3-dimensional conformal design submatrix X */
		array b2X[&n, 2]						/nosym; /* 2-dimensional conformal design submatrix X */
		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */

		array sumb3XtWX[3, 3] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb2XtWX[2, 2] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb3zbeta[3, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb2zbeta[2, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array aXt[1, &n] 						/nosym;	/* transpose of design matrix (intercept only) */
		array b3Xt[3, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b2Xt[2, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array aXtW[1, &n] 						/nosym;	/* matrix multiplication of Xt and Wg (intercept only) */
		array b3XtW[3, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

		array aXtWX[1, 1] 						/nosym;	/* precision matrix of WLS regression estimators (intercept only) */
	 	array b3XtWX[3, 3] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array aDXtWX[1, 1] 						/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (intercept only) */
		array b3DXtWX[3, 3]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b2DXtWX[2, 2]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array aprbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from prior */
		array b3prbeta[3, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b2prbeta[2, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array apbeta[1, 1]	 					/nosym;	/* vector (1x1) of intercepts from pooled posterior */
		array b3pbeta[3, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b2pbeta[2, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array azbeta[1, 1]	 					/nosym;	/* vector of intercepts from WLS */
		array b3zbeta[3, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b2zbeta[2, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array abeta[1, 1] 						/nosym;	/* sampled vector (1x1) of intercepts */
		array b3beta[3, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array b2beta[2, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */
		array b1beta[1, 1]		 	    		/nosym;	/* sampled vector of regression coefficients (excl. intercepts) */

		array aCC[1, 1] 						/nosym;	/* holds lower triangular matrix from Cholesky decomposition */
		array b3CC[3, 3]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array b2CC[2, 2]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */
		array b1CC[1, 1]     					/nosym; /* holds lower triangular matrix from Cholesky decomposition */

		array aCI[1, 1] 						/nosym;	/* inverse of CC */
		array b3CI[3, 3]     					/nosym; /* inverse of CC */
		array b2CI[2, 2]     					/nosym; /* inverse of CC */
		array b1CI[1, 1]     					/nosym; /* inverse of CC */

		/****************************************************************/
		/* Populate needed array structures depending on the model flag */	
		/****************************************************************/
		if flg = 1 then do;								/* indep cubic: q = 4 */
			do i = 1 to &n;
				do m = 1 to 4;
					q4X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q4Dbetag);
			do m = 1 to 4;
				q4mbetag[m, 1] = mbetag[m, 1];
				q4Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q4X, q4Xt);					/* transpose qX */
			call mult(q4Dbetag, q4mbetag, q4prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 2 then do;								/* indep quad: q = 3 */
			do i = 1 to &n;
				do m = 1 to 3;
					q3X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q3Dbetag);
			do m = 1 to 3;
				q3mbetag[m, 1] = mbetag[m, 1];
				q3Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q3X, q3Xt);					/* transpose qX */
			call mult(q3Dbetag, q3mbetag, q3prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 3 then do;								/* indep linear: q = 2 */
			do i = 1 to &n;
				do m = 1 to 2;
					q2X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q2Dbetag);
			do m = 1 to 2;
				q2mbetag[m, 1] = mbetag[m, 1];
				q2Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q2X, q2Xt);					/* transpose qX */
			call mult(q2Dbetag, q2mbetag, q2prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 7 then do;								/* dropped: q = 1 */
			do i = 1 to &n;
				do m = 1 to 1;
					q1X[i, m] = X[i, m];
				end;
			end;
			call zeromatrix(q1Dbetag);	
			do m = 1 to 1;
				q1mbetag[m, 1] = mbetag[m, 1];
				q1Dbetag[m, m] = Dbetag[m, m];
			end;
			call transpose(q1X, q1Xt);					/* transpose qX */
			call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */
		end;
		if flg = 4 then do;								/* common cubic: q = 4 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 4;
					b3X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b3Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 4;
				b3mbetag[m-1, 1]   = mbetag[m, 1];
			    b3Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b3X, b3Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b3Dbetag, b3mbetag, b3prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb3XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb3zbeta);	
		end;
		if flg = 5 then do;								/* common quad: q = 3 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 3;
					b2X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b2Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 3;
				b2mbetag[m-1, 1]   = mbetag[m, 1];
			    b2Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b2X, b2Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b2Dbetag, b2mbetag, b2prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb2XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb2zbeta);	
		end;
		if flg = 6 then do;								/* common linear: q = 2 */
			do i = 1 to &n;								
			  	aX[i, 1] = X[i, 1];
			  	do m = 2 to 2;
					b1X[i, m-1] = X[i, m];
			  	end;
		  	end;
			call zeromatrix(aDbetag);
			call zeromatrix(b1Dbetag); 
		  	ambetag[1,1] = mbetag[1,1];
		  	aDbetag[1,1] = Dbetag[1,1];	
			do m = 2 to 2;
				b1mbetag[m-1, 1]   = mbetag[m, 1];
			    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
			end;
			call transpose(aX, aXt);					/* transpose aX */
			call transpose(b1X, b1Xt);					/* transpose bX */
			call mult(aDbetag, ambetag, aprbeta);		/* contribution to posterior mean from prior */
			call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
			call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
			call zeromatrix(sumb1zbeta);	
		end;

		/*******************************/
		/* Group-specific trend models */
		/*******************************/
		if flg = 1 or flg = 2 or flg = 3 or flg = 7 then do;
			/**********************************/
			/* Update regression coefficients */
			/**********************************/
			do k = 1 to &g;									/* cycle through each group independently */
			    do i = 1 to &n;								
					Yvec[i,1] = Yarr[(k-1)*&n + i]; 		/* populate nx1 data vector Yvec */
				    Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
				    do j = i+1 to &n;
					    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
					    Vg[j,i] = Vg[i,j];
				    end; 
			  	end; 
				call inv(Vg, Wg);							/* Wg = Vg^{-1} */
				if flg = 1 then do;							/* group-specific cubic trend model (q=4)*/
					call mult(q4Xt, Wg, q4XtW);				/* multiply Xt and Wg */
					call mult(q4XtW, q4X, q4XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q4Dbetag,q4XtWX,q4DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q4XtW, Yvec, q4ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q4prbeta,q4ybeta,q4pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 4;
						q4beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q4DXtWX, q4CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q4CC, q4CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q4CI, q4pbeta, q4pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q4CI, q4CI);				/* transpose */
					call mult(q4CI, q4pbeta, q4pbeta);		/* re-scale pbeta (part 2) */
					call mult(q4CI, q4beta, q4beta);		/* re-scale beta */
					call addmatrix(q4pbeta,q4beta,q4beta);	/* re-center */
					call mult(q4X, q4beta, Xbeta);			/* updated vector Xb */
					a[k] = q4beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q4beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = q4beta[3,1];					/* 3rd output argument is 1-d array of group-specific quad coefficients */
					b3arr[k] = q4beta[4,1];					/* 4th output argument is 1-d array of group-specific cubic coefficients */
				end;
				if flg = 2 then do;							/* group-specific quad trend model (q=3)*/
					call mult(q3Xt, Wg, q3XtW);				/* multiply Xt and Wg */
					call mult(q3XtW, q3X, q3XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q3Dbetag,q3XtWX,q3DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q3XtW, Yvec, q3ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q3prbeta,q3ybeta,q3pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 3;
						q3beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q3DXtWX, q3CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q3CC, q3CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q3CI, q3pbeta, q3pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q3CI, q3CI);				/* transpose */
					call mult(q3CI, q3pbeta, q3pbeta);		/* re-scale pbeta (part 2) */
					call mult(q3CI, q3beta, q3beta);		/* re-scale beta */
					call addmatrix(q3pbeta,q3beta,q3beta);	/* re-center */
					call mult(q3X, q3beta, Xbeta);			/* updated vector Xb */
					a[k] = q3beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q3beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = q3beta[3,1];					/* 3rd output argument is 1-d array of group-specific quad coefficients */
					b3arr[k] = 0;							/* 4th output argument is 1-d array of group-specific cubic coefficients */
				end;
				if flg = 3 then do;							/* group-specific linear trend model (q=2)*/
					call mult(q2Xt, Wg, q2XtW);				/* multiply Xt and Wg */
					call mult(q2XtW, q2X, q2XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q2Dbetag,q2XtWX,q2DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q2XtW, Yvec, q2ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q2prbeta,q2ybeta,q2pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 2;
						q2beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q2DXtWX, q2CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q2CC, q2CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q2CI, q2CI);				/* transpose */
					call mult(q2CI, q2pbeta, q2pbeta);		/* re-scale pbeta (part 2) */
					call mult(q2CI, q2beta, q2beta);		/* re-scale beta */
					call addmatrix(q2pbeta,q2beta,q2beta);	/* re-center */
					call mult(q2X, q2beta, Xbeta);			/* updated vector Xb */
					a[k] = q2beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = q2beta[2,1];					/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = 0;							/* 3rd output argument is 1-d array of group-specific quad coefficients */
					b3arr[k] = 0;							/* 4th output argument is 1-d array of group-specific cubic coefficients */
				end;
				if flg = 7 then do;							/* group-specific intercept-only model (q=1)*/
					call mult(q1Xt, Wg, q1XtW);				/* multiply Xt and Wg */
					call mult(q1XtW, q1X, q1XtWX);			/* calculate XtWX, the precision matrix from WLS */
					call addmatrix(q1Dbetag,q1XtWX,q1DXtWX);/* posterior precision matrix for beta is qDbetag + XtWX */
					call mult(q1XtW, Yvec, q1ybeta);		/* contribution to posterior mean from WLS */
					call addmatrix(q1prbeta,q1ybeta,q1pbeta);/* sum of prior and WLS contributions */
					do m = 1 to 1;
						q1beta[m,1] = rand('normal');		/* sample from univariate standard normal */
					end;
					call chol(q1DXtWX, q1CC);				/* Cholesky decomposition for precision matrix (returns lower triangular) */
					call inv(q1CC, q1CI);					/* inverse of lower triangular matrix from Cholesky decomposition */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 1) */
					call transpose(q1CI, q1CI);				/* transpose */
					call mult(q1CI, q1pbeta, q1pbeta);		/* re-scale pbeta (part 2) */
					call mult(q1CI, q1beta, q1beta);		/* re-scale beta */
					call addmatrix(q1pbeta,q1beta,q1beta);	/* re-center */
					call mult(q1X, q1beta, Xbeta);			/* updated vector Xb */
					a[k] = q1beta[1,1];						/* 1st output argument is 1-d array of group-specific intercepts */
					b1arr[k] = 0;							/* 2nd output argument is 1-d array of group-specific linear coefficients */
					b2arr[k] = 0;							/* 3rd output argument is 1-d array of group-specific quad coefficients */
					b3arr[k] = 0;							/* 4th output argument is 1-d array of group-specific cubic coefficients */
				end;
				do i = 1 to &n; 								
			     	etamnarr[(k-1)*&n+i]= Xbeta[i,1]; 		/* updated predictions from regression */
				end;
			end;
			tmpb1 = 0;
			tmpb2 = 0;
			tmpb3 = 0;
			do k = 1 to &g;
				tmpb1 = tmpb1 + b1arr[k]; 		 			/* common coefficient(s) set to average of group-specific coefficients */
				tmpb2 = tmpb2 + b2arr[k];
				tmpb3 = tmpb3 + b3arr[k];
			end;
			b1 = tmpb1/&g;
			b2 = tmpb2/&g;
			b3 = tmpb3/&g;
		end;

		/***********************/
		/* Common trend models */
		/***********************/
		if flg = 4 or flg = 5 or flg = 6 then do;
			/********************************/
			/* Update common coefficient(s) */
			/********************************/
			do k = 1 to &g;										/* cycle through each group independently */
				abeta[1,1] = a[k];								/* group-specific abeta vector */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - aXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - a */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
				if flg = 4 then do;								/* common cubic trend model */
				  call mult(b3Xt, Wg, b3XtW);					/* multiply bXt and Wg */
				  call mult(b3XtW, b3X, b3XtWX);				/* calculate bXtWX */
				  call addmatrix(sumb3XtWX, b3XtWX, sumb3XtWX);	/* cumulative matrix sum */
				  call mult(b3XtW, Zvec, b3zbeta);			 	/* contributions to posterior mean from WLS */
				  call addmatrix(sumb3zbeta,b3zbeta,sumb3zbeta);/* cumulative matrix sum */
				end;
				if flg = 5 then do;								/* common quad trend model */
				  call mult(b2Xt, Wg, b2XtW);					/* multiply bXt and Wg */
				  call mult(b2XtW, b2X, b2XtWX);				/* calculate bXtWX */
				  call addmatrix(sumb2XtWX, b2XtWX, sumb2XtWX);	/* cumulative matrix sum */
				  call mult(b2XtW, Zvec, b2zbeta);			 	/* contributions to posterior mean from WLS */
				  call addmatrix(sumb2zbeta,b2zbeta,sumb2zbeta);/* cumulative matrix sum */
				end;
				if flg = 6 then do;								/* common linear trend model */
				  call mult(b1Xt, Wg, b1XtW);					/* multiply bXt and Wg */
				  call mult(b1XtW, b1X, b1XtWX);				/* calculate bXtWX */
				  call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
				  call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
				  call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);/* cumulative matrix sum */
				end;
			end;												/* end cycle through groups */
			if flg = 4 then do;									/* common cubic trend model */
				call addmatrix(b3Dbetag, sumb3XtWX, b3DXtWX); 	/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
				call addmatrix(b3prbeta, sumb3zbeta, b3pbeta);	/* sum of prior and the cumulative WLS contributions */
				do m = 2 to 4;
					b3beta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
				end;
				call chol(b3DXtWX, b3CC);						/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
				call inv(b3CC, b3CI);							/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(b3CI, b3pbeta, b3pbeta);				/* re-scale pbeta (part 1) */
				call transpose(b3CI, b3CI);						/* transpose */
				call mult(b3CI, b3pbeta, b3pbeta);				/* re-scale pbeta (part 2) */
				call mult(b3CI, b3beta, b3beta);				/* re-scale beta */
				call addmatrix(b3pbeta, b3beta, b3beta);		/* re-center */
				call mult(b3X, b3beta, bXbeta);					/* updated vector bX */
				b1 = b3beta[1,1];								/* output argument b1 is updated value of common linear coefficient */
				b2 = b3beta[2,1];								/* output argument b2 is updated value of common quad coefficient */
				b3 = b3beta[3,1];								/* output argument b3 is updated value of common cubic coefficient */
			end;
			if flg = 5 then do;									/* common quad trend model */
				call addmatrix(b2Dbetag, sumb2XtWX, b2DXtWX); 	/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
				call addmatrix(b2prbeta, sumb2zbeta, b2pbeta);	/* sum of prior and the cumulative WLS contributions */
				do m = 2 to 3;
					b2beta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
				end;
				call chol(b2DXtWX, b2CC);						/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
				call inv(b2CC, b2CI);							/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(b2CI, b2pbeta, b2pbeta);				/* re-scale pbeta (part 1) */
				call transpose(b2CI, b2CI);						/* transpose */
				call mult(b2CI, b2pbeta, b2pbeta);				/* re-scale pbeta (part 2) */
				call mult(b2CI, b2beta, b2beta);				/* re-scale beta */
				call addmatrix(b2pbeta, b2beta, b2beta);		/* re-center */
				call mult(b2X, b2beta, bXbeta);					/* updated vector bX */
				b1 = b2beta[1,1];								/* output argument b1 is updated value of common linear coefficient */
				b2 = b2beta[2,1];								/* output argument b2 is updated value of common quad coefficient */
				b3 = 0;											/* output argument b3 is updated value of common cubic coefficient */
			end;
			if flg = 6 then do;									/* common linear trend model */
				call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 	/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
				call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);	/* sum of prior and the cumulative WLS contributions */
				do m = 2 to 2;
					b1beta[m-1,1] = rand('normal');				/* sample from univariate standard normal(s) */
				end;
				call chol(b1DXtWX, b1CC);						/* Cholesky decomposition for (p-1)x(p-1) precision matrix (returns lower triangular) */
				call inv(b1CC, b1CI);							/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(b1CI, b1pbeta, b1pbeta);				/* re-scale pbeta (part 1) */
				call transpose(b1CI, b1CI);						/* transpose */
				call mult(b1CI, b1pbeta, b1pbeta);				/* re-scale pbeta (part 2) */
				call mult(b1CI, b1beta, b1beta);				/* re-scale beta */
				call addmatrix(b1pbeta, b1beta, b1beta);		/* re-center */
				call mult(b1X, b1beta, bXbeta);					/* updated vector bX */
				b1 = b1beta[1,1];								/* output argument b1 is updated value of common linear coefficient */
				b2 = 0;											/* output argument b2 is updated value of common quad coefficient */
				b3 = 0;											/* output argument b3 is updated value of common cubic coefficient */
			end;
			do k = 1 to &g;										/* arrays of group-specific coefficient also updated to reflect common value */
				b1arr[k] = b1;
				b2arr[k] = b2;
				b3arr[k] = b3;
			end;
			/************************************************/
			/* Update intercepts and regression predictions */
			/************************************************/
			do k = 1 to &g;
			    do i = 1 to &n;						
				  Zvec[i,1]= Yarr[(k-1)*&n+i] - bXbeta[i,1];	/* populate nx1 data vector Zvec = Yvec - bX */
				  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  	/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
				end; 
				do i = 1 to &n-1;
				  do j = i+1 to &n;
					Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
					Vg[j,i] = Vg[i,j];
				  end; 
			  	end; 
				call inv(Vg, Wg);								/* Wg = Vg^{-1} */
				call mult(aXt, Wg, aXtW);						/* multiply aXt and Wg */
				call mult(aXtW, aX, aXtWX);						/* calculate aXtWX, the precision matrix from WLS */
				call addmatrix(aDbetag, aXtWX, aDXtWX); 		/* posterior precision matrix is aDbetag + XtWX */
				call mult(aXtW, Zvec, azbeta);					/* contribution to posterior mean from WLS */
				call addmatrix(aprbeta, azbeta, apbeta); 		/* sum of prior and WLS contributions */
				abeta[1,1] = rand('normal');					/* sample intercept from univariate normal */
				call chol(aDXtWX, aCC);							/* Cholesky decomposition for precision matrix (returns lower triangular) */
				call inv(aCC, aCI);								/* inverse of lower triangular matrix from Cholesky decomposition */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 1) */
				call transpose(aCI, aCI);						/* transpose */
				call mult(aCI, apbeta, apbeta);					/* re-scale pbeta (part 2) */
				call mult(aCI, abeta, abeta);					/* re-scale beta */
				call addmatrix(apbeta, abeta, abeta);			/* re-center */
				a[k] = abeta[1,1];								/* output argument a is 1-dimensional array of group-specific intercepts */
				call mult(aX, abeta, aXbeta);					/* vector of intercepts a */
				do i = 1 to &n; 								
			      etamnarr[(k-1)*&n+i] = aXbeta[i,1] + bXbeta[i,1]; /* updated predictions from regression */
				end;
			end;
		end;

		endsub;
	run;
	quit;

%end;

%mend;

data _null_;
run;

/**********************************************************************/
/* eMKF: Gibbs samplers for model flags in the supported trend models */
/**********************************************************************/
%macro gibbs_uds_compile_FP(uvar=, g=, n=, loc=);

/* eMKF: return if no applicable model is indicated */
%if %upcase(&uvar) ^= BMA_CUBIC and %upcase(&uvar) ^= BMA_QUAD and %upcase(&uvar) ^= BMA_LINEAR
	%then %do;
		%put ERROR: No Gibbs sampler for model flag was found for the specified Bayesian model averaging &uvar: Please check!;
		proc iml;
		  print "  Error Note:";
		  print "  No Gibbs sampler for model flag was found for the specified Bayesian model averaging. Please check. ";
	    quit;
	%return;
%end;

%local p uloc;

/* eMKF: dimensionality (needed for UDS set up) */
%let p = 0;
%if %upcase(&uvar) = BMA_LINEAR %then %let p = 2;
%if %upcase(&uvar) = BMA_QUAD   %then %let p = 3;
%if %upcase(&uvar) = BMA_CUBIC  %then %let p = 4;
%let p = %eval(0+&p);

%let uloc = &loc..uds;

%if %upcase(&uvar) = BMA_LINEAR %then %do;
	/********************************************************************/
	/* eMKF: Gibbs sampler for model flag in the BMA linear trend model */
	/********************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine FP_bmal(flg,					/* model flag , with integer values 1 through 1+(p-1)*2 */
						   wts[*],				/* 1-dimensional array (length 1+(p-1)*2) of prior model probabilities */
						   a[*],				/* 1-dimensional array of intercepts to subtract off from y's */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						   );

		outargs flg;							/* argument that is updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array pwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */
		array qwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array abeta[1, 1] 						/nosym;	/* vector (1x1) of intercepts */

		/******************************************/
		/* Array structures for indep trend model */
		/******************************************/
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */

		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */

		array q1tmbetag[1, 1]					/nosym; /* transpose */

		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */

		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */

		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */

		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */

		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */

		array q1tpbeta[1, 1] 	    			/nosym;	/* transpose */

		array q1zbeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */

		array q1CI[1, 1] 					  	/nosym;	/* inverse */

		array q1exp[1, 1]						/nosym; /* holds quadratic form */

		/*******************************************/
		/* Array structures for common trend model */
		/*******************************************/
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for coefficients (excl. intercept) */

		array b1tmbetag[1, 1] 					/nosym;	/* transpose */

		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for coefficients (excl. intercept) */

		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */

		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array b1tpbeta[1, 1]            		/nosym; /* transpose */

		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array b1CI[1, 1]     					/nosym; /* inverse */

		array b1exp[1,1]						/nosym; /* holds quadratic form */

		/*******************************************************************/
		/* Populate required array structures for each possible model flag */	
		/*******************************************************************/

		/*****************************************/
		/* indep linear with no intercept: q = 1 */
		/*****************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 2;
				q1X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q1Dbetag);
		do m = 2 to 2;
			q1mbetag[m-1, 1]   = mbetag[m, 1];
			q1Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q1X, q1Xt);					/* transpose qX */
		call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */

		/******************************************/
		/* common linear with no intercept: q = 1 */
		/******************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 2;
				b1X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b1Dbetag); 
		do m = 2 to 2;
			b1mbetag[m-1, 1]   = mbetag[m, 1];
		    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b1X, b1Xt);					/* transpose bX */
		call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb1zbeta);	

		/**************************************************/
		/* Group-specific trend models with no intercepts */
		/**************************************************/

		pwts[3] = 0;									/* no contribution from marginal for group-specific intercept-only model */

		pwts[1] = 0;									/* group-specific linear trend model */
		call transpose(q1mbetag, q1tmbetag);
		call mult(q1tmbetag, q1prbeta, q1exp);			/* exponent from normal pdf features &g times */
		pwts[1] = pwts[1] - 0.5*&g*q1exp[1,1];			/* log scale */
		call det(q1Dbetag, q1wts);						/* determinant of prior precision matrix features &g times */
		pwts[1] = pwts[1] + 0.5*&g*log(q1wts);			/* log scale */
	
		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;								
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 		/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
			    do j = i+1 to &n;
				    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				    Vg[j,i] = Vg[i,j];
			    end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */

			/*************************************************************/
			/* group-specific linear trend model with no intercept (q=1) */
			/*************************************************************/
			call mult(q1Xt, Wg, q1XtW);					/* multiply Xt and Wg */
			call mult(q1XtW, q1X, q1XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q1Dbetag, q1XtWX, q1DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q1XtW, Zvec, q1zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q1prbeta, q1zbeta, q1pbeta);	/* sum of prior and WLS contributions */
			call transpose(q1pbeta, q1tpbeta);			/* transpose */
			call inv(q1DXtWX, q1CI);					/* inverse of posterior precision matrix */
			call mult(q1CI, q1pbeta, q1pbeta);			/* re-scale pbeta */
			call mult(q1tpbeta, q1pbeta, q1exp);		/* exponent from normal pdf */
			pwts[1] = pwts[1] + 0.5*q1exp[1,1];			/* log scale */
			call det(q1DXtWX, q1wts);
			pwts[1] = pwts[1] - 0.5*log(q1wts); 		/* determinant for normalizing constant */

		end;

		/******************************************/
		/* Common trend models with no intercepts */
		/******************************************/

		pwts[2] = 0;									/* common linear trend model */
		call transpose(b1mbetag, b1tmbetag);
		call mult(b1tmbetag, b1prbeta, b1exp);			/* exponent from normal pdf */
		pwts[2] = pwts[2] - 0.5*b1exp[1,1];				/* log scale */
		call det(b1Dbetag, q1wts);						/* determinant of prior precision matrix */
		pwts[2] = pwts[2] + 0.5*log(q1wts);				/* log scale */				

		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);								/* Wg = Vg^{-1} */

			/************************************************/
			/* common linear trend model with no intercepts */
			/************************************************/
		    call mult(b1Xt, Wg, b1XtW);						/* multiply bXt and Wg */
		    call mult(b1XtW, b1X, b1XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
		    call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);	/* cumulative matrix sum */

		end;												/* end cycle through groups */

		/************************************************/
		/* common linear trend model with no intercepts */
		/************************************************/
		call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b1pbeta, b1tpbeta);					/* transpose */
		call inv(b1DXtWX, b1CI);							/* inverse of posterior precision matrix  */
		call mult(b1CI, b1pbeta, b1pbeta);					/* re-scale pbeta */
		call mult(b1tpbeta, b1pbeta, b1exp);				/* exponent from normal pdf */
		pwts[2] = pwts[2] + 0.5*b1exp[1,1];					/* log scale */
		call det(b1DXtWX, b1wts);
		pwts[2] = pwts[2] - 0.5*log(b1wts); 				/* determinant for normalizing constant */

		/************************************/
		/* Posterior model weights and draw */
		/************************************/
		do m = 1 to (1+(&p-1)*2);							/* calculate using differences on log-scale (Bayes factors) for numerical stability */
			qwtsum = 0;
			do l = 1 to (1+(&p-1)*2);
				qwtsum = qwtsum + wts[l]*exp(pwts[l]-pwts[m]);
			end;
			qwts[m] = wts[m]/qwtsum;
		end;
		qwtsum = .;											/* check for any missing values */
		do m = 1 to (1+(&p-1)*2);
			if qwts[m] = . then qwtsum = 0;					/* set qwtsum to zero if any missings */
		end;
		if qwtsum = 0 then do;								/* replace missing values with zeroes and cumulate sum (if applicable) */
			do m = 1 to (1+(&p-1)*2);
				if qwts[m] = . then qwts[m] = 0;
				qwtsum = qwtsum + qwts[m];
			end;
		end;
		if qwtsum > 0 then do;								/* rescale to sum to 1 (if applicable) */
			do m = 1 to (1+(&p-1)*2);						
				qwts[m] = qwts[m]/qwtsum;
			end;
		end;
		if qwtsum = 0 then do;								/* in case all weights are zero, use prior weights (if applicable) */
			qwts[m] = wts[m];
		end;

		flg = rand('table', qwts[1], qwts[2], qwts[3]);		/* returned draw from conditional posterior for model flag */
		
		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = BMA_QUAD %then %do;
	/********************************************************************/
	/* eMKF: Gibbs sampler for model flag in the BMA quad trend model */
	/********************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine FP_bmaq(flg,					/* model flag , with integer values 1 through 1+(p-1)*2 */
						   wts[*],				/* 1-dimensional array (length 1+(p-1)*2) of prior model probabilities */
						   a[*],				/* 1-dimensional array of intercepts to subtract off from y's */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						   );

		outargs flg;							/* argument that is updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array pwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */
		array qwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array abeta[1, 1] 						/nosym;	/* vector (1x1) of intercepts */

		/*******************************************/
		/* Array structures for indep trend models */
		/*******************************************/
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */
		array q2X[&n, 2]						/nosym; /* 2-column version of the design matrix X */

		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */
		array q2mbetag[2, 1]					/nosym; /* 2-dimensional version of mbetag */

		array q1tmbetag[1, 1]					/nosym; /* transpose */
		array q2tmbetag[1, 2]					/nosym; /* transpose */

		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */
		array q2Dbetag[2, 2]					/nosym; /* 2-dimensional version of Dbetag */

		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */
		array q2Xt[2, &n]   					/nosym;	/* transpose of design matrix */

		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */
		array q2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators */

		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q2DXtWX[2, 2]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */
		array q2prbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from prior */

		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */
		array q2pbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from pooled posterior */

		array q1tpbeta[1, 1] 	    			/nosym;	/* transpose */
		array q2tpbeta[1, 2] 	    			/nosym;	/* transpose */

		array q1zbeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */
		array q2zbeta[2, 1] 	       			/nosym;	/* vector (2 x 1) of regression estimates from WLS */

		array q1CI[1, 1] 					  	/nosym;	/* inverse */
		array q2CI[2, 2] 					  	/nosym;	/* inverse */

		array q1exp[1, 1]						/nosym; /* holds quadratic form */
		array q2exp[1, 1]						/nosym; /* holds quadratic form */

		/********************************************/
		/* Array structures for common trend models */
		/********************************************/
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for coefficients (excl. intercept) */
		array b2mbetag[2, 1] 					/nosym;	/* prior mean vector (2 x 1) for coefficients (excl. intercept) */

		array b1tmbetag[1, 1] 					/nosym;	/* transpose */
		array b2tmbetag[1, 2] 					/nosym;	/* transpose */

		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for coefficients (excl. intercept) */
		array b2Dbetag[2, 2] 					/nosym;	/* diagonal matrix (2 x 2) of prior precisions for coefficients (excl. intercept) */

		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */
		array b2X[&n, 2]						/nosym; /* 2-dimensional conformal design submatrix X */

		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb2XtWX[2, 2] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb2zbeta[2, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b2Xt[2, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b2DXtWX[2, 2]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b2prbeta[2, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b2pbeta[2, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array b1tpbeta[1, 1]            		/nosym; /* transpose */
		array b2tpbeta[1, 2]            		/nosym; /* transpose */

		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b2zbeta[2, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array b1CI[1, 1]     					/nosym; /* inverse */
		array b2CI[2, 2]     					/nosym; /* inverse */

		array b1exp[1,1]						/nosym; /* holds quadratic form */
		array b2exp[1,1]						/nosym; /* holds quadratic form */

		/*******************************************************************/
		/* Populate required array structures for each possible model flag */	
		/*******************************************************************/

		/***************************************/
		/* indep quad with no intercept: q = 2 */
		/***************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 3;
				q2X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q2Dbetag);
		do m = 2 to 3;
			q2mbetag[m-1, 1]   = mbetag[m, 1];
			q2Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q2X, q2Xt);					/* transpose qX */
		call mult(q2Dbetag, q2mbetag, q2prbeta);	/* contribution to posterior mean from prior */

		/*****************************************/
		/* indep linear with no intercept: q = 1 */
		/*****************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 2;
				q1X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q1Dbetag);
		do m = 2 to 2;
			q1mbetag[m-1, 1]   = mbetag[m, 1];
			q1Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q1X, q1Xt);					/* transpose qX */
		call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */

		/****************************************/
		/* common quad with no intercept: q = 2 */
		/****************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 3;
				b2X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b2Dbetag); 
		do m = 2 to 3;
			b2mbetag[m-1, 1]   = mbetag[m, 1];
		    b2Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b2X, b2Xt);					/* transpose bX */
		call mult(b2Dbetag, b2mbetag, b2prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb2XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb2zbeta);	

		/******************************************/
		/* common linear with no intercept: q = 1 */
		/******************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 2;
				b1X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b1Dbetag); 
		do m = 2 to 2;
			b1mbetag[m-1, 1]   = mbetag[m, 1];
		    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b1X, b1Xt);					/* transpose bX */
		call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb1zbeta);	

		/**************************************************/
		/* Group-specific trend models with no intercepts */
		/**************************************************/

		pwts[5] = 0;									/* no contribution from marginal for group-specific intercept-only model */

		pwts[1] = 0;									/* group-specific quad trend model */
		call transpose(q2mbetag, q2tmbetag);
		call mult(q2tmbetag, q2prbeta, q2exp);			/* exponent from normal pdf features &g times */
		pwts[1] = pwts[1] - 0.5*&g*q2exp[1,1];			/* log scale */
		call det(q2Dbetag, q2wts);						/* determinant of prior precision matrix features &g times */
		pwts[1] = pwts[1] + 0.5*&g*log(q2wts);			/* log scale */

		pwts[2] = 0;									/* group-specific linear trend model */
		call transpose(q1mbetag, q1tmbetag);
		call mult(q1tmbetag, q1prbeta, q1exp);			/* exponent from normal pdf features &g times */
		pwts[1] = pwts[1] - 0.5*&g*q1exp[1,1];			/* log scale */
		call det(q1Dbetag, q1wts);						/* determinant of prior precision matrix features &g times */
		pwts[1] = pwts[1] + 0.5*&g*log(q1wts);			/* log scale */
	
		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;								
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 		/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
			    do j = i+1 to &n;
				    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				    Vg[j,i] = Vg[i,j];
			    end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */

			/***********************************************************/
			/* group-specific quad trend model with no intercept (q=2) */
			/***********************************************************/
			call mult(q2Xt, Wg, q2XtW);					/* multiply Xt and Wg */
			call mult(q2XtW, q2X, q2XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q2Dbetag, q2XtWX, q2DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q2XtW, Zvec, q2zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q2prbeta, q2zbeta, q2pbeta);	/* sum of prior and WLS contributions */
			call transpose(q2pbeta, q2tpbeta);			/* transpose */
			call inv(q2DXtWX, q2CI);					/* inverse of posterior precision matrix */
			call mult(q2CI, q2pbeta, q2pbeta);			/* re-scale pbeta */
			call mult(q2tpbeta, q2pbeta, q2exp);		/* exponent from normal pdf */
			pwts[1] = pwts[1] + 0.5*q2exp[1,1];			/* log scale */
			call det(q2DXtWX, q2wts);
			pwts[1] = pwts[1] - 0.5*log(q2wts); 		/* determinant for normalizing constant */

			/*************************************************************/
			/* group-specific linear trend model with no intercept (q=1) */
			/*************************************************************/
			call mult(q1Xt, Wg, q1XtW);					/* multiply Xt and Wg */
			call mult(q1XtW, q1X, q1XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q1Dbetag, q1XtWX, q1DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q1XtW, Zvec, q1zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q1prbeta, q1zbeta, q1pbeta);	/* sum of prior and WLS contributions */
			call transpose(q1pbeta, q1tpbeta);			/* transpose */
			call inv(q1DXtWX, q1CI);					/* inverse of posterior precision matrix */
			call mult(q1CI, q1pbeta, q1pbeta);			/* re-scale pbeta */
			call mult(q1tpbeta, q1pbeta, q1exp);		/* exponent from normal pdf */
			pwts[2] = pwts[2] + 0.5*q1exp[1,1];			/* log scale */
			call det(q1DXtWX, q1wts);
			pwts[2] = pwts[2] - 0.5*log(q1wts); 		/* determinant for normalizing constant */

		end;

		/******************************************/
		/* Common trend models with no intercepts */
		/******************************************/

		pwts[3] = 0;									/* common quad trend model */
		call transpose(b2mbetag, b2tmbetag);
		call mult(b2tmbetag, b2prbeta, b2exp);			/* exponent from normal pdf */
		pwts[3] = pwts[3] - 0.5*b2exp[1,1];				/* log scale */
		call det(b2Dbetag, q2wts);						/* determinant of prior precision matrix */
		pwts[3] = pwts[3] + 0.5*log(q2wts);				/* log scale */		
	
		pwts[4] = 0;									/* common linear trend model */
		call transpose(b1mbetag, b1tmbetag);
		call mult(b1tmbetag, b1prbeta, b1exp);			/* exponent from normal pdf */
		pwts[4] = pwts[4] - 0.5*b1exp[1,1];				/* log scale */
		call det(b1Dbetag, q1wts);						/* determinant of prior precision matrix */
		pwts[4] = pwts[4] + 0.5*log(q1wts);				/* log scale */				

		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);								/* Wg = Vg^{-1} */

			/**********************************************/
			/* common quad trend model with no intercepts */
			/**********************************************/
		    call mult(b2Xt, Wg, b2XtW);						/* multiply bXt and Wg */
		    call mult(b2XtW, b2X, b2XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb2XtWX, b2XtWX, sumb2XtWX);	/* cumulative matrix sum */
		    call mult(b2XtW, Zvec, b2zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb2zbeta,b2zbeta,sumb2zbeta);	/* cumulative matrix sum */

			/************************************************/
			/* common linear trend model with no intercepts */
			/************************************************/
		    call mult(b1Xt, Wg, b1XtW);						/* multiply bXt and Wg */
		    call mult(b1XtW, b1X, b1XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
		    call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);	/* cumulative matrix sum */

		end;												/* end cycle through groups */

		/**********************************************/
		/* common quad trend model with no intercepts */
		/**********************************************/
		call addmatrix(b2Dbetag, sumb2XtWX, b2DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b2prbeta, sumb2zbeta, b2pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b2pbeta, b2tpbeta);					/* transpose */
		call inv(b2DXtWX, b2CI);							/* inverse of posterior precision matrix  */
		call mult(b2CI, b2pbeta, b2pbeta);					/* re-scale pbeta */
		call mult(b2tpbeta, b2pbeta, b2exp);				/* exponent from normal pdf */
		pwts[3] = pwts[3] + 0.5*b2exp[1,1];					/* log scale */
		call det(b2DXtWX, b2wts);
		pwts[3] = pwts[3] - 0.5*log(b2wts); 				/* determinant for normalizing constant */

		/************************************************/
		/* common linear trend model with no intercepts */
		/************************************************/
		call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b1pbeta, b1tpbeta);					/* transpose */
		call inv(b1DXtWX, b1CI);							/* inverse of posterior precision matrix  */
		call mult(b1CI, b1pbeta, b1pbeta);					/* re-scale pbeta */
		call mult(b1tpbeta, b1pbeta, b1exp);				/* exponent from normal pdf */
		pwts[4] = pwts[4] + 0.5*b1exp[1,1];					/* log scale */
		call det(b1DXtWX, b1wts);
		pwts[4] = pwts[4] - 0.5*log(b1wts); 				/* determinant for normalizing constant */

		/************************************/
		/* Posterior model weights and draw */
		/************************************/
		do m = 1 to (1+(&p-1)*2);							/* calculate using differences on log-scale (Bayes factors) for numerical stability */
			qwtsum = 0;
			do l = 1 to (1+(&p-1)*2);
				qwtsum = qwtsum + wts[l]*exp(pwts[l]-pwts[m]);
			end;
			qwts[m] = wts[m]/qwtsum;
		end;
		qwtsum = .;											/* check for any missing values */
		do m = 1 to (1+(&p-1)*2);
			if qwts[m] = . then qwtsum = 0;					/* set qwtsum to zero if any missings */
		end;
		if qwtsum = 0 then do;								/* replace missing values with zeroes and cumulate sum (if applicable) */
			do m = 1 to (1+(&p-1)*2);
				if qwts[m] = . then qwts[m] = 0;
				qwtsum = qwtsum + qwts[m];
			end;
		end;
		if qwtsum > 0 then do;								/* rescale to sum to 1 (if applicable) */
			do m = 1 to (1+(&p-1)*2);						
				qwts[m] = qwts[m]/qwtsum;
			end;
		end;
		if qwtsum = 0 then do;								/* in case all weights are zero, use prior weights (if applicable) */
			qwts[m] = wts[m];
		end;

		flg = rand('table', qwts[1], qwts[2], qwts[3], 
									 qwts[4], qwts[5]);		/* returned draw from conditional posterior for model flag */
		
		endsub;
	run;
	quit;

%end;

%if %upcase(&uvar) = BMA_CUBIC %then %do;
	/********************************************************************/
	/* eMKF: Gibbs sampler for model flag in the BMA cubic trend model  */
	/********************************************************************/
	proc fcmp outlib=&uloc; 			

		subroutine FP_bmac(flg,					/* model flag , with integer values 1 through 1+(p-1)*2 */
						   wts[*],				/* 1-dimensional array (length 1+(p-1)*2) of prior model probabilities */
						   a[*],				/* 1-dimensional array of intercepts to subtract off from y's */
						   mbetag[*,*], 		/* prior mean vector (p x 1) for regression coefficients */
						   Dbetag[*,*], 		/* diagonal matrix (p x p) of prior precisions for regression coefficients */
						   rhoarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances rho */
						   nuarr[*],			/* 1-dimensional array (length g) of current values of group-specific AR variances nu */
						   rts[*],				/* 1-dimensional array (length n) of real times */
						   X[*,*], 				/* design matrix (n x p) using real times */
						   Yarr[*], 			/* 1-dimensional array (length gn) for _y from dataset */
						   Sarr[*]				/* 1-dimensional array (length gn) for _var from dataset */
						   );

		outargs flg;							/* argument that is updated after execution */

		/****************************/
		/* General array structures */
		/****************************/
		array pwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */
		array qwts[%eval(1+(&p-1)*2)]			/nosym; /* holds posterior weights for model flags */

		array Zvec[&n, 1]		 				/nosym;	/* de-trended group-specific observations */
		array Vg[&n, &n]  						/nosym;	/* Vgamma + sampling variances */
		array Wg[&n, &n]   						/nosym;	/* (Vgamma + sampling variances)^{-1} */

		array aX[&n, 1]							/nosym; /* 1-dimensional conformal design submatrix X */
		array aXbeta[&n, 1]						/nosym;	/* holds matrix multiplication */
		array abeta[1, 1] 						/nosym;	/* vector (1x1) of intercepts */

		/*******************************************/
		/* Array structures for indep trend models */
		/*******************************************/
		array q1X[&n, 1]						/nosym; /* 1-column version of the design matrix X */
		array q2X[&n, 2]						/nosym; /* 2-column version of the design matrix X */
		array q3X[&n, 3]						/nosym; /* 3-column version of the design matrix X */

		array q1mbetag[1, 1]					/nosym; /* 1-dimensional version of mbetag */
		array q2mbetag[2, 1]					/nosym; /* 2-dimensional version of mbetag */
		array q3mbetag[3, 1]					/nosym; /* 3-dimensional version of mbetag */

		array q1tmbetag[1, 1]					/nosym; /* transpose */
		array q2tmbetag[1, 2]					/nosym; /* transpose */
		array q3tmbetag[1, 3]					/nosym; /* transpose */

		array q1Dbetag[1, 1]					/nosym; /* 1-dimensional version of Dbetag */
		array q2Dbetag[2, 2]					/nosym; /* 2-dimensional version of Dbetag */
		array q3Dbetag[3, 3]					/nosym; /* 3-dimensional version of Dbetag */

		array q1Xt[1, &n]   					/nosym;	/* transpose of design matrix */
		array q2Xt[2, &n]   					/nosym;	/* transpose of design matrix */
		array q3Xt[3, &n]   					/nosym;	/* transpose of design matrix */

		array q1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg */
		array q3XtW[3, &n]						/nosym; /* matrix multiplication of Xt and Wg */

		array q1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators */
		array q2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators */
		array q3XtWX[3, 3] 						/nosym; /* precision matrix of WLS regression estimators */

		array q1DXtWX[1, 1]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q2DXtWX[2, 2]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */
		array q3DXtWX[3, 3]						/nosym; /* Dbetag + XtWX = posterior precision matrix for beta */

		array q1prbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from prior */
		array q2prbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from prior */
		array q3prbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from prior */

		array q1pbeta[1, 1] 	    			/nosym;	/* vector (1 x 1) of regression estimates from pooled posterior */
		array q2pbeta[2, 1] 	    			/nosym;	/* vector (2 x 1) of regression estimates from pooled posterior */
		array q3pbeta[3, 1] 	    			/nosym;	/* vector (3 x 1) of regression estimates from pooled posterior */

		array q1tpbeta[1, 1] 	    			/nosym;	/* transpose */
		array q2tpbeta[1, 2] 	    			/nosym;	/* transpose */
		array q3tpbeta[1, 3] 	    			/nosym;	/* transpose */

		array q1zbeta[1, 1] 	       			/nosym;	/* vector (1 x 1) of regression estimates from WLS */
		array q2zbeta[2, 1] 	       			/nosym;	/* vector (2 x 1) of regression estimates from WLS */
		array q3zbeta[3, 1] 	       			/nosym;	/* vector (3 x 1) of regression estimates from WLS */

		array q1CI[1, 1] 					  	/nosym;	/* inverse */
		array q2CI[2, 2] 					  	/nosym;	/* inverse */
		array q3CI[3, 3] 					  	/nosym;	/* inverse */

		array q1exp[1, 1]						/nosym; /* holds quadratic form */
		array q2exp[1, 1]						/nosym; /* holds quadratic form */
		array q3exp[1, 1]						/nosym; /* holds quadratic form */

		/********************************************/
		/* Array structures for common trend models */
		/********************************************/
		array b1mbetag[1, 1] 					/nosym;	/* prior mean vector (1 x 1) for coefficients (excl. intercept) */
		array b2mbetag[2, 1] 					/nosym;	/* prior mean vector (2 x 1) for coefficients (excl. intercept) */
		array b3mbetag[3, 1] 					/nosym;	/* prior mean vector (3 x 1) for coefficients (excl. intercept) */

		array b1tmbetag[1, 1] 					/nosym;	/* transpose */
		array b2tmbetag[1, 2] 					/nosym;	/* transpose */
		array b3tmbetag[1, 3] 					/nosym;	/* transpose */

		array b1Dbetag[1, 1] 					/nosym;	/* diagonal matrix (1 x 1) of prior precisions for coefficients (excl. intercept) */
		array b2Dbetag[2, 2] 					/nosym;	/* diagonal matrix (2 x 2) of prior precisions for coefficients (excl. intercept) */
		array b3Dbetag[3, 3] 					/nosym;	/* diagonal matrix (3 x 3) of prior precisions for coefficients (excl. intercept) */

		array b1X[&n, 1]						/nosym; /* 1-dimensional conformal design submatrix X */
		array b2X[&n, 2]						/nosym; /* 2-dimensional conformal design submatrix X */
		array b3X[&n, 3]						/nosym; /* 3-dimensional conformal design submatrix X */

		array sumb1XtWX[1, 1] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb2XtWX[2, 2] 					/nosym;	/* cumulative sum of group-specific precision matrices */
		array sumb3XtWX[3, 3] 					/nosym;	/* cumulative sum of group-specific precision matrices */

		array sumb1zbeta[1, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb2zbeta[2, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */
		array sumb3zbeta[3, 1] 	    			/nosym;	/* cumulative sum of group-specific vector of regression estimates */

		array b1Xt[1, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b2Xt[2, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */
		array b3Xt[3, &n] 						/nosym; /* transpose of design matrix (excl. intercept) */

		array b1XtW[1, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b2XtW[2, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */
		array b3XtW[3, &n]						/nosym; /* matrix multiplication of Xt and Wg (excl. intercept) */

	 	array b1XtWX[1, 1] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b2XtWX[2, 2] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */
	 	array b3XtWX[3, 3] 						/nosym; /* precision matrix of WLS regression estimators (excl. intercept) */

		array b1DXtWX[1, 1]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b2DXtWX[2, 2]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */
		array b3DXtWX[3, 3]  					/nosym;	/* Dbetag + XtWX = posterior precision matrix for beta (excl. intercept) */

		array b1prbeta[1, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b2prbeta[2, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */
		array b3prbeta[3, 1]            		/nosym; /* vector of regression estimates (excl. intercept) from prior */

		array b1pbeta[1, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b2pbeta[2, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */
		array b3pbeta[3, 1]            			/nosym; /* vector of regression estimates (excl. intercept) from pooled posterior */

		array b1tpbeta[1, 1]            		/nosym; /* transpose */
		array b2tpbeta[1, 2]            		/nosym; /* transpose */
		array b3tpbeta[1, 3]            		/nosym; /* transpose */

		array b1zbeta[1, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b2zbeta[2, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */
		array b3zbeta[3, 1]	 					/nosym;	/* vector of regression estimates (excl. intercepts) from WLS */

		array b1CI[1, 1]     					/nosym; /* inverse */
		array b2CI[2, 2]     					/nosym; /* inverse */
		array b3CI[3, 3]     					/nosym; /* inverse */

		array b1exp[1,1]						/nosym; /* holds quadratic form */
		array b2exp[1,1]						/nosym; /* holds quadratic form */
		array b3exp[1,1]						/nosym; /* holds quadratic form */

		/*******************************************************************/
		/* Populate required array structures for each possible model flag */	
		/*******************************************************************/

		/****************************************/
		/* indep cubic with no intercept: q = 3 */
		/****************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 4;
				q3X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q3Dbetag);
		do m = 2 to 4;
			q3mbetag[m-1, 1]   = mbetag[m, 1];
			q3Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q3X, q3Xt);					/* transpose qX */
		call mult(q3Dbetag, q3mbetag, q3prbeta);	/* contribution to posterior mean from prior */

		/***************************************/
		/* indep quad with no intercept: q = 2 */
		/***************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 3;
				q2X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q2Dbetag);
		do m = 2 to 3;
			q2mbetag[m-1, 1]   = mbetag[m, 1];
			q2Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q2X, q2Xt);					/* transpose qX */
		call mult(q2Dbetag, q2mbetag, q2prbeta);	/* contribution to posterior mean from prior */

		/*****************************************/
		/* indep linear with no intercept: q = 1 */
		/*****************************************/
		do i = 1 to &n;
		  	aX[i, 1] = X[i, 1];
			do m = 2 to 2;
				q1X[i, m-1] = X[i, m];
			end;
		end;
		call zeromatrix(q1Dbetag);
		do m = 2 to 2;
			q1mbetag[m-1, 1]   = mbetag[m, 1];
			q1Dbetag[m-1, m-1] = Dbetag[m, m];
		end;
		call transpose(q1X, q1Xt);					/* transpose qX */
		call mult(q1Dbetag, q1mbetag, q1prbeta);	/* contribution to posterior mean from prior */

		/*****************************************/
		/* common cubic with no intercept: q = 3 */
		/*****************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 4;
				b3X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b3Dbetag); 
		do m = 2 to 4;
			b3mbetag[m-1, 1]   = mbetag[m, 1];
		    b3Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b3X, b3Xt);					/* transpose bX */
		call mult(b3Dbetag, b3mbetag, b3prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb3XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb3zbeta);	

		/****************************************/
		/* common quad with no intercept: q = 2 */
		/****************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 3;
				b2X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b2Dbetag); 
		do m = 2 to 3;
			b2mbetag[m-1, 1]   = mbetag[m, 1];
		    b2Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b2X, b2Xt);					/* transpose bX */
		call mult(b2Dbetag, b2mbetag, b2prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb2XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb2zbeta);	

		/******************************************/
		/* common linear with no intercept: q = 1 */
		/******************************************/
		do i = 1 to &n;								
		  	aX[i, 1] = X[i, 1];
		  	do m = 2 to 2;
				b1X[i, m-1] = X[i, m];
		  	end;
	  	end;
		call zeromatrix(b1Dbetag); 
		do m = 2 to 2;
			b1mbetag[m-1, 1]   = mbetag[m, 1];
		    b1Dbetag[m-1, m-1] = Dbetag[m, m];	
		end;
		call transpose(b1X, b1Xt);					/* transpose bX */
		call mult(b1Dbetag, b1mbetag, b1prbeta);	/* contribution to posterior mean from prior */
		call zeromatrix(sumb1XtWX);					/* initialize applicable cumulative sums to all zeroes */
		call zeromatrix(sumb1zbeta);	

		/**************************************************/
		/* Group-specific trend models with no intercepts */
		/**************************************************/

		pwts[7] = 0;									/* no contribution from marginal for group-specific intercept-only model */

		pwts[1] = 0;									/* group-specific cubic trend model */
		call transpose(q3mbetag, q3tmbetag);
		call mult(q3tmbetag, q3prbeta, q3exp);			/* exponent from normal pdf features &g times */
		pwts[1] = pwts[1] - 0.5*&g*q3exp[1,1];			/* log scale */
		call det(q3Dbetag, q3wts);						/* determinant of prior precision matrix features &g times */
		pwts[1] = pwts[1] + 0.5*&g*log(q3wts);			/* log scale */

		pwts[2] = 0;									/* group-specific quad trend model */
		call transpose(q2mbetag, q2tmbetag);
		call mult(q2tmbetag, q2prbeta, q2exp);			/* exponent from normal pdf features &g times */
		pwts[2] = pwts[2] - 0.5*&g*q2exp[1,1];			/* log scale */
		call det(q2Dbetag, q2wts);						/* determinant of prior precision matrix features &g times */
		pwts[2] = pwts[2] + 0.5*&g*log(q2wts);			/* log scale */

		pwts[3] = 0;									/* group-specific linear trend model */
		call transpose(q1mbetag, q1tmbetag);
		call mult(q1tmbetag, q1prbeta, q1exp);			/* exponent from normal pdf features &g times */
		pwts[3] = pwts[3] - 0.5*&g*q1exp[1,1];			/* log scale */
		call det(q1Dbetag, q1wts);						/* determinant of prior precision matrix features &g times */
		pwts[3] = pwts[3] + 0.5*&g*log(q1wts);			/* log scale */
	
		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;								
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i]= nuarr[k]+Sarr[(k-1)*&n+i]; 		/* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;							/* off-diagonal elements are those of AR matrix Vgamma */
			    do j = i+1 to &n;
				    Vg[i,j] = (rhoarr[k]**(rts[j]-rts[i]))*nuarr[k];
				    Vg[j,i] = Vg[i,j];
			    end; 
		  	end; 
			call inv(Vg, Wg);							/* Wg = Vg^{-1} */

			/************************************************************/
			/* group-specific cubic trend model with no intercept (q=3) */
			/************************************************************/
			call mult(q3Xt, Wg, q3XtW);					/* multiply Xt and Wg */
			call mult(q3XtW, q3X, q3XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q3Dbetag, q3XtWX, q3DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q3XtW, Zvec, q3zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q3prbeta, q3zbeta, q3pbeta);	/* sum of prior and WLS contributions */
			call transpose(q3pbeta, q3tpbeta);			/* transpose */
			call inv(q3DXtWX, q3CI);					/* inverse of posterior precision matrix */
			call mult(q3CI, q3pbeta, q3pbeta);			/* re-scale pbeta */
			call mult(q3tpbeta, q3pbeta, q3exp);		/* exponent from normal pdf */
			pwts[1] = pwts[1] + 0.5*q3exp[1,1];			/* log scale */
			call det(q3DXtWX, q3wts);
			pwts[1] = pwts[1] - 0.5*log(q3wts); 		/* determinant for normalizing constant */

			/***********************************************************/
			/* group-specific quad trend model with no intercept (q=2) */
			/***********************************************************/
			call mult(q2Xt, Wg, q2XtW);					/* multiply Xt and Wg */
			call mult(q2XtW, q2X, q2XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q2Dbetag, q2XtWX, q2DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q2XtW, Zvec, q2zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q2prbeta, q2zbeta, q2pbeta);	/* sum of prior and WLS contributions */
			call transpose(q2pbeta, q2tpbeta);			/* transpose */
			call inv(q2DXtWX, q2CI);					/* inverse of posterior precision matrix */
			call mult(q2CI, q2pbeta, q2pbeta);			/* re-scale pbeta */
			call mult(q2tpbeta, q2pbeta, q2exp);		/* exponent from normal pdf */
			pwts[2] = pwts[2] + 0.5*q2exp[1,1];			/* log scale */
			call det(q2DXtWX, q2wts);
			pwts[2] = pwts[2] - 0.5*log(q2wts); 		/* determinant for normalizing constant */

			/*************************************************************/
			/* group-specific linear trend model with no intercept (q=1) */
			/*************************************************************/
			call mult(q1Xt, Wg, q1XtW);					/* multiply Xt and Wg */
			call mult(q1XtW, q1X, q1XtWX);				/* calculate XtWX, the precision matrix from WLS */
			call addmatrix(q1Dbetag, q1XtWX, q1DXtWX);	/* posterior precision matrix for beta is qDbetag + XtWX */
			call mult(q1XtW, Zvec, q1zbeta);			/* contribution to posterior mean from WLS */
			call addmatrix(q1prbeta, q1zbeta, q1pbeta);	/* sum of prior and WLS contributions */
			call transpose(q1pbeta, q1tpbeta);			/* transpose */
			call inv(q1DXtWX, q1CI);					/* inverse of posterior precision matrix */
			call mult(q1CI, q1pbeta, q1pbeta);			/* re-scale pbeta */
			call mult(q1tpbeta, q1pbeta, q1exp);		/* exponent from normal pdf */
			pwts[3] = pwts[3] + 0.5*q1exp[1,1];			/* log scale */
			call det(q1DXtWX, q1wts);
			pwts[3] = pwts[3] - 0.5*log(q1wts); 		/* determinant for normalizing constant */

		end;

		/******************************************/
		/* Common trend models with no intercepts */
		/******************************************/

		pwts[4] = 0;									/* common cubic trend model */
		call transpose(b3mbetag, b3tmbetag);
		call mult(b3tmbetag, b3prbeta, b3exp);			/* exponent from normal pdf */
		pwts[4] = pwts[4] - 0.5*b3exp[1,1];				/* log scale */
		call det(b3Dbetag, q3wts);						/* determinant of prior precision matrix */
		pwts[4] = pwts[4] + 0.5*log(q3wts);				/* log scale */		
	
		pwts[5] = 0;									/* common quad trend model */
		call transpose(b2mbetag, b2tmbetag);
		call mult(b2tmbetag, b2prbeta, b2exp);			/* exponent from normal pdf */
		pwts[5] = pwts[5] - 0.5*b2exp[1,1];				/* log scale */
		call det(b2Dbetag, q2wts);						/* determinant of prior precision matrix */
		pwts[5] = pwts[5] + 0.5*log(q2wts);				/* log scale */		
	
		pwts[6] = 0;									/* common linear trend model */
		call transpose(b1mbetag, b1tmbetag);
		call mult(b1tmbetag, b1prbeta, b1exp);			/* exponent from normal pdf */
		pwts[6] = pwts[6] - 0.5*b1exp[1,1];				/* log scale */
		call det(b1Dbetag, q1wts);						/* determinant of prior precision matrix */
		pwts[6] = pwts[6] + 0.5*log(q1wts);				/* log scale */				

		do k = 1 to &g;									/* cycle through each group independently */

			abeta[1,1] = a[k];							/* group-specific abeta vector */
			call mult(aX, abeta, aXbeta);				/* vector of intercepts a */
		    do i = 1 to &n;						
			  Zvec[i,1]=Yarr[(k-1)*&n+i] - aXbeta[i,1]; /* populate nx1 data vector Zvec = Yvec - a */
			  Vg[i,i] = nuarr[k] + Sarr[(k-1)*&n + i];  /* Vg = sum of Vgamma and the diagonal matrix of sampling variances */
			end; 
			do i = 1 to &n-1;
			  do j = i+1 to &n;
				Vg[i,j] = (rhoarr[k]**(rts[j] - rts[i]))*nuarr[k];
				Vg[j,i] = Vg[i,j];
			  end; 
		  	end; 
			call inv(Vg, Wg);								/* Wg = Vg^{-1} */

			/***********************************************/
			/* common cubic trend model with no intercepts */
			/***********************************************/
		    call mult(b3Xt, Wg, b3XtW);						/* multiply bXt and Wg */
		    call mult(b3XtW, b3X, b3XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb3XtWX, b3XtWX, sumb3XtWX);	/* cumulative matrix sum */
		    call mult(b3XtW, Zvec, b3zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb3zbeta,b3zbeta,sumb3zbeta);	/* cumulative matrix sum */

			/**********************************************/
			/* common quad trend model with no intercepts */
			/**********************************************/
		    call mult(b2Xt, Wg, b2XtW);						/* multiply bXt and Wg */
		    call mult(b2XtW, b2X, b2XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb2XtWX, b2XtWX, sumb2XtWX);	/* cumulative matrix sum */
		    call mult(b2XtW, Zvec, b2zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb2zbeta,b2zbeta,sumb2zbeta);	/* cumulative matrix sum */

			/************************************************/
			/* common linear trend model with no intercepts */
			/************************************************/
		    call mult(b1Xt, Wg, b1XtW);						/* multiply bXt and Wg */
		    call mult(b1XtW, b1X, b1XtWX);					/* calculate bXtWX */
		    call addmatrix(sumb1XtWX, b1XtWX, sumb1XtWX);	/* cumulative matrix sum */
		    call mult(b1XtW, Zvec, b1zbeta);			 	/* contributions to posterior mean from WLS */
		    call addmatrix(sumb1zbeta,b1zbeta,sumb1zbeta);	/* cumulative matrix sum */

		end;												/* end cycle through groups */

		/***********************************************/
		/* common cubic trend model with no intercepts */
		/***********************************************/
		call addmatrix(b3Dbetag, sumb3XtWX, b3DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b3prbeta, sumb3zbeta, b3pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b3pbeta, b3tpbeta);					/* transpose */
		call inv(b3DXtWX, b3CI);							/* inverse of posterior precision matrix  */
		call mult(b3CI, b3pbeta, b3pbeta);					/* re-scale pbeta */
		call mult(b3tpbeta, b3pbeta, b3exp);				/* exponent from normal pdf */
		pwts[4] = pwts[4] + 0.5*b3exp[1,1];					/* log scale */
		call det(b3DXtWX, b3wts);
		pwts[4] = pwts[4] - 0.5*log(b3wts); 				/* determinant for normalizing constant */

		/**********************************************/
		/* common quad trend model with no intercepts */
		/**********************************************/
		call addmatrix(b2Dbetag, sumb2XtWX, b2DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b2prbeta, sumb2zbeta, b2pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b2pbeta, b2tpbeta);					/* transpose */
		call inv(b2DXtWX, b2CI);							/* inverse of posterior precision matrix  */
		call mult(b2CI, b2pbeta, b2pbeta);					/* re-scale pbeta */
		call mult(b2tpbeta, b2pbeta, b2exp);				/* exponent from normal pdf */
		pwts[5] = pwts[5] + 0.5*b2exp[1,1];					/* log scale */
		call det(b2DXtWX, b2wts);
		pwts[5] = pwts[5] - 0.5*log(b2wts); 				/* determinant for normalizing constant */

		/************************************************/
		/* common linear trend model with no intercepts */
		/************************************************/
		call addmatrix(b1Dbetag, sumb1XtWX, b1DXtWX); 		/* posterior precision matrix for common beta is bDbetag + sumbXtWX */
		call addmatrix(b1prbeta, sumb1zbeta, b1pbeta);		/* sum of prior and the cumulative WLS contributions */
		call transpose(b1pbeta, b1tpbeta);					/* transpose */
		call inv(b1DXtWX, b1CI);							/* inverse of posterior precision matrix  */
		call mult(b1CI, b1pbeta, b1pbeta);					/* re-scale pbeta */
		call mult(b1tpbeta, b1pbeta, b1exp);				/* exponent from normal pdf */
		pwts[6] = pwts[6] + 0.5*b1exp[1,1];					/* log scale */
		call det(b1DXtWX, b1wts);
		pwts[6] = pwts[6] - 0.5*log(b1wts); 				/* determinant for normalizing constant */

		/************************************/
		/* Posterior model weights and draw */
		/************************************/
		do m = 1 to (1+(&p-1)*2);							/* calculate using differences on log-scale (Bayes factors) for numerical stability */
			qwtsum = 0;
			do l = 1 to (1+(&p-1)*2);
				qwtsum = qwtsum + wts[l]*exp(pwts[l]-pwts[m]);
			end;
			qwts[m] = wts[m]/qwtsum;
		end;
		qwtsum = .;											/* check for any missing values */
		do m = 1 to (1+(&p-1)*2);
			if qwts[m] = . then qwtsum = 0;					/* set qwtsum to zero if any missings */
		end;
		if qwtsum = 0 then do;								/* replace missing values with zeroes and cumulate sum (if applicable) */
			do m = 1 to (1+(&p-1)*2);
				if qwts[m] = . then qwts[m] = 0;
				qwtsum = qwtsum + qwts[m];
			end;
		end;
		if qwtsum > 0 then do;								/* rescale to sum to 1 (if applicable) */
			do m = 1 to (1+(&p-1)*2);						
				qwts[m] = qwts[m]/qwtsum;
			end;
		end;
		if qwtsum = 0 then do;								/* in case all weights are zero, use prior weights (if applicable) */
			qwts[m] = wts[m];
		end;

		flg = rand('table', qwts[1], qwts[2], qwts[3], 
				   qwts[4], qwts[5], qwts[6], qwts[7]);		/* returned draw from conditional posterior for model flag */
		
		endsub;
	run;
	quit;

%end;

%mend;

data _null_;
run;
