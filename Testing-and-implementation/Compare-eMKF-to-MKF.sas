/*  
 * This file tests the enhanced Modified Kalman Filter (MKF) macro (eMKF, v12) and compares its output against the original MKF macro
 * Additional test cases for eMKF are also provided
 *
 * Main methodological differences between the enhanced and earlier MKF macros are described in README.md
 */

/* Specify the directory path:
 *  The macro is assumed to be saved in a user directory called: ..\eMKF\MKFmacro
 *  The data is assumed to be in user directory: ..\eMKF\MKFdata
 */
%let user_path = \\..; /* Replace \\.. with the applicable path */

/* Specify the data library */
libname sdata "&user_path\eMKF\MKFdata";

/* Read in sample NHIS data with 8 timepoints included with the original release for the MKF macro */
data NHISdata; 
 set sdata.nhis;
run;

data NHISdata2;
 set NHISdata;
 /* calculate effective sample sizes for use with eMKF macro */
 stroke_neff = stroke*(1-stroke)/stroke_se**2;
 diabetes_neff = diabetes*(1-diabetes)/diabetes_se**2;
 keep race_ethnicity year stroke stroke_se stroke_neff diabetes diabetes_se diabetes_neff;
 rename race_ethnicity=race;
run;

/* Print 20-line preview of the NHISdata2 dataset */
proc print data=NHISdata2(obs=20);
run;

/* Simulated data provided with the release of the original MKF macro */
data withgender;
	length Race $20 Gender$20;
	label 	Gender =’Gender subset’
			Race =’Race group surveyed’
			Year=’Year of the survey’
			Disease = ’Prevalence of the Disease’
			SE = ’Prevalence Standard Error’
	;
	input Gender $ Race $ Year Disease SE @@;
	datalines;
	Male White 2000 0.2455 0.0028 Male White 2001 0.2454 0.0030
	Male White 2002 0.2461 0.0031 Male White 2003 0.2464 0.0031
	Female White 2000 0.1197 0.0231 Female White 2001 0.1239 0.0273
	Female White 2002 0.1109 0.0252 Female White 2003 0.1240 0.0259
	Male Black 2000 0.2527 0.0031 Male Black 2001 0.2650 0.0032
	Male Black 2002 0.2744 0.0033 Male Black 2003 0.2794 0.0033
	Female Black 2000 0.1309 0.0255 Female Black 2001 0.1914 0.0293
	Female Black 2002 0.1475 0.0270 Female Black 2003 0.1712 0.0309
	Male Chinese 2000 0.3138 0.0068 Male Chinese 2001 0.3221 0.0075
	Male Chinese 2002 0.3063 0.0074 Male Chinese 2003 0.3141 0.0072
	Female Chinese 2000 0.1666 0.0076 Female Chinese 2001 0.1804 0.0079
	Female Chinese 2002 0.1810 0.0085 Female Chinese 2003 0.1634 0.0082
	Male Indian 2000 0.3232 0.0073 Male Indian 2001 0.3215 0.0076
	Male Indian 2002 0.3439 0.0079 Male Indian 2003 0.3334 0.0077
	Female Indian 2000 0.1249 0.0086 Female Indian 2001 0.1227 0.0096
	Female Indian 2002 0.1161 0.0084 Female Indian 2003 0.1300 0.0086
	Male Hispanic 2000 0.2530 0.0373 Male Hispanic 2001 0.2263 0.0388
	Male Hispanic 2002 0.2100 0.0355 Male Hispanic 2003 0.2852 0.0383
	Female Hispanic 2000 0.1173 0.0289 Female Hispanic 2001 0.1083 0.0276
	Female Hispanic 2002 0.1055 0.0245 Female Hispanic 2003 0.1056 0.0226
	Male Other 2000 0.2696 0.0314 Male Other 2001 0.2686 0.0353
	Male Other 2002 0.2597 0.0354 Male Other 2003 0.3326 0.0364
	Female Other 2000 0.0649 0.0194 Female Other 2001 0.1304 0.0269
	Female Other 2002 0.1993 0.0178 Female Other 2003 0.2232 0.0184
	;
run;

data withgender2; /* to test behavior with only one group */
  set withgender;
  where Race = "White";
run;

/* Read NHANES data from CSV file: obesity data with years 1999-2000 through 2017-March 2020, used in:
 * Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical Guidance for Using the Modified Kalman Filter 
 * in Small Domain Estimation at the National Center for Health Statistics. National Center for Health Statistics. 
 * Vital Health Stat 2(209). 2024. DOI: https://dx.doi.org/10.15620/cdc:xxxxxx.
 */
data NHANESobesity    ; 
	%let _EFIERR_ = 0;
	infile "&user_path\eMKF\MKFdata\nhanes9920obesity.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;

	informat year best32. ;				/* Numeric timepoint variable to use in calculations */
	informat cycle $16. ;				/* Character timepoint variable (here: survey cycle year range) */
	informat population $32. ;			/* Population subgroups of interest (here: race and Hispanic origin) */
	informat age $6. ;					/* Age group variable (here: 18-24, 25-44, 45-64, 65+) */
	informat obesity best32. ;			/* Proportion of adults with obesity */
	informat se_obesity best32. ;		/* Design-based standard error */
	informat neff_obesity best32. ;		/* Design-based effective sample size */
	informat ll_obesity best32. ;		/* Lower limit for 95 percent Korn-Graubard confidence interval */
	informat ul_obesity best32. ;		/* Upper limit for 95 percent Korn-Graubard confidence interval */
	informat rel_obesity $4. ;			/* Yes if proportion meets NCHS data presentation standards, No if not */

	format year best12. ;
	format cycle $16. ;
	format population $32. ;
	format age $6. ;
	format obesity best16. ;
	format se_obesity best16. ;
	format neff_obesity best16. ;
	format ll_obesity best16. ;
	format ul_obesity best16. ;
	format rel_obesity $4. ;

	input
		year
		cycle $
		population  $
	 	age $
	 	obesity
	 	se_obesity
	 	neff_obesity
	 	ll_obesity
	 	ul_obesity
	 	rel_obesity $
	;
	if _ERROR_ then call symputx('_EFIERR_',1); 
run;

data NHANESobesity2; /* to test behavior with more than 15 groups (cross race by age) */
 set NHANESobesity;
 popbyage = catx(":", population, age);
 drop population age;
run;

proc sort data=NHANESobesity2;
  by popbyage year;
run;

/* Print 20-line preview of the obesity datasets */

proc print data=NHANESobesity(obs=20);
run;

proc print data=NHANESobesity2(obs=20);
run;

/**************************Comparison runs using extended MKF (eMKF) ********************************/

/* Compile extended MKF macro (eMKF v12) */
%include "&user_path\eMKF\MKFmacro\emkf_macro_v12.sas";

/***************/
/* One outcome */
/***************/

title "Enhanced MKF. One outcome: MLE-based independent linear trends by gender"; 	
/* results identical to those from original MKF */
%mkf(data= withgender, group= race, by= gender,
	 time= year, outcome= disease, se= se, 
	 slopes= indep_linear, 				 	/* note: independent in original MKF --> indep_linear in eMKF */
	 Bayesmodel= , 						 	/* As in original MKF, bayesmodel needs to be explicitly blank to override default*/
	 orpoly=NO,							 	/* NO to over-ride orthogonal polynomial default in eMKF */
	 checkSampleSize=NO,					/* NO to over-ride recommended minimum number of timepoints in eMKF */
	 modelprint= YES,					 	/* YES to over-ride eMKF default to not print intermediate model output */
	 out= emkf1by_cl
      );

title "Enhanced MKF. One outcome: MLE-based independent linear trends";				
/* results identical to those from original MKF */
%mkf(data= NHISdata2, group= race, 
  	 time= year, outcome= stroke, se= stroke_se, 
	 slopes= indep_linear, 
	 Bayesmodel= ,	
	 orpoly=NO,
	 checkSampleSize=NO,
 	 out= emkf1_il	 
      );

title "Enhanced MKF. One outcome: MLE-based common linear trend";					
/* results identical to those from original MKF */
%mkf(data= NHISdata2, group= race, 
	 time= year, outcome= stroke, se= stroke_se, 
	 slopes= common_linear , 				/* note: common in original MKF --> common_linear in eMKF */
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
	 out= emkf1_cl
      );

title "Enhanced MKF. One outcome: MLE model averaging estimates using independent and common linear trend(s) only"; 
/* results identical to those from original MKF */
%mkf(data= NHISdata2,  group= race, 
	 time= year, outcome= stroke, se= stroke_se,
	 slopes= indep_linear common_linear ,
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
	 out= emkf1_ic 
      );

title "Enhanced MKF. One outcome: MLE-based with no trend"; 	
/* estimates differ slightly from those in original MKF (different initial values for NLMIXED) */
/* original MSEs significantly larger due to keeping both intercept and linear term in the X matrix */
%mkf(data= NHISdata2,  group= race,
	 time= year, outcome= stroke, se= stroke_se, 
	 slopes= dropped ,
	 Bayesmodel=,  
	 orpoly=NO, 
	 checkSampleSize=NO,	
	 out= emkf1_d 
      );

title "Enhanced MKF. One outcome: MLE-based model averaging estimates using independent, common, and no linear trend(s)"; 
/* results differ from those in original MKF, as above */
%mkf(data= NHISdata2,  group= race,
	 time= year, outcome= stroke, se= stroke_se, 
	 slopes= indep_linear common_linear dropped,
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
	 out= emkf1_icd
      );

/****************/
/* Two outcomes */
/****************/

title "Enhanced MKF. Two outcomes: MLE-based independent linear trend"; 						
/* results identical to those from original MKF */
%mkf(data= NHISdata2, group= race,
	 time= year, outcome= stroke, se= stroke_se, outcome2= diabetes, se2= diabetes_se, 
	 slopes= indep_linear, 
	 Bayesmodel= , 
	 orpoly=NO,
	 checkSampleSize=NO,	
     out= emkf2_i
      );

title "Enhanced MKF. Two outcomes: MLE-based common linear trend"; 							
/* results identical to those from original MKF */
%mkf(data= NHISdata2,  group= race,
  	 time= year, outcome= stroke, se= stroke_se, outcome2= diabetes, se2= diabetes_se, 
	 slopes= common_linear, 
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
  	 out= emkf2_c
      );

title "Enhanced MKF. Two outcomes: MLE-based model averaging estimates using independent and common linear trend(s) only"; 
/* results identical to those from original MKF */
%mkf(data= NHISdata2,  group= race, 
	 time= year, outcome= stroke, se= stroke_se, outcome2= diabetes, se2= diabetes_se, 
	 slopes= indep_linear common_linear , 
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
	 out= emkf2_ic
      );

title "Enhanced MKF. Two outcomes: MLE-based with no trend";										
/* results differ from those in original MKF, as explained in single-outcome case */
%mkf(data= NHISdata2,  group= race, 
     time= year, outcome= stroke, se= stroke_se, outcome2= diabetes, se2= diabetes_se, 
	 slopes= dropped, 
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
  	 out= emkf2_d
      );

title "Enhanced MKF. Two outcomes: MLE-based model averaging estimates using independent, common, and no linear trend(s)"; 
/* results differ from those in original MKF, as explained in single-outcome case */
%mkf(data= NHISdata2,  group= race, 
  	 time= year, outcome= stroke, se= stroke_se, outcome2= diabetes, se2= diabetes_se, 
 	 slopes= indep_linear common_linear dropped , 
	 Bayesmodel=, 
	 orpoly=NO,
	 checkSampleSize=NO,	
 	 out= emkf2_icd
      );


/**********************************/
/* Additional test cases for eMKF */
/**********************************/

/* Special cases */

title "Enhanced MKF. One outcome: MLE-based independent linear trends with just 1 group per replication"; 	
%mkf(data= withgender2, group= race, by= gender,
     time= year, outcome= disease, se= se, 
	 slopes= indep_linear, 
	 Bayesmodel= , 
	 checkSampleSize = NO,	
 	 out= emkf1by_il01
          );

title "Enhanced MKF. One outcome: MLE-based independent linear trends with > 15 groups";
/* Number of groups was limited to 15 in original MKF due to the use IML function block to create block diagonal matrices */	
%mkf(data= NHANESobesity2, group= popbyage, 
     time= year, outcome= obesity, se= se_obesity, 
	 slopes= indep_linear, 
	 Bayesmodel= ,	
	 out= emkf2_il15
          );

/* Quadratic and cubic trend models */

title "Enhanced MKF. One outcome: MLE-based common quadratic trend with > 15 groups";				
%mkf(data= NHANESobesity2, group= popbyage, 
     time= year, outcome= obesity, se= se_obesity, 
	 slopes= common_quad, 
	 Bayesmodel= ,	
	 out= emkf2_cq15
          );

title "Enhanced MKF. One outcome: MLE model averaging estimates across all 5 models up to quadratic"; 
 %mkf(data= NHISdata2,  group= race, 
      time= year, outcome= stroke, se= stroke_se, 
	  slopes= indep_quad indep_linear common_quad common_linear dropped, 
	  Bayesmodel=,
	  out= emkf1_icd_q
          );

title "Enhanced MKF. One outcome: MLE model averaging estimates across all 7 models up to cubic"; 
 %mkf(data= NHISdata2,  group= race, 
      time= year, outcome= stroke, se= stroke_se, 
	  slopes= indep_cubic indep_quad indep_linear common_cubic common_quad common_linear dropped, 
	  Bayesmodel=,
	  out= emkf1_icd_c
          );

title "Enhanced MKF. Two outcomes: MLE model averaging estimates across all 5 models up to quadratic"; 
 %mkf(data= NHISdata2,  group= race, 
      time= year, outcome= stroke, se= stroke_se, outcome2=diabetes , se2=diabetes_se, 
	  slopes= indep_quad indep_linear common_quad common_linear dropped, 
	  Bayesmodel=,
	  out= emkf2_icd_q
          );

title "Enhanced MKF. Two outcomes: MLE model averaging estimates across all 7 models up to cubic"; 
 %mkf(data= NHISdata2,  group= race, 
      time= year, outcome= stroke, se= stroke_se, outcome2=diabetes , se2=diabetes_se, 
	  slopes= indep_cubic indep_quad indep_linear common_cubic common_quad common_linear dropped, 
	  Bayesmodel=,
	  out= emkf2_icd_c
          );


/*****************************/
/* Bayesian model test cases */
/*****************************/

/* Fully Bayesian linear trend model with disparities calculations */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Fully Bayesian linear trends with disparities calculations";
 %mkf(data= NHISdata2, group= race, 
	  comparedto=min, 						/* eMKF allows for comparisons to min or max and includes ratios */
      time=year, outcome=stroke, se = stroke_se,
	  Bayesmodel=full_linear,
	  randomVars = NO,						/* over-ride eMKF default of random sampling variances */ 
	  nbi=1000, nmc=5000, 					/* over-ride eMKF default values of nbi=10000 and nmc=50000 for testing purposes */
	  out=emkf1_bfl 
          );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Fully-Bayesian linear model with random sampling variances */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Fully Bayesian linear trend model with random sampling variances ";
 %mkf(data= NHISdata2,  group= race, 
      time= year, outcome= stroke, se = stroke_se ,
      neff = stroke_neff,         	  		/* (Effective) samples sizes required for random sampling variances (= default) */  				
	  Bayesmodel= full_linear, 
	  nbi=1000, nmc=5000, 					/* over-ride eMKF default values of nbi=10000 and nmc=50000 for testing purposes */
 	  out = emkf1_bfl_rv
     );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Common quadratic model with random sampling variances and group-specific AR parameters */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Common Bayesian quadratic trend with random sampling variances and group-specific AR parameters";
 %mkf(data= NHISdata2,  group= race, 
      time=year, outcome=stroke, se = stroke_se , 
      neff = stroke_neff,         	  		/* (Effective) samples sizes required for random sampling variances (= default) */  				
	  Bayesmodel=common_quad,
	  ARmodel = indep_AR, 					/* over-ride eMKF default of assuming common AR(1) parameters across groups */
	  nbi=1000, nmc=10000, 					/* over-ride eMKF default values of nbi=10000 and nmc=50000 for testing purposes */
	  out=emkf1_b1q_rv_ar
);
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Bayesian model averaging estimation with random sampling variances */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends with random sampling variances";
 %mkf(data= NHISdata2,  group= race,
      time=year, outcome=stroke, se = stroke_se ,  
      neff = stroke_neff,         	  		/* (Effective) samples sizes required for random sampling variances (= default) */  				
	  Bayesmodel=bma_cubic, 				/* bma_cubic is the eMKF default */
	  nbi=2500, nmc=25000, 					/* over-ride eMKF default values of nbi=10000 and nmc=50000 for testing purposes */
	  out=emkf1_bmac_rv 
          );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;


/*******************************************************/
/* Test cases with irregularly spaced data from NHANES */
/*******************************************************/

%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome (NHANES): Bayesian model averaging up to quadratic with random sampling variances";
 %mkf(data= NHANESobesity, group= Population, by=Age,
      time=Year,  outcome=obesity, se = se_obesity , 
	  neff=neff_obesity, 					/* (Effective) samples sizes required for random sampling variances (= default) */
	  Bayesmodel=bma_quad, 					/* over-ride eMKF default of bma_cubic */
	  nbi=2500, nmc=25000, 					/* increase nbi and nmc to 10000 and 50000 (default) or even higher to improve mixing */
	  out=emkf3_bmaq_rv
          );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;


/* Compare to MLE-based model averaging */
title "Enhanced MKF. One outcome (NHANES): MLE model averaging estimates across all 5 models up to quadratic"; 
 %mkf(data= NHANESobesity, group= Population, by=Age,
      time= year,  outcome= obesity, se= se_obesity, 
	  slopes= indep_quad common_quad indep_linear common_linear dropped, 
	  Bayesmodel=,
	  out= emkf3_icd_q
          );

