/*  
 * This file illustrates the use of the enhanced Modified Kalman Filter (MKF) macro (eMKF, v12) using example from:
 *
 * Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical Guidance for Using the Modified Kalman Filter 
 * in Small Domain Estimation at the National Center for Health Statistics. National Center for Health Statistics. 
 * Vital Health Stat 2(209). 2024. DOI: https://dx.doi.org/10.15620/cdc:xxxxxx.
 *
 * Main methodological differences between the enhanced and earlier MKF macros are described in README.md
 */

/* Specify the directory path:
 * The macro is assumed to be saved in a user directory called: ..\eMKF\MKFmacro
 * The data is assumed to be in user directory: ..\eMKF\MKFdata
 */
%let user_path = \\..; /* Replace \\.. with the applicable path */

options nonotes;
OPTIONS FORMCHAR="|----|+|---+=|-/\<>*";

/* To output the log and results to a .log and .html file */

ods listing close;
proc printto log="&user_path.\logfile.log" new;
ODS HTML5 path="&user_path.\"
(url=none)
body="Output_emkf.html";
run;

/* Define the data library */
libname sdata "&user_path\eMKF\MKFdata";

/* Read data from CSV file: NHANES obesity data with years 1999-2000 through 2017-March 2020 */
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

/* Print 20-line preview of the dataset */ 
proc print data=NHANESobesity(obs=20);
run;

/* Compile the enhanced Modified Kalman Filter macro (eMKF, v12) */
%include "&user_path\eMKF\MKFmacro\emkf_macro_v12.sas";

/* Conduct Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends ";
 %mkf(data		 = NHANESobesity,  
	  group		 = Population,
      time		 = Year, 
	  by 		 = Age,
	  outcome	 = Obesity, 
	  se 		 = SE_obesity,
	  neff		 = NEFF_obesity,
	  Bayesmodel = bma_cubic,
	  comparedto = MIN,
	  out		 = bmac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Conduct Maximum likelihood-based model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Maximum likelihood-based model averaging with up to cubic trends ";
 %mkf(data		 = NHANESobesity,  
	  group		 = Population,
      time		 = Year, 
	  by 		 = Age,
	  outcome	 = Obesity, 
	  se 		 = SE_obesity,
	  Bayesmodel = ,
	  slopes 	 = indep_cubic indep_quad indep_linear 
				   common_cubic common_quad common_linear 
				   dropped,
	  out		 = mac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Run using group-specific first-order autoregressive coefficients with fully Bayesian cubic trend */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Group-specific AR(1) coefficients with fully Bayesian cubic trend ";
 %mkf(data		 = NHANESobesity,  
	  group		 = Population,
      time		 = Year, 
	  by 		 = Age,
	  outcome	 = Obesity, 
	  se 		 = SE_obesity,
	  neff		 = NEFF_obesity,
	  ARmodel 	 = indep_ar,
	  Bayesmodel = full_cubic,
	  out		 = bfc
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

