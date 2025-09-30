/*  
 * This file illustrates the use of the enhanced Modified Kalman Filter (eMKF) macro using data from the National Survey of Family 
 * Growth (NSFG) in:
 *
 * Forrest SE, Rossen LM, Ahrens KA. Trends in Risk of Pregnancy Loss Among US Women by Metropolitan Status, 2000-2018. 
 * Paediatric and Perinatal Epidemiology. 2025. DOI:10.1111/ppe.70066.
 *
 * This file reads in the eMKF input dataset created in Example-eMKF-NSFG-Data.R
 * After running, the eMKF output from this file is read back into Example-eMKF-NSFG-Data.R to complete the analysis
 *
 * Pregnancies reported in the 2006-2010, 2011-2013, 2013-2015, 2015-2017, and 2017-2019 NSFG survey periods--excluding induced 
 * abortions, ongoing pregnancies, pregnancies conceived before 2000 and after 2018, and pregnancies with outcomes occurring at ages 
 * younger than 15 or older than 44 years--are used in this analysis. Pregnancy-level data is tabulated by:
 * - Maternal age group (15-19, 20-24, 25-29, 30-34, and 35-44 years)
 * - Metropolitan (metro) status (metropolitan and nonmetropolitan)
 * - Conception year interval (2000-2002, 2003-2004, 2005-2006, 2007-2008, 2009-2010, 2011-2012, 2013-2014, 2015-2016, and 2017-2018)
 *
 * Technical guidance for using the enhanced MKF macro is available from:
 *
 * Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical Guidance for Using the Modified Kalman Filter 
 * in Small Domain Estimation at the National Center for Health Statistics. National Center for Health Statistics. 
 * Vital Health Stat 2(209). 2024. DOI: 10.15620/cdc:157496.
 */

/* Specify the directory path:
 * The macro is assumed to be saved in a user directory called: ..\eMKF\MKFmacro
 * The data is assumed to be in user directory: ..\eMKF\MKFdata
 */
%let user_path = \\..; /* Replace \\.. with the applicable path */

/* Define the data library */
libname data "&user_path\eMKF\MKFdata";

/* Read in CSV file: NSFG pregnancy loss data from 2006 to 2019, tabulated maternal age group, metro status, and 
 * conception year interval, created in Example-eMKF-NSFG-Data.R 
 */
proc import datafile = "&user_path\eMKF\MKFdata\nsfg_2006-2019_pregnancy-loss.csv"
	out = data.nsfg
	dbms = csv 
	replace;
	getnames = yes;
run;

/* Print 20-line preview of the tabulated 2006-2019 NSFG pregnancy loss dataset */
proc print data = data.nsfg(obs=20);
run;

/* Prepare data for the eMKF macro */
data emkf_input; 
	set data.nsfg;
	keep agegrp_con metro_int year loss se DEff_loss neff;
run;

proc sort data = emkf_input;
	by agegrp_con metro_int year;
run;

/* Print 20-line preview of the eMKF input dataset */
proc print data = emkf_input(obs=20);
run;

/* Compile the eMKF macro */
%include "&user_path\eMKF\MKFmacro\emkf_macro.sas";

/* Conduct Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); /* Start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends";
 %mkf(data 		  = emkf_input,
	  group 	    = metro_int, 
	  by 	        = agegrp_con,
    time 		    = year, 
	  outcome 	  = loss, 
	  se 		      = se,  
	  neff 		    = neff,   
	  randomVars  = YES,
	  Bayesmodel  = bma_cubic, 
	  nbi 		    = 20000, 
	  nmc 		    = 50000,
	  thin 		    = 2,
	  mcmcplot 	  = NO,
    GRthreshold = 1.1,
	  modelprint  = NO,
	  seed 		    = 44,
    out 		    = bmac 
    );
data _null_;
 	dur = datetime() - &_timer_start; /* Stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Export eMKF output CSV files */
proc export data = bmac_pred
	outfile = "&user_path\eMKF\MKFdata\bmac_pred_nsfg_2006-2019_pregnancy-loss.csv" 
	dbms = csv
	replace;
run;

proc export data = work.bmac_bayeslogmod
	outfile = "&user_path\eMKF\MKFdata\bmac_bayeslogmod_nsfg_2006-2019_pregnancy-loss.csv" 
	dbms = csv
	replace;
run;

proc export data = work.bmac_bayeslogGR
	outfile = "&user_path\eMKF\MKFdata\bmac_bayeslogGR_nsfg_2006-2019_pregnancy-loss.csv" 
	dbms = csv
	replace;
run;