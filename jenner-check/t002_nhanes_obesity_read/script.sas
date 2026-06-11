/* Adapted from Testing-and-implementation/Example-eMKF-NHANES-Data.sas
 * The DATA step that reads the NHANES adult obesity prevalence file
 * (1999-2000 through 2017-March 2020) into the eMKF stacked input layout,
 * followed by the 20-row PROC PRINT preview from the same example.
 * The upstream script reads nhanes9920obesity.csv with INFILE ... DELIMITER=','
 * DSD; the same comma-delimited DSD parsing is kept here, reading a sample of
 * the repo's real obesity rows inline via DATALINES so the bundle is fully
 * self-contained. The INFORMAT / FORMAT / INPUT statements are unchanged.
 */

/* Read NHANES obesity data into the eMKF stacked (long) input layout */
data NHANESobesity    ;
	%let _EFIERR_ = 0;
	infile datalines delimiter = ',' MISSOVER DSD lrecl=32767 ;

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
	datalines;
1999.5,1999-2000,"Black, non-Hispanic",18-24,0.26555875,0.041440195,113.572739,0.187084579,0.356647469,Yes
1999.5,1999-2000,"Black, non-Hispanic",25-44,0.42735,0.026901661,338.1542087,0.373985642,0.481998005,Yes
1999.5,1999-2000,"Black, non-Hispanic",45-64,0.426657608,0.024977453,392.1003781,0.367670718,0.487230712,Yes
1999.5,1999-2000,"Black, non-Hispanic",65+,0.417530684,0.034928817,199.3396492,0.348252418,0.489299383,Yes
1999.5,1999-2000,"White, non-Hispanic",18-24,0.177681538,0.03083708,153.6511608,0.120772182,0.247447617,Yes
1999.5,1999-2000,"White, non-Hispanic",25-44,0.252721585,0.020425115,452.6846502,0.213312703,0.2953959,Yes
1999.5,1999-2000,"White, non-Hispanic",45-64,0.354981814,0.034358996,193.9532809,0.287752199,0.42672031,Yes
1999.5,1999-2000,"White, non-Hispanic",65+,0.318314376,0.022186944,440.803899,0.275040962,0.364046908,Yes
1999.5,1999-2000,"Other race, non-Hispanic",18-24,0.170656575,0.065100333,33.39573528,0.063220131,0.3403857,No
1999.5,1999-2000,"Other race, non-Hispanic",25-44,0.446906139,0.071244255,48.69848538,0.304333389,0.596169443,Yes
1999.5,1999-2000,"Other race, non-Hispanic",45-64,0.272709901,0.07170283,38.57766967,0.142170274,0.439815335,Yes
1999.5,1999-2000,"Other race, non-Hispanic",65+,0.196128619,0.086458425,21.09176383,0.057788375,0.424855386,No
1999.5,1999-2000,Mexican American,18-24,0.266460281,0.026814662,271.8388087,0.21486899,0.32322295,Yes
1999.5,1999-2000,Mexican American,25-44,0.333944165,0.029551841,254.6920662,0.276297911,0.395508385,Yes
1999.5,1999-2000,Mexican American,45-64,0.373632533,0.03874706,155.8821451,0.297557902,0.45460983,Yes
1999.5,1999-2000,Mexican American,65+,0.361517743,0.048122385,99.67452338,0.267645863,0.463854979,Yes
1999.5,1999-2000,Other Hispanic,18-24,0.181052225,0.043943688,76.78335169,0.095348254,0.298361332,Yes
1999.5,1999-2000,Other Hispanic,25-44,0.275961151,0.019973061,500.8648699,0.19703902,0.366672712,Yes
1999.5,1999-2000,Other Hispanic,45-64,0.389136032,0.063094757,59.71169671,0.265594033,0.524063401,Yes
1999.5,1999-2000,Other Hispanic,65+,0.26662934,0.041270303,114.8037697,0.166639688,0.387777996,Yes
2001.5,2001-2002,"Black, non-Hispanic",18-24,0.286406368,0.031748249,202.7655022,0.225258401,0.353919012,Yes
2001.5,2001-2002,"Black, non-Hispanic",25-44,0.394721541,0.026554928,338.80966,0.342321186,0.44897918,Yes
2001.5,2001-2002,"Black, non-Hispanic",45-64,0.445159287,0.026593566,349.2453272,0.390221035,0.501109247,Yes
2001.5,2001-2002,"Black, non-Hispanic",65+,0.483392362,0.036901992,183.3838041,0.40912974,0.558200491,Yes
2001.5,2001-2002,"White, non-Hispanic",18-24,0.244576334,0.028922361,220.8704039,0.189403573,0.306742971,Yes
2001.5,2001-2002,"White, non-Hispanic",25-44,0.301951132,0.014582171,991.2383472,0.271496461,0.33376176,Yes
2001.5,2001-2002,"White, non-Hispanic",45-64,0.372202566,0.01585599,929.4203784,0.338625651,0.406723947,Yes
2001.5,2001-2002,"White, non-Hispanic",65+,0.363404047,0.026889997,319.9426411,0.310619903,0.418743502,Yes
2001.5,2001-2002,"Other race, non-Hispanic",18-24,0.10338461,0.034932991,75.96097474,0.035252025,0.222413415,No
2001.5,2001-2002,"Other race, non-Hispanic",25-44,0.214989576,0.036937499,123.696628,0.131894755,0.319559886,Yes
;
run;

/* Print 20-line preview of the dataset */
proc print data=NHANESobesity(obs=20);
run;
