/*  
 * This file illustrates the enhanced Modified Kalman Filter (MKF) macro using simulated state-level mortality data.
 *
 * Annual state-level data were queried using CDC WONDER for the 21-year period 1999-2000. Data were tabulated by: 
 * - Age group (< 1 year, 1-4, 5-14, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75-84, 85 years and over)
 * - Race (American Indian or Alaska Native; Asian or Pacific Islander; Black or African American; and White)
 * - Hispanic origin (Hispanic or Latino and not Hispanic or Latino)
 *
 * States with numerator case counts <10 for selected combinations of year, age, and race and Hispanic origin,
 * are suppressed in CDC WONDER due to NCHS confidentiality protection rules. Missing cell case counts were
 * simulated/imputed, holding fixed the marginal count by state, year, age, and race and Hispanic origin.
 * A similar strategy was used with county-level data in Talih et al (2022), Population Health Metrics DOI: 10.1186/s12963-022-00288-1
 *
 * Technical guidance for using the enhanced MKF macro is available from:
 *
 * Talih M, Rossen LM, Patel P, Earp M, Parker JD. Technical Guidance for Using the Modified Kalman Filter 
 * in Small Domain Estimation at the National Center for Health Statistics. National Center for Health Statistics. 
 * Vital Health Stat 2(209). 2024. DOI: 10.15620/cdc:157496.
 *
 * Main methodological differences between the enhanced and earlier MKF macros are described in README.md.
 */

/* Specify the directory path:
 * The macro is assumed to be saved in a user directory called: ..\eMKF\MKFmacro
 * The data is assumed to be in user directory: ..\eMKF\MKFdata
 */
%let user_path = \\..; /* Replace \\.. with the applicable path */

/* Define the data library */
libname sdata "&user_path\eMKF\MKFdata";

/* State-level crude and age-adjusted mortality data from external causes (V01-Y89) */
data ExternalCauses    ; 
	%let _EFIERR_ = 0;
	infile "&user_path\eMKF\MKFdata\Underlying Cause of Death, 1999-2020 States V01-Y89 by Race and Ethnicity.csv" 
	delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;

	informat Cause $7. ;				/* Cause of death code  */
	informat Year_Code best32.;			/* Year (numeric) */
	informat State $20. ;				/* State name */
	informat State_Code best32. ;		/* State code (numeric) */
	informat Population_Group $56. ;	/* Race and Hispanic origin */
	informat Deaths best32. ;			/* Number of deaths (simulated/imputed if < 10) */
	informat Population_Size best32. ;	/* Population at risk */
	informat Crude_Rate best32. ;		/* Crude rate, per 100,000 */
	informat Crude_SE best32. ;			/* SE for crude rate, calculating assuming Poisson number of deaths */
	informat Age_Adjusted_Rate best32. ;/* Age-adjusted rate, per 100,000 - see https://wonder.cdc.gov/wonder/help/ucd-expanded.html */
	informat Age_Adjusted_SE best32. ;	/* SE for age-adjusted rate, calculating assuming Poisson number of deaths */

	format Cause $7. ;
	format Year_Code best16. ;
	format State $20. ;
	format State_Code best16. ;
	format Population_Group $56. ; 
	format Deaths best16. ;
	format Population_Size best16. ;
	format Crude_Rate best16. ;
	format Crude_SE best16. ;
	format Age_Adjusted_Rate best16. ;
	format Age_Adjusted_SE best16. ;

	input
		Cause $
		Year_Code
		State  $
	 	State_Code
	 	Population_Group $
	 	Deaths
	 	Population_Size
	 	Crude_Rate
	 	Crude_SE
	 	Age_Adjusted_Rate
		Age_Adjusted_SE
	;
	if _ERROR_ then call symputx('_EFIERR_',1); 
run;

/* State-level age-specific mortality data from external causes (V01-Y89) */
data ExternalCausesByAge    ; 
	%let _EFIERR_ = 0;
	infile "&user_path\eMKF\MKFdata\Underlying Cause of Death, 1999-2020 States V01-Y89 by Race, Ethnicity, and Age.csv" 
	delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;

	informat Cause $7. ;				/* Cause of death code  */
	informat Year_Code best32.;			/* Year (numeric) */
	informat State $20. ;				/* State name */
	informat State_Code best32. ;		/* State code (numeric) */
	informat Population_Group $56. ;	/* Race and Hispanic origin */
	informat Age_Group $11. ;			/* Age group */
	informat Deaths best32. ;			/* Number of deaths (simulated/imputed if < 10) */
	informat Population_Size best32. ;	/* Population at risk */
	informat Age_Specific_Rate best32. ;/* Age-specific rate, per 100,000 */
	informat Age_Specific_SE best32. ;	/* SE for age-specific rate, calculating assuming Poisson number of deaths */

	format Cause $7. ;
	format Year_Code best16. ;
	format State $20. ;
	format State_Code best16. ;
	format Population_Group $56. ; 
	format Age_Group $11. ;
 	format Deaths best16. ;
	format Population_Size best16. ;
	format Age_Specific_Rate best16. ;
	format Age_Specific_SE best16. ;

	input
		Cause $
		Year_Code
		State  $
	 	State_Code
	 	Population_Group $
		Age_Group $
	 	Deaths
	 	Population_Size
	 	Age_Specific_Rate
		Age_Specific_SE
	;
	if _ERROR_ then call symputx('_EFIERR_',1); 
run;

/* Remove "All" from the population groups and abbreviate */
data ExternalCauses;
  set ExternalCauses;
  if Population_Group ^= "All";
  if Population_Group = "Hispanic or Latino" then Population_Group = "Hispanic";
  if Population_Group = "American Indian or Alaska Native, Not Hispanic or Latino" then Population_Group = "AIAN, NH";
  if Population_Group = "Asian or Pacific Islander, Not Hispanic or Latino" then Population_Group = "API, NH";
  if Population_Group = "Black or African American, Not Hispanic or Latino" then Population_Group = "B, NH";
  if Population_Group = "White, Not Hispanic or Latino" then Population_Group = "W, NH";
  rename Population_Group = Population Year_Code = Year;
run;
data ExternalCausesByAge;
  set ExternalCausesByAge;
  if Population_Group ^= "All";
  if Population_Group = "Hispanic or Latino" then Population_Group = "Hispanic";
  if Population_Group = "American Indian or Alaska Native, Not Hispanic or Latino" then Population_Group = "AIAN, NH";
  if Population_Group = "Asian or Pacific Islander, Not Hispanic or Latino" then Population_Group = "API, NH";
  if Population_Group = "Black or African American, Not Hispanic or Latino" then Population_Group = "B, NH";
  if Population_Group = "White, Not Hispanic or Latino" then Population_Group = "W, NH";
  rename Population_Group = Population Year_Code = Year Age_Group = Age;
run;

/* For age-specific rates, create r-e by age grouping to use in model */
data ExternalCausesByAge;
  set ExternalCausesByAge;
  PopByAge = catx(":", Population, Age);
run;

/* Print 20-line previews of the datasets */ 
proc print data=ExternalCauses(obs=20);
run;
proc print data=ExternalCausesByAge(obs=20);
run;

/* Compile the enhanced Modified Kalman Filter macro (eMKF) */
%include "&user_path\eMKF\MKFmacro\emkf_macro.sas";

/*******************************************************************************/
/* Age-adjusted data stratified by state (borrowing strength across r-e groups */
/*******************************************************************************/

/* Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends ";
 %mkf(data		 = ExternalCauses,  
	  group		 = Population,
      time		 = Year, 
	  by 		 = State,
	  outcome	 = Age_Adjusted_Rate, 
	  se 		 = Age_Adjusted_SE,
	  randomVars = NO,		/* over-ride default due to underlying Poisson model for deaths whereby mean = variance */
	  Bayesmodel = bma_cubic,
	  out		 = bmac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Maximum likelihood-based model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Maximum likelihood-based model averaging with up to cubic trends ";
 %mkf(data		 = ExternalCauses,  
	  group		 = Population,
      time		 = Year, 
	  by 		 = State,
	  outcome	 = Age_Adjusted_Rate, 
	  se 		 = Age_Adjusted_SE,
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


/*******************************************************************************/
/* Age-adjusted data stratified by r-e group (borrowing strength across states */
/*******************************************************************************/

/* Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends ";
 %mkf(data		 = ExternalCauses,  
	  group		 = State,
      time		 = Year, 
	  by 		 = Population,
	  outcome	 = Age_Adjusted_Rate, 
	  se 		 = Age_Adjusted_SE,
	  randomVars = NO,		/* over-ride default due to underlying Poisson model for deaths whereby mean = variance */
	  Bayesmodel = bma_cubic,
	  out		 = bmac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Maximum likelihood-based model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Maximum likelihood-based model averaging with up to cubic trends ";
 %mkf(data		 = ExternalCauses,  
	  group		 = State,
      time		 = Year, 
	  by 		 = Population,
	  outcome	 = Age_Adjusted_Rate, 
	  se 		 = Age_Adjusted_SE,
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

/************************************************************************************/
/* Age-specific data stratified by state (borrowing strength across PopByAge groups */
/* Running the below code is not recommended. Each model may take 40 hours or more. */
/************************************************************************************/

/* 

/* Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends ";
 %mkf(data		 = ExternalCausesByAge,  
	  group		 = PopByAge,
      time		 = Year, 
	  by 		 = State,
	  outcome	 = Age_Specific_Rate, 
	  se 		 = Age_Specific_SE,
	  randomVars = NO,		/* over-ride default due to underlying Poisson model for deaths whereby mean = variance */
	  Bayesmodel = bma_cubic,
	  out		 = bmac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Maximum likelihood-based model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Maximum likelihood-based model averaging with up to cubic trends ";
 %mkf(data		 = ExternalCausesByAge,  
	  group		 = PopByAge,
      time		 = Year, 
	  by 		 = State,
	  outcome	 = Age_Specific_Rate, 
	  se 		 = Age_Specific_SE,
	  Bayesmodel = ,
	  slopes 	 = /*indep_cubic indep_quad indep_linear 
				   common_cubic common_quad common_linear */
				   dropped,
	  out		 = mac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;


/************************************************************************************/
/* Age-specific data stratified by PopByAge group (borrowing strength across states */
/************************************************************************************/

/* Bayesian model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Bayesian model averaging for up to cubic trends ";
 %mkf(data		 = ExternalCausesByAge,  
	  group		 = State,
      time		 = Year, 
	  by 		 = PopByAge,
	  outcome	 = Age_Specific_Rate, 
	  se 		 = Age_Specific_SE,
	  randomVars = NO,		/* over-ride default due to underlying Poisson model for deaths whereby mean = variance */
	  Bayesmodel = bma_cubic,
	  out		 = bmac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* Maximum likelihood-based model averaging estimation */
%let _timer_start = %sysfunc(datetime()); 	/* start timer */
title "Enhanced MKF. One outcome: Maximum likelihood-based model averaging with up to cubic trends ";
 %mkf(data		 = ExternalCausesByAge,  
	  group		 = State,
      time		 = Year, 
	  by 		 = PopByAge,
	  outcome	 = Age_Specific_Rate, 
	  se 		 = Age_Specific_SE,
	  Bayesmodel = ,
	  slopes 	 = indep_cubic indep_quad indep_linear 
				   common_cubic common_quad common_linear 
				   dropped,
	  out		 = mac
	  );
data _null_;
	dur = datetime() - &_timer_start;		/* stop timer */

 */
	put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;
