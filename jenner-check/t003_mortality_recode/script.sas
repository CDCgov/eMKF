/* Adapted from Testing-and-implementation/Example-eMKF-Simulated-Mortality-Data.sas
 * The DATA step that reads state-level mortality data from external causes
 * (V01-Y89) into the eMKF layout, then the population-group recode/rename and
 * the catx() race-by-state grouping the example performs before modeling.
 * Upstream reads "Underlying Cause of Death, 1999-2020 States V01-Y89 by Race
 * and Ethnicity.csv" via INFILE ... DELIMITER=',' DSD; the same comma-delimited
 * DSD parsing and the INFORMAT/FORMAT/INPUT block are kept, with a sample of
 * the repo's real rows embedded inline via DATALINES. The recode IF/THEN,
 * RENAME, subsetting IF, and CATX logic are exactly as written upstream.
 */

/* State-level crude and age-adjusted mortality data from external causes (V01-Y89) */
data ExternalCauses    ;
	%let _EFIERR_ = 0;
	infile datalines delimiter = ',' MISSOVER DSD lrecl=32767 ;

	informat Cause $7. ;				/* Cause of death code  */
	informat Year_Code best32.;			/* Year (numeric) */
	informat State $20. ;				/* State name */
	informat State_Code best32. ;		/* State code (numeric) */
	informat Population_Group $56. ;	/* Race and Hispanic origin */
	informat Deaths best32. ;			/* Number of deaths (simulated/imputed if < 10) */
	informat Population_Size best32. ;	/* Population at risk */
	informat Crude_Rate best32. ;		/* Crude rate, per 100,000 */
	informat Crude_SE best32. ;			/* SE for crude rate, calculating assuming Poisson number of deaths */
	informat Age_Adjusted_Rate best32. ;/* Age-adjusted rate, per 100,000 */
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
	datalines;
"V01-Y89",1999,"Alabama",1,"All",3405,4430141,76.8598561535626,1.31716750874477,76.7187164697253,1.31580417419667
"V01-Y89",1999,"Alabama",1,"Hispanic or Latino",32,65183,49.092554807235,8.67841960249203,46.1610801526505,9.71004046228264
"V01-Y89",1999,"Alabama",1,"American Indian or Alaska Native, Not Hispanic or Latino",7,22132,31.6284113500813,11.9544158280526,29.5873539650331,11.284561300193
"V01-Y89",1999,"Alabama",1,"Asian or Pacific Islander, Not Hispanic or Latino",8,32953,24.2770005765788,8.58321586728428,19.9736770454933,7.15508700744688
"V01-Y89",1999,"Alabama",1,"Black or African American, Not Hispanic or Latino",956,1156898,82.634769875996,2.67259945712419,86.5880126831126,2.85221567706966
"V01-Y89",1999,"Alabama",1,"White, Not Hispanic or Latino",2402,3152975,76.1820185697635,1.55441140571592,74.0101569846148,1.51661927633679
"V01-Y89",1999,"Alaska",2,"All",473,624779,75.7067699138415,3.48100098929886,87.4985186415768,4.55444102458647
"V01-Y89",1999,"Alaska",2,"Hispanic or Latino",13,23997,54.17343834646,15.0250084404883,95.0420549481561,51.1714183825618
"V01-Y89",1999,"Alaska",2,"American Indian or Alaska Native, Not Hispanic or Latino",167,103680,161.072530864198,12.4641666505788,193.128733239025,16.4116635678658
"V01-Y89",1999,"Alaska",2,"Asian or Pacific Islander, Not Hispanic or Latino",12,29625,40.5063291139241,11.6931700089038,38.6962057006436,11.2573346017626
"V01-Y89",1999,"Alaska",2,"Black or African American, Not Hispanic or Latino",14,24886,56.2565297757775,15.0351900135576,48.1131296608722,12.90047105158
"V01-Y89",1999,"Alaska",2,"White, Not Hispanic or Latino",267,442591,60.3265769073479,3.69192655032935,66.2625511810003,4.48256416643767
"V01-Y89",1999,"Arizona",4,"All",3502,5023823,69.7078698831547,1.17794154973322,71.0737558482057,1.20491751557282
"V01-Y89",1999,"Arizona",4,"Hispanic or Latino",791,1235716,64.0114718915997,2.27598592401899,72.9122768883473,3.09600249112444
"V01-Y89",1999,"Arizona",4,"American Indian or Alaska Native, Not Hispanic or Latino",360,242308,148.571239909537,7.83039188182407,179.978494808379,10.3838060263242
"V01-Y89",1999,"Arizona",4,"Asian or Pacific Islander, Not Hispanic or Latino",35,100257,34.9102805789122,5.90091443300679,34.6762827175276,6.68583244979163
"V01-Y89",1999,"Arizona",4,"Black or African American, Not Hispanic or Latino",128,159387,80.307678794381,7.09826303210724,88.4899472102733,8.46032845930545
"V01-Y89",1999,"Arizona",4,"White, Not Hispanic or Latino",2188,3286155,66.5823736251029,1.4234283609296,62.0486484388202,1.34484188199984
"V01-Y89",1999,"Arkansas",5,"All",1929,2651860,72.7413966046473,1.65621045273866,71.8494793784833,1.64006089905598
"V01-Y89",1999,"Arkansas",5,"Hispanic or Latino",25,74627,33.4999397001085,6.69998794002171,27.6702284566086,6.57089258606535
"V01-Y89",1999,"Arkansas",5,"American Indian or Alaska Native, Not Hispanic or Latino",3,17293,17.3480599086336,10.0159070581673,20.1654044232787,11.6425016720977
"V01-Y89",1999,"Arkansas",5,"Asian or Pacific Islander, Not Hispanic or Latino",5,20980,23.8322211630124,10.658093315061,19.289516775464,8.90752342130367
"V01-Y89",1999,"Arkansas",5,"Black or African American, Not Hispanic or Latino",335,422065,79.3716607631526,4.33653707787263,85.7776933268488,4.82125812091401
"V01-Y89",1999,"Arkansas",5,"White, Not Hispanic or Latino",1561,2116895,73.7400768578508,1.86638886341915,71.3092324693971,1.81860249887442
"V01-Y89",1999,"California",6,"All",14681,33499204,43.8249219294882,0.361695688397154,45.4958689928676,0.37783021951262
"V01-Y89",1999,"California",6,"Hispanic or Latino",3702,10682606,34.6544653991732,0.569561986847638,39.797501220987,0.770841498119362
"V01-Y89",1999,"California",6,"American Indian or Alaska Native, Not Hispanic or Latino",93,221880,41.914548404543,4.34633620019513,43.125115433545,4.62996855382741
"V01-Y89",1999,"California",6,"Asian or Pacific Islander, Not Hispanic or Latino",947,3892494,24.3288750091843,0.790582210450632,26.5198177658207,0.898752726972822
"V01-Y89",1999,"California",6,"Black or African American, Not Hispanic or Latino",1427,2312110,61.7185168525719,1.63381731557098,64.9918074486694,1.76220868179805
"V01-Y89",1999,"California",6,"White, Not Hispanic or Latino",8512,16390114,51.9337449391749,0.562903349893517,48.4078523755103,0.530607506166995
;
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

/* Create race-by-state grouping to use in model */
data ExternalCausesGrouped;
  set ExternalCauses;
  PopByState = catx(":", Population, State);
run;

proc sort data=ExternalCausesGrouped;
  by Population Year;
run;

/* Print 20-line preview of the dataset */
proc print data=ExternalCausesGrouped(obs=20);
run;
