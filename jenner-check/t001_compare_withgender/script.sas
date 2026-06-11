/* Adapted from Testing-and-implementation/Compare-eMKF-to-MKF.sas
 * The simulated multi-group survey dataset distributed with the original MKF
 * macro, used here to demonstrate the eMKF input layout: one row per
 * group-by-time estimate, with a standard error column.
 * Smart-quote glyphs in the original LABEL statement are rendered as plain
 * double quotes; the labels themselves are unchanged.
 */

/* Simulated data provided with the release of the original MKF macro */
data withgender;
	length Race $20 Gender$20;
	label 	Gender ="Gender subset"
			Race ="Race group surveyed"
			Year="Year of the survey"
			Disease = "Prevalence of the Disease"
			SE = "Prevalence Standard Error"
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

proc print data=withgender;
run;

proc print data=withgender2;
run;
