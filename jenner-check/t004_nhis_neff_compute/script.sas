/* Adapted from Testing-and-implementation/Compare-eMKF-to-MKF.sas
 * The NHISdata2 step that derives effective sample sizes (NEFFs) from a
 * prevalence estimate and its standard error, then keeps and renames the
 * columns the eMKF macro expects. Upstream this reads sdata.nhis (the NHIS
 * sample dataset shipped in the repo) into NHISdata via SET; here NHISdata is
 * built inline from a small set of race-by-year stroke and diabetes estimates
 * matching the same column shape, so the bundle runs self-contained. The NEFF
 * derivation, KEEP, and RENAME are exactly as written upstream.
 */

/* Stand-in for sdata.nhis: race-by-year stroke and diabetes prevalence + SEs */
data NHISdata;
  length race_ethnicity $20;
  input race_ethnicity $ year stroke stroke_se diabetes diabetes_se;
  datalines;
White 1999 0.0214 0.00100 0.0721 0.00190
White 2000 0.0221 0.00104 0.0743 0.00201
White 2001 0.0218 0.00102 0.0768 0.00205
White 2002 0.0226 0.00107 0.0795 0.00211
Black 1999 0.0302 0.00313 0.1099 0.00540
Black 2000 0.0318 0.00328 0.1142 0.00561
Black 2001 0.0331 0.00340 0.1188 0.00583
Black 2002 0.0345 0.00352 0.1231 0.00604
Hispanic 1999 0.0188 0.00271 0.0967 0.00612
Hispanic 2000 0.0195 0.00280 0.1014 0.00638
Hispanic 2001 0.0203 0.00289 0.1062 0.00665
Hispanic 2002 0.0211 0.00298 0.1109 0.00691
;
run;

data NHISdata2;
 set NHISdata;
 /* calculate effective sample sizes for use with eMKF macro */
 /* Please note that the calculation provided here is for illustration purposes only.
    Refinements include truncating effective sample size at the nominal sample size to avoid neff > n.
    See: Korn EL, Graubard BI. Analysis of Health Surveys. 1999. Wiley: New York.
  */
 stroke_neff = stroke*(1-stroke)/stroke_se**2;
 diabetes_neff = diabetes*(1-diabetes)/diabetes_se**2;
 keep race_ethnicity year stroke stroke_se stroke_neff diabetes diabetes_se diabetes_neff;
 rename race_ethnicity=race;
run;

/* Print 20-line preview of the NHISdata2 dataset */
proc print data=NHISdata2(obs=20);
run;
