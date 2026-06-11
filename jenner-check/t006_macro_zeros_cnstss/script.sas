/* Exercises three list-building utility macros defined in
 * SAS-macro/emkf_macro_v24.sas:
 *   %zeros   - prints a number of comma-separated zeros
 *   %zeross  - prints a number of space-separated zeros
 *   %cnstss  - prints a constant the specified number of times
 * The eMKF macro uses these to build the parameter mean vectors passed to the
 * nonlinear mixed model (PROC NLMIXED) and to PROC MCMC. The three macro
 * definitions below are copied verbatim from the library; the caller drives
 * them and surfaces the generated vectors in a dataset.
 */

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
%mend zeros;

%macro zeross(t);
0
%do j = 1 %to %eval(&t-1);
  0
%end;
%mend zeross;

/* eMKF: Just like %zeross, this macro prints the specified constant a number of times */
%macro cnstss(v, t);
%sysevalf(&v)
%do j = 1 %to %eval(&t-1);
  %sysevalf(&v)
%end;
%mend cnstss;

/* Caller: build mean vectors of the kind passed to PROC NLMIXED / PROC MCMC.
 * The macros emit their elements on separate lines; translate the line breaks
 * to spaces and COMPBL so each generated vector lands in a single cell. */
data param_vectors;
  length kind $24 vector $80;
  kind = "zeros(5), comma sep";   vector = compbl(translate("%zeros(5)",      "  ", '0D0A'x)); output;
  kind = "zeross(5), space sep";  vector = compbl(translate("%zeross(5)",     "  ", '0D0A'x)); output;
  kind = "cnstss(0.5, 4)";        vector = compbl(translate("%cnstss(0.5, 4)", "  ", '0D0A'x)); output;
  kind = "cnstss(1, 3)";          vector = compbl(translate("%cnstss(1, 3)",   "  ", '0D0A'x)); output;
run;

proc print data=param_vectors noobs;
run;
