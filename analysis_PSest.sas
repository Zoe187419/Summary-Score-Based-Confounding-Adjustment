***************************************************************************************************
* NAME: analysis_PSest.sas
* PROGRAM PURPOSE: Confounder summary score estimation:
*                             -Propensity score
* CREATED BY: XLI   
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org 
*
* INPUT DATA: individual-level data sets with baseline characteristics from participating data partners
*
* OUTPUT: Data sets with estimated propensity scores;
*         Distribution of propensity scores; 
*
* NOTE: The program was built for a method study where the analytical core had access to individual-
*       level data sets from all participating data partners in this 3-site distributed network study.
*       In a typical distributed data network study, the analytical core will not have access to these 
*       individual-level data sets, and this step will be carried out at each site of 
*       data partner. Thus some parts of the program will not be applicable.  
*    
**************************************************************************************************/;
dm 'log;clear;output;clear;odsresults;select all;clear;';

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

data bsln_f;
  set final.bsln_f;
      if bslnbmicat <=2 then bslnbmicat1 = 1;
 else bslnbmicat1 = bslnbmicat-1;

      if pcombined_score_num <0 then pccs=1;
 else if pcombined_score_num =0 then pccs=2;
 else if pcombined_score_num >0 then pccs=3;

      if icombined_score_num <0 then iccs=1;
 else if icombined_score_num =0 then iccs=2;
 else if icombined_score_num >0 then iccs=3;

      if comorbidscore <0 then ccs=1;
 else if comorbidscore =0 then ccs=2;
 else if comorbidscore >0 then ccs=3;

      if       year <=2007 then yearcat=1;
 else if 2008<=year <=2010 then yearcat=2;
 else if 2010< year        then yearcat=3;

	  numIPcat=(numIP >=1); 	
      numIScat=(numIS >=1);
      numEDcat=(numED >=1); 
run;

*estimate PS by site;
%macro ps(site=, out=, dp=);
proc logistic data=bsln_f descending;
  where outc="LT05" and fu='1Y';
  where also site="&site" ;
  class agecat yearcat bslnbmicat1 iccs /param=ref; 
  effect s_numav=spline(numav/naturalcubic);
  effect s_numoa=spline(numoa/naturalcubic);
  effect s_numGeneric=spline(numGeneric/naturalcubic);
  effect s_numClass=spline(numClass/naturalcubic);
  effect s_numRX=spline(numRX/naturalcubic); 
  model AGB = agecat female yearcat bslnbmicat1 
              iccs 
			  covar1--covar3 covar6--covar8 covar10 covar11 covar13 covar17--covar21
			  s_numav s_numoa numIPcat numIScat numEDcat
              s_numGeneric s_numClass s_numRX;
  output out=&out predicted=ps;
  title "Estimated PS for &dp in LT05-1Y";
run;

*distribution of PS;
proc univariate data = &out;
   class agb;
   var ps;
   histogram ps;
run;

*Creating PS treatment groups for plotting;
data &out;
  set &out;
  if agb=1 then agb_ps=ps;  else agb_ps=.;
  if agb=0 then rygb_ps=ps; else rygb_ps=.;
run;
ods graphics on;
proc kde data = &out;
   univar rygb_ps agb_ps/plots=densityoverlay;
   title "Propensity score distributions by treatment group";
run;
ods graphics off;
proc sgplot data=&out;
   title "Propensity score distributions by treatment group";
   title2 "Site &site";
   histogram agb_ps  / binstart=0 binwidth=0.005 transparency=0.5;
   histogram rygb_ps / binstart=0 binwidth=0.005 transparency=0.5;
   *density agb_ps  / type=kernel legendlabel='Kernel' lineattrs=(pattern=solid);
   *density rygb_ps / type=kernel legendlabel='Kernel' lineattrs=(pattern=dash);
   xaxis display=(nolabel);
   yaxis grid;
   keylegend / location=inside position=topright across=1;
run;
%mend;
%ps(site=01, out=final.ps01, dp=dp1);
%ps(site=02, out=final.ps02, dp=dp2);
%ps(site=03, out=final.ps03, dp=dp3);





