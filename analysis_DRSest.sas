***************************************************************************************************
* NAME: analysis_DRSest.sas
* PROGRAM PURPOSE: Confounder summary score estimation:
*                             -Disease Risk score: separate models for binary and survival outcomes
* CREATED BY: XLI     
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org  
*
* INPUT DATA: individual-level data sets with baseline characteristics from participating data partners
*
* OUTPUT: Data sets with estimated disease risk scores;
*         Distribution of disease risk scores; 
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

*1.estimate DRS within each data site;
%macro drs(outc=, fu=, site=, out=, dp=);
*for binary outcome;
data &out;
  set bsln_f;
  where outc=&outc and fu=&fu;
  where also site="&site" ;
run;

proc logistic data=bsln_f descending;
  where outc=&outc and fu=&fu;
  where also site="&site";
  where also agb=0;              *for estimation, specify the untreated;
  class agecat yearcat bslnbmicat1 iccs /param=ref; 
  effect s_numav=spline(numav/naturalcubic);
  effect s_numoa=spline(numoa/naturalcubic);
  effect s_numGeneric=spline(numGeneric/naturalcubic);
  effect s_numClass=spline(numClass/naturalcubic);
  effect s_numRX=spline(numRX/naturalcubic); 
  model event_itt (event="1") =   agecat female yearcat bslnbmicat1 
					              iccs 
								  covar1--covar3 covar6--covar8 covar10 covar11 covar13 covar17--covar21
								  s_numav s_numoa numipcat numiscat numedcat
					              s_numGeneric s_numClass s_numRX;
  score data=&out out=&out(rename=(p_1=drs_bin));         *for prediction, specify the entire dataset;
  title "Estimated binary DRS for &dp";
run;

*for survival outcome;
data test;
  set &out;
  event_itt1=event_itt;
  if agb=1 then event_itt1=.; *to trick proc phreg to output estimate of linear predictors for agb users;
run;

proc phreg data=test;
  class agecat yearcat bslnbmicat1 iccs /param=ref; 
  effect s_numav=spline(numav/naturalcubic);
  effect s_numoa=spline(numoa/naturalcubic);
  effect s_numGeneric=spline(numGeneric/naturalcubic);
  effect s_numClass=spline(numClass/naturalcubic);
  effect s_numRX=spline(numRX/naturalcubic); 
  model followuptime_itt*event_itt1(0) =  agecat female yearcat bslnbmicat1 
							              iccs 
										  covar1--covar3 covar6--covar8 covar10 covar11 covar13 covar17--covar21
										  s_numav s_numoa numipcat numiscat numedcat
							              s_numGeneric s_numClass s_numRX;
  output out=&out xbeta=drs_sur;
  title "Estimated survival DRS for &dp";
run;
%mend;


*call macros to perform DRS estimation for BMI loss <=5% outocme ;
%drs(outc="LT05", fu='1Y', site=01, out=drs01_lt05,  dp=dp1); 
%drs(outc="LT05", fu='1Y', site=02, out=drs02_lt05,  dp=dp2);
%drs(outc="LT05", fu='1Y', site=03, out=drs03_lt05,  dp=dp3);



