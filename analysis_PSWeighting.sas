***************************************************************************************************
* NAME: analysis_PSWeighting.sas
* PURPOSE: Weighted analysis 
*          Confounder summary score: Propensity Score
*          Confounding adjustment methods: inverse probability weighting (IPTW) and matching weighting (MW)
*          Data-sharing approaches: reference: pooled individual-level data
                                    comparison 1: risk-set data-sharing 
*                                   comparison 2: summary-table data-sharing (not established)
*                                   comparison 3: meta-analysis data-sharing
*          Outcomes: binary and survival
*
* CREATED BY: XLI
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org 
*
* INPUT DATA: data sets with estimated propensity scores
*
* OUTPUT: Results from privacy-protecting analytical methods where confounding is adjusted through 
*         PS weighting: inverse probability weighting and matching weighting; 
*
*
**************************************************************************************************/;

dm 'log;clear;output;clear;odsresults;select all;clear;';

%include "myproject/Program/metaanal.sas";                *macro for meta-analysis - available at 
                                                           https://www.hsph.harvard.edu/donna-spiegelman/software/metaanal/;

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

*1. estimate PS within each data site - see program analysis_PSest.sas;
*Merge back so every outcome group has PS;
%macro psdt(in=, out=);
data &in.est;
  set &in;
  keep custompatid site proc agb ps;
run;
proc sql;
  create table &out as
  select a.*, b.ps
  from bsln_f as a, &in.est as b
  where a.custompatid=b.custompatid and a.site=b.site and a.proc=b.proc and a.agb=b.agb;
quit;
%mend;
%psdt(in=final.ps01, out=ps01all);
%psdt(in=final.ps02, out=ps02all);
%psdt(in=final.ps03, out=ps03all);

*2. do PS weighting;
%macro PS_siptw(in=, outc=, fu=, outfile=, dp=, site=);
*prepare data for analysis;
*restrict to an outcome;
data &in._grp1;
  set &in;
  where outc=&outc and fu=&fu;
run;

proc univariate data = &in._grp1; 
   class agb;
   var ps;
   output out=ps_pct min=min max=max p1=p1 p99=p99 pctlpts=0.5 99.5 pctlpre=p;
run;

/* Labeling the percentiles at the lower extremes of the treated in macro variables 
   which can be called later.
   Defining the minimum, 0.5th percentiles, and 1st percentile of the treated */
data _Null_;
  set ps_pct;
  where agb=1;
  call symput ("treated_min", min);
  call symput ("treated_05", p0_5);
  call symput ("treated_1", p1);
run;

/* Labeling the percentiles at the upper extremes of the untreated in macro variables
   which can be called later. 
   Defining the maximum, 99th, and 99.5th percentile of the untreated. */
data _Null_;
  set ps_pct;
  where agb=0;
  call symput ("untreated_max", max);
  call symput ("untreated_99", p99);
  call symput ("untreated_995", p99_5);
run;

/*to trim non-overlapping regions of the PS distribution*/
*first calculate the marginal probability of treatment for the stablized IPTW;

proc means data = &in._grp1 (keep=ps);
 where &treated_min <=ps<= &untreated_max;
 var ps;
 output out=ps_mean mean=marg_prob;
run;

data _null_;
 set ps_mean;
 call symput("marg_prob",marg_prob);
run;

DATA &outfile;
  set &in._grp1;
  where &treated_min <=ps<= &untreated_max;
  *calculate IPTW;
       if agb=1 then iptw=1/ps;
  else if agb=0 then iptw=1/(1-ps);

  *calculating stabilized iptw;
  	   if agb=1 then siptw=&marg_prob/ps;
  else if agb=0 then siptw=(1-&marg_prob)/(1-ps);

  *calculate matching weight;
       if agb=1 then mw=(min(ps,1-ps))/ps;
  else if agb=0 then mw=(min(ps,1-ps))/(1-ps);

RUN;
  
proc means data = &outfile;
   var iptw siptw mw;
run;

*Getting effect estimate at each site;
*For survival outcome;
*IPTW survival analysis;
ods output HazardRatios=HR_iptw_&dp ParameterEstimates=HRparmest_iptw_&dp;
proc phreg data=&outfile;
   model followuptime_itt*event_itt(0)= agb;
   weight siptw;
   hazardratio agb;
   title 'siptw-weighted HR';
run;
data HRparmest_iptw_&dp;
  set HRparmest_iptw_&dp;
  site=&site;
run;

*MW survival analysis;
ods output HazardRatios=HR_mw_&dp ParameterEstimates=HRparmest_mw_&dp;
proc phreg data=&outfile;
   model followuptime_itt*event_itt(0)= agb;
   weight mw;
   hazardratio agb;
   title 'mw-weighted HR';
run;
data HRparmest_mw_&dp;
  set HRparmest_mw_&dp;
  site=&site;
run;

*For binary outcome;
*IPTW logistic regression;
*with proc logistic;
ods output ParameterEstimates=condparmest_iptw_&dp;
proc logistic data=&outfile;
    weight siptw;
	model event_itt (event="1")=agb;
	title 'siptw-weighted OR';
run;
data condparmest_iptw_&dp;
  set condparmest_iptw_&dp;
  site=&site;
  if variable="Intercept" then delete;
run;

*MW logistic regression;
*with proc logistic;
ods output ParameterEstimates=condparmest_mw_&dp;
proc logistic data=&outfile;
    weight mw;
	model event_itt (event="1")=agb;
    title 'mw-weighted OR';
run;
data condparmest_mw_&dp;
  set condparmest_mw_&dp;
  site=&site;
  if variable="Intercept" then delete;
run;

*IPTW;
*with proc genmod;
proc genmod data=&outfile descending;
    weight siptw;
	model event_itt (event="1")=agb/dist=binomial link=log;
	estimate 'AGB vs RYGB' agb 1 /exp;
	title 'siptw-weighted';
run;


*MW;
*with proc genmod;
proc genmod data=&outfile descending;
    weight mw;
	model event_itt (event="1")=agb/dist=binomial link=log;
	estimate 'AGB vs RYGB' agb 1 /exp;
    title 'mw-weighted';
run;

*Estimating subclass-specific and overall effect estimates;
* Binary outcome: Mantel-Haenszel stratified analysis ;
proc freq data= &outfile;
	table agb*event_Itt / nocol cmh;
    weight siptw;
	output out=mh_iptw_&dp(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
    title 'siptw-weighted';
run;

proc freq data= &outfile;
	table agb*event_Itt / nocol cmh;
    weight mw;
	output out=mh_mw_&dp(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
    title 'mw-weighted';
run;
title;
%mend;

*3. carry out analysis with each privacy-protecting analytical method;
%macro analysis();
*survival outcomes;
data hrparmest_iptw;
  set hrparmest_iptw_dp1  
      hrparmest_iptw_dp2
      hrparmest_iptw_dp3;
proc print data = hrparmest_iptw;
run;

data hrparmest_mw;
  set hrparmest_mw_dp1  
      hrparmest_mw_dp2
      hrparmest_mw_dp3;
proc print data = hrparmest_mw;
run;

data hr_iptw;
  set hr_iptw_dp1  
      hr_iptw_dp2
      hr_iptw_dp3;
proc print data = hr_iptw;
 title 'site-specific iptw HR';
run;

data hr_mw;
  set hr_mw_dp1  
      hr_mw_dp2
      hr_mw_dp3;
proc print data = hr_mw;
 title 'site-specific mw HR';
run;
title;

*Binary outcomes;
data condparmest_iptw;
  set condparmest_iptw_dp1  
      condparmest_iptw_dp2
      condparmest_iptw_dp3;
  exp=exp(estimate);
proc print data = condparmest_iptw;
 title 'site-specific iptw OR';
run;

data condparmest_mw;
  set condparmest_mw_dp1  
      condparmest_mw_dp2
      condparmest_mw_dp3;
  exp=exp(estimate);
proc print data = condparmest_mw;
 title 'site-specific mw OR';
run;

*MH;
data mh_iptw;
  set mh_iptw_dp1  
      mh_iptw_dp2
      mh_iptw_dp3;
proc print data = mh_iptw;
title 'site-specific iptw-MH OR';
run;

data mh_mw;
  set mh_mw_dp1  
      mh_mw_dp2
      mh_mw_dp3;
proc print data = mh_iptw;
title 'site-specific mw-MH OR';
run;

*********************************
*REFERENCE ANALYSIS 
*********************************;
*the result is used as the benchmark to evaluate the performance of the combinations 
*of data-sharing approaches and confonding adjustment methods;
*pooled data across sites;
data ps01_trimshare;
   set ps01_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps iptw siptw mw;
data ps02_trimshare;
   set ps02_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps iptw siptw mw;
data ps03_trimshare;
   set ps03_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps iptw siptw mw;
run;

*stack them together;
data ps_trimshare;
  set ps01_trimshare
      ps02_trimshare
	  ps03_trimshare;
run;

*survival outcome - IPW, site-stratified Cox PH;
proc phreg data=ps_trimshare;
   model followuptime_itt*event_itt(0)=agb;
   strata  site;
   weight siptw;
   hazardratio agb;
   title 'PS-siptw pooled HR';
run;

*survival outcome - MW, site-stratified Cox PH;
proc phreg data=ps_trimshare;
   model followuptime_itt*event_itt(0)=agb;
   strata  site;
   weight mw;
   hazardratio agb;
   title 'PS-mw pooled HR';
run;


/**using the exact conditional logistic;-can't use strata when using weight statement*/
/*proc logistic data=ps_trimshare exactonly;*/
/*    strata site;*/
/*	model event_itt (event="1")=agb;*/
/*    weight siptw;*/
/*	exact agb / estimate=both;*/
/*	title 'pooled OR';*/
/*run;*/

*binary outcome - IPW, site-stratified logistic;
proc freq data= ps_trimshare;
	table site*agb*event_Itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	weight siptw;
	title 'PS-siptw pooled MH OR';
run;

*binary outcome - MW, site-stratified logistic;
proc freq data= ps_trimshare;
	table site*agb*event_Itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	weight mw;
	title 'PS-mw pooled MH OR';
run;

title;
%mend;


*********************************
*COMPARISON 1-risk-set 
*********************************;
*first export ps_trimshare data to R;
*then use the R program to create the datasets first -https://github.com/kaz-yos/distributed/tree/master/R;
*then call the following macro for analysis;
%macro analyse_risksetbin(infile1=, infile2=);
DATA  &infile1;
	LENGTH outcome $ 4 method $ 5;
	INFILE  "myproject/data/Analysis/&infile1..txt" 
	DSD 
	LRECL= 59 ;
	INPUT site outcome method strata start_time stop_time A event W count;
DATA  &infile2;
	LENGTH outcome $ 4 method $ 5;
	INFILE  "myproject/data/Analysis/&infile2..txt" 
	DSD 
	LRECL= 59 ;
	INPUT site outcome method strata start_time stop_time A event W count;
RUN;

*survival outcomes - IPW, site-stratified Cox PH;
proc phreg data=&infile1 covs outest =res covout;
  freq count;
  weight w;
  model (start_time, stop_time)*event(0)=A /ties=efron;
  strata site;
  hazardratio a;
  title "PS-siptw risk-set HR";
run;

*survival outcomes - MW, site-stratified Cox PH;
proc phreg data=&infile2 covs outest =res covout;
  freq count;
  weight w;
  model (start_time, stop_time)*event(0)=A /ties=efron;
  strata site;
  hazardratio a;
  title "PS-mw risk-set HR";
run;

*binary outcomes;
*first need to summarize the weighted riskset data;
*need overall weights at individual site;
proc sql;
  create table w_all as
  select unique site, input(site,5.) as site1, AGB, sum(siptw) as sum_siptw, sum(mw) as sum_mw
  from ps_trimshare
  group by site, agb;
quit;

**FOR SIPTW analysis;
proc sql; 
*summarize weights from risk-set data for cases;
  create table w_case as 
  select unique site, A, sum(count) as case, sum(count*w) as sum_w, 1 as event
  from &infile1
  where event=1
  group by site, a;

*merge to calculate the weights for non-cases;
  create table wt as
  select a.site, a.A as AGB, a.sum_w as sum_siptw_1, b.sum_siptw as sum_iptw_all, b.sum_siptw - a.sum_w as sum_siptw_0
  from w_case as a left join w_all as b
  on a.site=b.site1 and a.A=b.AGB
  order by site, A;
quit; 
proc transpose data = wt out=wtlong1 prefix=sum_siptw;
   by site agb;
   var sum_siptw_1 sum_siptw_0;
run;
data wtlong;
   set wtlong1 ;
   event=scan(_name_, -1, '_');
   rename sum_siptw1=sum_siptw;
   drop _name_;
run;
*to get OR for binary outcomes;
proc freq data= wtlong;
	table site*agb*event / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	weight sum_siptw;
	title 'PS-siptw risk-set MH-OR';
run;

**FOR MW analysis;
proc sql; 
*summarize weights from risk-set data for cases;
  create table w_case as 
  select unique site, A, sum(count) as case, sum(count*w) as sum_w, 1 as event
  from &infile2 
  where event=1
  group by site, a;

*merge to calculate the weights for non-cases;
  create table wt as
  select a.site, a.A as AGB, a.sum_w as sum_mw_1, b.sum_mw as sum_mw_all, b.sum_mw - a.sum_w as sum_mw_0
  from w_case as a left join w_all as b
  on a.site=b.site1 and a.A=b.AGB
  order by site, A;
quit; 
proc transpose data = wt out=wtlong1 prefix=sum_mw;
   by site agb;
   var sum_mw_1 sum_mw_0;
run;
data wtlong;
   set wtlong1 ;
   event=scan(_name_, -1, '_');
   rename sum_mw1=sum_mw;
   drop _name_;
run;
*to get OR for binary outcomes;
proc freq data= wtlong;
	table site*agb*event / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	weight sum_mw;
	title 'PS-mw risk-set MH-OR';
run;
%mend;


*********************************
*COMPARISON2-summary-table 
*********************************;
* - method not established yet;


*********************************
*COMPARISON3-meta-analysis 
*********************************;
%macro meta(infile=, var=, outfile=, para=);
data &infile;
  set &infile;
  invar=1/(stderr**2);
run;

proc sql;
  create table sum as 
  select &var, sum(invar) as sum_invar
  from &infile
  group by &var;

  create table &infile.1 as
  select a.*, b.*
  from &infile as a join sum as b
  on a.&var=b.&var;
quit;

proc sql;
  create table &outfile as
  select unique &var, sum(estimate*invar/sum_invar) as ma_log&para,  1/sum_invar as ma_var, exp(calculated ma_log&para) as ma_&para,
  exp(calculated ma_log&para - 1.96*sqrt(calculated ma_var)) as ma_ll, exp(calculated ma_log&para + 1.96*sqrt(calculated ma_var)) as ma_ul
  from &infile.1
  group by &var;
quit;
  
proc print data = &outfile;
 title 'PS-weighted fixed-effect meta-analysis';
run;
title;

*both fixed-effect and random-effects meta-analysis;
%metaanal(data=&infile, beta=estimate, se_or_var=s, se=stderr, studylab=site, loglinear=t, coeff=t, printcoeff=t,
          explabel=AGB, outcomelabel=outcome, pooltype=both);
%mend;


*call macros to perform analysis on BMI loss <=5% outocme ;
%ps_siptw(in=ps01all, outc="LT05", fu='1Y', outfile=ps01_trim, dp=dp1, site="01");
%ps_siptw(in=ps02all, outc="LT05", fu='1Y', outfile=ps02_trim, dp=dp2, site="02");
%ps_siptw(in=ps03all, outc="LT05", fu='1Y', outfile=ps03_trim, dp=dp3, site="03");
%analysis();

*weighted analysis of risk-set data;
proc export data = ps_trimshare outfile ="myproject/data/Analysis/lt05.dta"
   replace;
run;
*use R program to create the datasets first -https://github.com/kaz-yos/distributed/tree/master/R;
*call the following macro to use datasets created in the R program for weighted analysis of the riskset data;
%analyse_risksetbin(infile1=lt05wrs, infile2=lt05wrs_mw);

*for meta-analysis;
%meta(infile=hrparmest_iptw,   var=parameter, outfile=ma_LT05, para=HR_iptw);
%meta(infile=hrparmest_mw,     var=parameter, outfile=ma_LT05, para=HR_mw);
%meta(infile=condparmest_iptw, var=variable,  outfile=ma_LT05, para=OR_iptw);
%meta(infile=condparmest_mw,   var=variable,  outfile=ma_LT05, para=OR_mw);

