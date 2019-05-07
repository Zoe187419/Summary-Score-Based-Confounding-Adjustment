***************************************************************************************************
* NAME: analysis_PSstrat.sas
* PURPOSE: Stratified analysis 
*          Confounder summary score: Propensity Score
*          Confounding adjustment methods: Stratification
*          Data-sharing approaches: reference: pooled individual-level data
                                    comparison 1: risk-set data-sharing 
*                                   comparison 2: summary-table data-sharing
*                                   comparison 3: meta-analysis data-sharing
*          Outcomes: binary and survival
*
* CREATED BY: XLI
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org 
*
* INPUT DATA: data sets with estimated propensity scores
*
* OUTPUT: Results from privacy-protecting analytical methods where confounding is adjusted through 
*         PS stratification; 
*
*
**************************************************************************************************/;

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;

%include "myproject/Program/ms_createriskset_test1.sas"; *macro to create risk-set data;
%include "myproject/Program/metaanal.sas";               *macro for meta-analysis- available at 
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

*2. do PS stratification;
%macro PS_strat(in=, outc=, fu=, outfile=, dp=, site=);
*prepare data for analysis;
*restrict to an outcome;
data &in._grp1;
  set &in;
  where outc=&outc and fu=&fu;
run;

*percentiles;
*check for nonoverlapping regions;
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

* can check the macro variable;
%put &untreated_max;

/*to trim non-overlapping regions of the PS distribution*/
proc rank data = &in._grp1 
        groups = 5         /*5 for quintiles, 10 for deciles*/
           out = &outfile;
    ranks ps_rank;
	var ps;
	where &treated_min <=ps<= &untreated_max;
run;

data &outfile;
    set &outfile;
	quintile=ps_rank + 1;
run;

proc means data = &outfile;
   class quintile agb;
   var ps;
run;

*Getting effect estimate at each site;
*For survival outcomes;
*PSstratified survival analysis;
ods output HazardRatios=HR_&dp ParameterEstimates=HRparmest_&dp;
proc phreg data=&outfile;
   model followuptime_itt*event_itt(0)=agb;
   strata quintile;
   hazardratio agb;
run;
data HRparmest_&dp;
  set HRparmest_&dp;
  site=&site;
run;

*For binary outcomes;
*PSstratified logistic regression;
*getting conditional estimate;
ods output ParameterEstimates=condparmest_&dp;
proc logistic data=&outfile;
    strata quintile;
	model event_itt (event="1")=agb;
run;
data condparmest_&dp;
  set condparmest_&dp;
  site=&site;
run;

*getting exact estimate;
ods output ExactParmEst=exactparmest_&dp;
proc logistic data=&outfile exactonly;
    strata quintile;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
run;
data exactparmest_&dp;
  set exactparmest_&dp;
  site=&site;
run;

*logistic regression using proc genmod;
proc sort data = &outfile;
   by quintile;
proc genmod data=&outfile descending;
    by quintile;
	model event_itt (event="1")=agb/dist=binomial link=log;
	estimate 'AGB vs RYGB' agb 1 /exp;
run;

*Estimating subclass-specific and overall effect estimates;
* Binary outcome: Mantel-Haenszel stratified analysis ;
proc freq data= &outfile;
	table quintile*agb*event_Itt / nocol cmh;
	output out=mh_&dp(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
run;
%mend;

*3. carry out analysis with each privacy-protecting analytical method;
%macro analysis();
*survival outcomes;
data hrparmest;
  set hrparmest_dp1  
      hrparmest_dp2
      hrparmest_dp3;
proc print data = hrparmest;
data hr;
  set hr_dp1  
      hr_dp2
      hr_dp3;
proc print data = hr;
  title 'site-specific PS-strat HR';
run;
title;

*Binary outcomes;
*conditional;
data condparmest;
  set condparmest_dp1  
      condparmest_dp2
      condparmest_dp3;
  exp=exp(estimate);
proc print data = condparmest;
  title 'site-specific PS-strat conditinal OR';
run;

*exact estimate;
data exactparmest;
  set exactparmest_dp1  
      exactparmest_dp2
      exactparmest_dp3;
  exp=exp(estimate);
proc print data = exactparmest;
  title 'site-specific PS-strat exact OR';
run;

*MH estimate;
data mh;
  set mh_dp1  
      mh_dp2
      mh_dp3;
proc print data = mh;
	title 'site-specific PS-strat MH OR';
run;
title;

*********************************
*REFERENCE ANALYSIS 
*********************************;
*the result is used as the benchmark to evaluate the performance of the combinations 
*of data-sharing approaches and confonding adjustment methods;
*pooled data across sites;
data ps01_trimshare;
   set ps01_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps ps_rank quintile;
data ps02_trimshare;
   set ps02_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps ps_rank quintile;
data ps03_trimshare;
   set ps03_trim;
   keep custompatid group site proc agb event_itt followuptime_itt ps ps_rank quintile;
*stack them together;
data ps_trimshare;
   set ps01_trimshare
       ps02_trimshare
	   ps03_trimshare;
run;

*survival outcomes - PS- and site-stratified Cox PH;
proc phreg data=ps_trimshare;
   model followuptime_itt*event_itt(0)=agb;
   strata quintile site;
   hazardratio agb;
   title 'pooled PS-strat HR';
run;

*binary outcomes - PS- and site-stratified logistic;
*using the exact conditional logistic;
proc logistic data=ps_trimshare exactonly;
    strata quintile site;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
  	title 'pooled PS-strat OR';
run;

*using the MH estimate;
proc freq data= ps_trimshare;
	table site*quintile*agb*event_itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
    title 'PS-strat MH OR';
run;

*Privacy-protecting data-sharing approaches;
*********************************
*COMPARISON 1-risk-set 
*********************************;
data ps01_forriskset;
  set ps01_trimshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=1;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);
data r01_risksetdata_1_1; set r01_risksetdata_1_1; matchid=1;
data r01_risksetdata_1_2; set r01_risksetdata_1_2; matchid=2;
data r01_risksetdata_1_3; set r01_risksetdata_1_3; matchid=3;
data r01_risksetdata_1_4; set r01_risksetdata_1_4; matchid=4;
data r01_risksetdata_1_5; set r01_risksetdata_1_5; matchid=5;
data ps01_riskset;
  set r01_risksetdata_1_1-r01_risksetdata_1_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='01';
run;

data ps02_forriskset;
  set ps02_trimshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=2;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);

data r01_risksetdata_2_1; set r01_risksetdata_2_1; matchid=1;
data r01_risksetdata_2_2; set r01_risksetdata_2_2; matchid=2;
data r01_risksetdata_2_3; set r01_risksetdata_2_3; matchid=3;
data r01_risksetdata_2_4; set r01_risksetdata_2_4; matchid=4;
data r01_risksetdata_2_5; set r01_risksetdata_2_5; matchid=5;
data ps02_riskset;
  set r01_risksetdata_2_1-r01_risksetdata_2_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='02';
run;

data ps03_forriskset;
  set ps03_trimshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=3;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);

data r01_risksetdata_3_1; set r01_risksetdata_3_1; matchid=1;
data r01_risksetdata_3_2; set r01_risksetdata_3_2; matchid=2;
data r01_risksetdata_3_3; set r01_risksetdata_3_3; matchid=3;
data r01_risksetdata_3_4; set r01_risksetdata_3_4; matchid=4;
data r01_risksetdata_3_5; set r01_risksetdata_3_5; matchid=5;
data ps03_riskset;
  set r01_risksetdata_3_1-r01_risksetdata_3_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='03';
run;

data ps_riskset;
   set ps01_riskset ps02_riskset ps03_riskset;
run;

*survival outcome - case-centered logistic regression;
ods output ParameterEstimates=OR;
proc logistic data =ps_riskset;
	model case_exposure (event='1')= /offset=lnodds;
	*strata site;
run;
data HR;
  set OR;
  where variable="Intercept";
  HR =exp(estimate);
  HR_ll=exp(estimate - 1.96*stderr);
  HR_ul=exp(estimate + 1.96*stderr);
run;
proc print data = HR;
   title 'PS-strat riskset - HR';
run;
title;

proc datasets library=work;
 delete r01_risksetdata_1_1-r01_risksetdata_1_5 r01_risksetdata_2_1-r01_risksetdata_2_5 r01_risksetdata_3_1-r01_risksetdata_3_5;
run;

*binary outcome - site-stratified logistic regression;
proc sql;
  create table summary01 as
  select unique site, quintile as ps_stratum, exposure as AGB, count(exposure) as count 
  from ps01_forriskset
  group by site,ps_stratum, exposure;

  create table summary02 as
  select unique site, quintile as ps_stratum, exposure as AGB, count(exposure) as count 
  from ps02_forriskset
  group by site,ps_stratum, exposure;

  create table summary03 as
  select unique site, quintile as ps_stratum, exposure as AGB, count(exposure) as count
  from ps03_forriskset
  group by site,ps_stratum, exposure;
quit;
data summary;
  set summary01-summary03;
run;

*keep the first exposure prob and lnodds of each PS stratum which contains the exposure prob;
data prob;
  set ps_riskset;
  by site matchid;
  if first.matchid then output;
  keep site matchid exposureprobability lnodds;
run;

proc sql;
  create table ps_risksetbin as
  select unique site, matchid, case_exposure as agb, count(case_exposure) as case
  from ps_riskset
  group by site, matchid, case_exposure;
quit;

proc sql;
  create table ps_risksetbin1 as
  select a.*, b.exposureprobability, b.lnodds 
  from ps_risksetbin as a left join prob as b
  on a.site=b.site and a.matchid=b.matchid
  order by site, matchid, agb;
quit;

data ps_risksetbin1;
   set ps_risksetbin1;
   ps_stratum=1*matchid;
run;

proc sort data = ps_risksetbin1;
  by site ps_stratum agb;
proc sort data = summary;
  by site ps_stratum agb;
data ps_risksetbin2;
  merge ps_risksetbin1 summary;
  by site ps_stratum agb;
  case1=case; if case=. then case1=0;
run;
 
/*proc logistic data =ps_risksetbin1;*/
/*	model agb (event='1')= /offset=lnodds;*/
/*	*strata site;*/
/*	title 'PS-strat riskset - OR';*/
/*run;*/
/*title;*/

proc logistic data = ps_risksetbin2;
    strata site ps_stratum;
    model case1/count=agb;
	exact agb / estimate=both;
	title 'PS-strat riskset -OR';
run;
title;


*********************************
*COMPARISON2-summary-table 
*********************************;
proc sql;
  create table ps_stratST as
  select unique site, quintile, agb, sum(event_itt) as event, count(CustomPatid) as total, 
         sum(followuptime_itt) as pt, log(calculated pt) as lnpt
  from ps_trimshare
  group by site, quintile, agb;
quit;

*binary outcomes - site-stratified exact conditional logistic regression; ; 
proc logistic data=ps_stratST exactonly;
    strata site quintile;
	model event/total=agb;
	exact agb / estimate=both;
	title 'PS-strat summary-table -OR';
run;
title;

*survival outcomes;
*if proc genmod doesnot work, use MH;
data ps_stratSTs;
  set ps_stratST;
  rygb=(agb=0); *hard to trick stdrate to calculate the agb(1vs0) so rygb(0vs1). Alternatively, could just take the inverse of 
                what is calcuated using AGB;
proc stdrate data=ps_stratSTs
             method=mh
			 stat=rate
			 effect=ratio;
 population group=rygb event=event total=pt ;
 strata site quintile/order=data ;
 title 'PS-strat summary-table -HR using MH';
run;

*survival outcome - site-stratified conditional poisson regression;
proc genmod data=ps_stratST;
    strata site quintile;
	model event=agb /offset=lnpt dist=poisson link=log;
	exact agb / estimate=both;
    title 'PS-strat summary-table -HR';
run;
title;
%mend;

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
  select unique &var, sum(estimate*invar/sum_invar) as ma_log&para, exp(calculated ma_log&para) as ma_&para, 1/sum_invar as ma_var,
  exp(calculated ma_log&para - 1.96*sqrt(calculated ma_var)) as ma_ll, exp(calculated ma_log&para + 1.96*sqrt(calculated ma_var)) as ma_ul
  from &infile.1
  group by &var;
quit;
  
proc print data = &outfile;
   title 'PS-strat fixed-effect meta-analysis';
run;
title;

*both fixed-effect and random-effects meta-analysis;
%metaanal(data=&infile, beta=estimate, se_or_var=s, se=stderr, studylab=site, loglinear=t, coeff=t, printcoeff=t,
          explabel=AGB, outcomelabel=outcome, pooltype=both);
%mend;

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;

*call macros to perform analysis on BMI loss <=5% outocme ;
%ps_strat(in=ps01all, outc="LT05", fu='1Y', outfile=ps01_trim, dp=dp1, site="01");
%ps_strat(in=ps02all, outc="LT05", fu='1Y', outfile=ps02_trim, dp=dp2, site="02");
%ps_strat(in=ps03all, outc="LT05", fu='1Y', outfile=ps03_trim, dp=dp3, site="03");
%analysis();
%meta(infile=condparmest, var=variable, outfile=ma_LT05, para=OR);
%meta(infile=hrparmest, var=parameter, outfile=ma_LT05, para=HR);
