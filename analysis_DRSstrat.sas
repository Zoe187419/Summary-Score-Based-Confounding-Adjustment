***************************************************************************************************
* NAME: analysis_DRSstrat.sas
* PURPOSE: Stratified analysis 
*          Confounder summary score: Disease Risk Score (DRS)
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
* INPUT DATA: data sets with estimated disease risk scores
*
* OUTPUT: Results from privacy-protecting analytical methods where confounding is adjusted through 
*         DRS stratification; 
*
*
**************************************************************************************************/;
dm 'log;clear;output;clear;odsresults;select all;clear;';

%include "myproject/Program/metaanal.sas"; *macro for meta-analysis- available at 
                                            https://www.hsph.harvard.edu/donna-spiegelman/software/metaanal/;

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

*1. estimate DRS within each data site - see program analysis_DRSest.sas;
*2. do DRS stratification;
%macro DRS_strat(in=, outc=, fu=, var=, outfile=, dp=, site=);
*prepare data for analysis;
*restrict to a outcome;
data &in._grp1;
  set &in;
  where outc=&outc and fu=&fu;
run;

*percentiles;
*check for nonoverlapping regions;
proc univariate data = &in._grp1;
   class agb;
   var &var;
   output out=drs_pct min=min max=max p1=p1 p99=p99 pctlpts=0.5 99.5 pctlpre=p;
run;

/* Labeling the percentiles at the lower extremes of the treated in macro variables which can be called later.
   Defining the minimum, 0.5th percentiles, and 1st percentile of the treated */
data _Null_;
  set drs_pct;
  where agb=1;
  call symput ("treated_min", min);
  call symput ("treated_05", p0_5);
  call symput ("treated_1", p1);
run;

/* Labeling the percentiles at the upper extremes of the untreated in macro variables which can be called later. 
   Defining the maximum, 99th, and 99.5th percentile of the untreated. */
data _Null_;
  set drs_pct;
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
    ranks drs_rank;
	var &var;
	where &treated_min <= &var <= &untreated_max;
run;

data &outfile;
    set &outfile;
	quintile=drs_rank + 1;
run;

proc means data = &outfile;
   class quintile agb;
   var &var;
run;

*Getting effect estimate at each site;
*For survival outcomes;
*DRSstratified survival analysis;
ods output HazardRatios=HR_&dp ParameterEstimates=HRparmest_&dp;
proc phreg data=&outfile;
   model followuptime_itt*event_itt(0)=agb;
   strata quintile;
   hazardratio agb;
data HRparmest_&dp;
  set HRparmest_&dp;
  site=&site;
run;

*For binary outcomes;
*DRSstratified logistic regression;
*getting conditional estimate;
ods output ParameterEstimates=condparmest_&dp;
proc logistic data=&outfile;
    strata quintile;
	model event_itt (event="1")=agb;
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
%macro analysis(in1=, in2=, in3=, out=);
*survival outcomes;
data hrparmest;
  set hrparmest_dp1  
      hrparmest_dp2
      hrparmest_dp3;
proc print data = hrparmest;
run;
data hr;
  set hr_dp1  
      hr_dp2
      hr_dp3;
proc print data = hr;
 title 'site-specific DRS-strat HR';
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
  title 'site-specific DRS-strat conditinal OR';
run;

*exact;
data exactparmest;
  set exactparmest_dp1  
      exactparmest_dp2
      exactparmest_dp3;
  exp=exp(estimate);
proc print data = exactparmest;
  title 'site-specific DRS-strat exact OR';
run;

*exact;
data mh;
  set mh_dp1  
      mh_dp2
      mh_dp3;
proc print data = mh;
	title 'site-specific DRS-strat MH OR';
run;

*********************************
*REFERENCE ANALYSIS 
*********************************;
*the result is used as the benchmark to evaluate the performance of the combinations 
*of data-sharing approaches and confonding adjustment methods;
*pooled data across sites;
data &in1.share;
   set &in1;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur drs_rank  quintile;
data &in2.share;
   set &in2;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur drs_rank  quintile;
data &in3.share;
   set &in3;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur drs_rank  quintile;
run;

*stack them together;
data &out;
  set &in1.share
      &in2.share
	  &in3.share;
run;

*survival outcome - DRS- and site-stratified Cox PH;
proc phreg data=&out;
   model followuptime_itt*event_itt(0)=agb;
   strata quintile site;
   hazardratio agb;
   title 'pooled DRS-strat HR';
run;

*binary outcomes - DRS- and site-stratified logistic;
*using the exact conditional logistic;
proc logistic data=&out exactonly;
    strata quintile site;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
  	title 'pooled DRS-strat OR';
run;

*using the MH estimate;
proc freq data= &out;
	table site*quintile*agb*event_Itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
    title 'DRS-strat MH';
run;

*Privacy-protecting data-sharing approaches;

*********************************
*COMPARISON 1-risk-set 
*********************************;
data drs01_forriskset;
  set &in1.share;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=1;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs01_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs01_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs01_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs01_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=1;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs01_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);
data r01_risksetdata_1_1; set r01_risksetdata_1_1; matchid=1;
data r01_risksetdata_1_2; set r01_risksetdata_1_2; matchid=2;
data r01_risksetdata_1_3; set r01_risksetdata_1_3; matchid=3;
data r01_risksetdata_1_4; set r01_risksetdata_1_4; matchid=4;
data r01_risksetdata_1_5; set r01_risksetdata_1_5; matchid=5;
data drs01_riskset;
  set r01_risksetdata_1_1-r01_risksetdata_1_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='01';
run;

data drs02_forriskset;
  set &in2.share;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=2;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs02_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs02_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs02_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs02_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=2;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs02_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);

data r01_risksetdata_2_1; set r01_risksetdata_2_1; matchid=1;
data r01_risksetdata_2_2; set r01_risksetdata_2_2; matchid=2;
data r01_risksetdata_2_3; set r01_risksetdata_2_3; matchid=3;
data r01_risksetdata_2_4; set r01_risksetdata_2_4; matchid=4;
data r01_risksetdata_2_5; set r01_risksetdata_2_5; matchid=5;
data drs02_riskset;
  set r01_risksetdata_2_1-r01_risksetdata_2_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='02';
run;

data drs03_forriskset;
  set &in3.share;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=3;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs03_forriskset(where=(quintile=1)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=2;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs03_forriskset(where=(quintile=2)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=3;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs03_forriskset(where=(quintile=3)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=4;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs03_forriskset(where=(quintile=4)),var=subgroup,COVARNUM=0);
%let RUNID=r01;
%let COMP=3;
%let LOOK=5;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=drs03_forriskset(where=(quintile=5)),var=subgroup,COVARNUM=0);

data r01_risksetdata_3_1; set r01_risksetdata_3_1; matchid=1;
data r01_risksetdata_3_2; set r01_risksetdata_3_2; matchid=2;
data r01_risksetdata_3_3; set r01_risksetdata_3_3; matchid=3;
data r01_risksetdata_3_4; set r01_risksetdata_3_4; matchid=4;
data r01_risksetdata_3_5; set r01_risksetdata_3_5; matchid=5;
data drs03_riskset;
  set r01_risksetdata_3_1-r01_risksetdata_3_5;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='03';
run;

data drs_riskset;
   set drs01_riskset drs02_riskset drs03_riskset;
run;

*survival outcome - case-centered logistic regression;
ods output ParameterEstimates=OR;
proc logistic data =drs_riskset;
	model case_exposure (event='1')= /offset=lnodds;
	*strata site;
run;
data HR;
  set OR;
  where variable="Intercept";
  HR =exp(estimate);
  HR_ll=exp(estimate -1.96*stderr);
  HR_ul=exp(estimate +1.96*stderr);
proc print data = HR;
   title 'DRS-strat riskset - HR';
run;
title;

proc datasets library=work;
 delete r01_risksetdata_1_1-r01_risksetdata_1_5 r01_risksetdata_2_1-r01_risksetdata_2_5 r01_risksetdata_3_1-r01_risksetdata_3_5;
run;

*binary outcome - site-stratified logistic regression;
proc sql;
  create table summary01 as
  select unique site, quintile as drs_stratum, exposure as AGB, count(exposure) as count 
  from drs01_forriskset
  group by site,drs_stratum, exposure;

  create table summary02 as
  select unique site, quintile as drs_stratum, exposure as AGB, count(exposure) as count 
  from drs02_forriskset
  group by site,drs_stratum, exposure;

  create table summary03 as
  select unique site, quintile as drs_stratum, exposure as AGB, count(exposure) as count
  from drs03_forriskset
  group by site,drs_stratum, exposure;
quit;
data summary;
  set summary01-summary03;
run;

*keep the first exposure prob and lnodds of each PS stratum which contains the exposure prob;
data prob;
  set drs_riskset;
  by site matchid;
  if first.matchid then output;
  keep site matchid exposureprobability lnodds;
run;

proc sql;
  create table drs_risksetbin as
  select unique site, matchid, case_exposure as agb, count(case_exposure) as case
  from drs_riskset
  group by site, matchid, case_exposure;
quit;

proc sql;
  create table drs_risksetbin1 as
  select a.*, b.exposureprobability, b.lnodds 
  from drs_risksetbin as a left join prob as b
  on a.site=b.site and a.matchid=b.matchid
  order by site, matchid, agb;
quit;

data drs_risksetbin1;
   set drs_risksetbin1;
   drs_stratum=1*matchid;
run;

proc sort data = drs_risksetbin1;
  by site drs_stratum agb;
proc sort data = summary;
  by site drs_stratum agb;
data drs_risksetbin2;
  merge drs_risksetbin1 summary;
  by site drs_stratum agb;
  case1=case; if case=. then case1=0;
run;
 
proc logistic data = drs_risksetbin2;
    strata site drs_stratum;
    model case1/count=agb;
	exact agb / estimate=both;
	title 'DRS-strat riskset -OR';
run;
title;


*********************************
*COMPARISON2-summary-table 
*********************************;
proc sql;
  create table drssur_stratST as
  select unique site, quintile, agb, sum(event_itt) as event, count(CustomPatid) as total, sum(followuptime_itt) as pt, log(calculated pt) as lnpt
  from &out
  group by site, quintile, agb;
quit;

*binary outcomes - site-stratified exact conditional logistic regression; 
proc logistic data=drssur_stratST exactonly;
    strata site quintile;
	model event/total=agb;
	exact agb / estimate=both;
	title 'DRS-strat summary table -OR';
run;
title;

*survival outcome ;
*if proc genmod doesnot work, use MH;
data drssur_stratSTs;
  set drssur_stratST;
  rygb=(agb=0); *hard to trick stdrate to calculate the agb(1vs0) so rygb(0vs1). Alternatively, could just take the inverse of 
                what calcuated using AGB;
proc stdrate data=drssur_stratSTs
             method=mh
			 stat=rate
			 effect=ratio;
 population group=rygb event=event total=pt ;
 strata site quintile/order=data ;
 title 'DRS-strat summary table -HR using MH';
run;

*survival outcome - site-stratified conditional poisson regression;
proc genmod data=drssur_stratST;
    strata site quintile;
	model event=agb /offset=lnpt dist=poisson link=log;
	exact agb / estimate=both;
    title 'DRS-strat summary table -HR';
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
  select unique &var, sum(estimate*invar/sum_invar) as ma_log&para, 1/sum_invar as ma_var, exp(calculated ma_log&para) as ma_&para, 
  exp(calculated ma_log&para - 1.96*sqrt(calculated ma_var)) as ma_ll, exp(calculated ma_log&para + 1.96*sqrt(calculated ma_var)) as ma_ul
  from &infile.1
  group by &var;
quit;
  
proc print data = &outfile;
 title 'DRS-strat fixed-effect meta-analysis';
run;
title;

*both fixed-effect and random-effects meta-analysis;
%metaanal(data=&infile, beta=estimate, se_or_var=s, se=stderr, studylab=site, loglinear=t, coeff=t, printcoeff=t,
          explabel=AGB, outcomelabel=outcome, pooltype=both);
%mend;

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;

*call macros to perform analysis on BMI loss <=5% outocme;
*using DRS estimated from binary outcome model;
%drs_strat(in=drs01_lt05, outc="LT05", fu='1Y', var=drs_bin, outfile=drsbin01_lt05trim, dp=dp1, site="01");
%drs_strat(in=drs02_lt05, outc="LT05", fu='1Y', var=drs_bin, outfile=drsbin02_lt05trim, dp=dp2, site="02");
%drs_strat(in=drs03_lt05, outc="LT05", fu='1Y', var=drs_bin, outfile=drsbin03_lt05trim, dp=dp3, site="03");
%analysis(in1=drsbin01_lt05trim, in2=drsbin02_lt05trim, in3=drsbin03_lt05trim, out=drsbin_lt05trimshare);
%meta(infile=condparmest, var=variable,  outfile=ma_LT05, para=OR_drss);
%meta(infile=hrparmest,   var=parameter, outfile=ma_LT05, para=HR_drss);

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;
*using DRS estimated from survival outcome model;
%drs_strat(in=drs01_lt05, outc="LT05", fu='1Y', var=drs_sur, outfile=drssur01_lt05trim, dp=dp1, site="01");
%drs_strat(in=drs02_lt05, outc="LT05", fu='1Y', var=drs_sur, outfile=drssur02_lt05trim, dp=dp2, site="02");
%drs_strat(in=drs03_lt05, outc="LT05", fu='1Y', var=drs_sur, outfile=drssur03_lt05trim, dp=dp3, site="03");
%analysis(in1=drssur01_lt05trim, in2=drssur02_lt05trim, in3=drssur03_lt05trim, out=drssur_lt05trimshare);
%meta(infile=condparmest, var=variable,  outfile=ma_LT05, para=OR_drss);
%meta(infile=hrparmest,   var=parameter, outfile=ma_LT05, para=HR_drss);

