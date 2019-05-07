***************************************************************************************************
* NAME: analysis_DRSmatching.sas
* PURPOSE: Matched analysis 
*          Confounder summary score: Disease Risk Score (DRS)
*          Confounding adjustment methods: Matching
*          Data-sharing approaches: reference: pooled individual-level data
                                    comparison 1: risk-set data-sharing 
*                                   comparison 2: summary-table data-sharing
*                                   comparison 3: meta-analysis data-sharing
*          DRS matching method: Nearest neighbor matching
*          Outcomes: binary and survival
*
* CREATED BY: XLI
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org 
*
* INPUT DATA: data sets with estimated disease risk scores
*
* OUTPUT: Results from privacy-protecting analytical methods where confounding is adjusted through 
*         DRS matching; 
*
*
**************************************************************************************************/;

dm 'log;clear;output;clear;odsresults;select all;clear;';

%include "myproject/Program/ms_nearestneighbormatch.sas";  *macro to do DRS matching;
%include "myproject/Program/ms_createriskset_test1.sas";   *macro to create risk-set data;
%include "myproject/Program/metaanal.sas";                 *macro for meta-analysis - available at 
                                                            https://www.hsph.harvard.edu/donna-spiegelman/software/metaanal/;

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

*1. estimate DRS within each data site - see program analysis_DRSest.sas;
*2. do DRS matching;
%macro drsmatch(in=, out=, var=);
data drstest;
  set &in;
  exposure=agb;
  patid=custompatid;
run;

%ms_NearestNeighborMatch(infile=drstest,
						 outfile=&out(rename=(trt_ps=trt_&var ctrl_ps=ctrl_&var)),  
						 matchvars=patid,
						 psvar=&var,
						 matchratio=1,
						 withreplacement=N,
						 caliper=0.2,
						 digits=0.0000000000000001);
%mend;
*Switching binary outcome;
%drsmatch(in=drs01_lt05, out=drsbin01_lt05match, var=drs_bin);
%drsmatch(in=drs02_lt05, out=drsbin02_lt05match, var=drs_bin);
%drsmatch(in=drs03_lt05, out=drsbin03_lt05match, var=drs_bin);
*Switching survival outcome;
%drsmatch(in=drs01_lt05, out=drssur01_lt05match, var=drs_sur);
%drsmatch(in=drs02_lt05, out=drssur02_lt05match, var=drs_sur);
%drsmatch(in=drs03_lt05, out=drssur03_lt05match, var=drs_sur);

title;
%macro drs_match(in1=, in2=, outc=, fu=, outfile=, dp=, site=);
*prepare data for analysis;
proc sql;
  create table matchtrt as 
  select a.*, b.set_num
  from &in1 as a inner join &in2 as b
  on a.custompatid=b.patid;

  create table matchctrl as 
  select a.*, b.set_num
  from &in1 as a inner join &in2 as b
  on a.custompatid=b.ctrl_patid;
quit;

data matchall;
  set matchtrt matchctrl;
run;

*restrict to an outcome;
data &outfile;
  set matchall;
  where outc=&outc and fu=&fu;
run;

*getting effect estimate at each site;
*For survival outcome;
*DRSmatched survival analysis;
ods output HazardRatios=HR_drsm_&dp ParameterEstimates=HRparmest_drsm_&dp;
proc phreg data=&outfile;  
   model followuptime_itt*event_itt(0)= agb;
   hazardratio agb;
data HRparmest_drsm_&dp;
  set HRparmest_drsm_&dp;
  site=&site;
run;

*For binary outcomes;
*DRSmatched logistic regression;
*with proc logistic;
ods output ParameterEstimates=condparmest_drsm_&dp;
proc logistic data=&outfile;
	model event_itt (event="1")=agb;
data condparmest_drsm_&dp;
  set condparmest_drsm_&dp;
  site=&site;
  if variable="Intercept" then delete;
run;

*getting exact estimate;
ods output ExactParmEst=exactparmest_drsm_&dp;
proc logistic data=&outfile exactonly;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
data exactparmest_drsm_&dp;
  set exactparmest_drsm_&dp;
  site=&site;
run;

*logistic regression using proc genmod;
proc genmod data=&outfile descending;
	model event_itt (event="1")=agb/dist=binomial link=log;
	estimate 'AGB vs RYGB' agb 1 /exp;
run;

* Mantel-Haenszel stratified analysis ;
proc sort data=&outfile;
   by descending AGB descending event_Itt;
run;
proc freq data= &outfile order=data;
	table agb*event_Itt / nocol cmh;
	output out=mh_drsm_&dp(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
run;
%mend;

*3. carry out analysis with each privacy-protecting analytical method;
%macro analysis(in1=, in2=, in3=, out=);
*survival outcomes;
data hrparmest_drsm;
  set hrparmest_drsm_dp1  
      hrparmest_drsm_dp2
      hrparmest_drsm_dp3;
proc print data = hrparmest_drsm;
run;

data hr_drsm;
  set hr_drsm_dp1  
      hr_drsm_dp2
      hr_drsm_dp3;
proc print data = hr_drsm;
  title 'site-specific DRS-matched HR';
run;
title;

*Binary outcomes;
*conditional;
data condparmest_drsm;
  set condparmest_drsm_dp1  
      condparmest_drsm_dp2
      condparmest_drsm_dp3;
  exp=exp(estimate);
proc print data = condparmest_drsm;
  title 'site-specific DRS-matched conditinal OR';
run;

*exact;
data exactparmest_drsm;
  set exactparmest_drsm_dp1  
      exactparmest_drsm_dp2
      exactparmest_drsm_dp3;
  exp=exp(estimate);
proc print data = exactparmest_drsm;
  title 'site-specific DRS-matched exact OR';
run;

*MH;
data mh_drsm;
  set mh_drsm_dp1  
      mh_drsm_dp2
      mh_drsm_dp3;
proc print data = mh_drsm;
  title 'site-specific DRS-matched MH OR';
run;

*********************************
*REFERENCE ANALYSIS 
*********************************;
*the result is used as the benchmark to evaluate the performance of the combinations 
*of data-sharing approaches and confonding adjustment methods;
*pooled data across sites;
data &in1.share;
   set &in1;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur set_num;
data &in2.share;
   set &in2;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur set_num;
data &in3.share;
   set &in3;
   keep custompatid group site proc agb event_itt followuptime_itt drs_bin drs_sur set_num;
run;

*stack them together;
data &out;
  set &in1.share
      &in2.share
	  &in3.share;
run;

*survival outcome -DRS-matched, site-stratified Cox PH;
ods output HazardRatios=HR_pool;
proc phreg data=&out;
   model followuptime_itt*event_itt(0)=agb;
   strata  site;
   hazardratio agb;
   title 'pooled DRS-matched HR';
run;
proc print data=hr_pool label;
   format hazardratio waldlower waldupper 9.5;
   var _numeric_;
   title 'pooled DRS-matched HR more precision';
run;

*binary outcomes -DRS-matched, site-stratified logistic;
*using the exact conditional logistic-can't use strata when using weight statement;
proc logistic data=&out exactonly;
    strata site;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
	title 'pooled DRS-matched OR';
run;

*using the MH estimate;
proc sort data=&out;
   by site descending AGB descending event_Itt;
proc freq data= &out;
	table site*agb*event_Itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	title 'DRS-matched MH OR';
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
%ms_CreateRiskSet(file=drs01_forriskset,var=subgroup,COVARNUM=0);

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
%ms_CreateRiskSet(file=drs02_forriskset,var=subgroup,COVARNUM=0);

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
%ms_CreateRiskSet(file=drs03_forriskset,var=subgroup,COVARNUM=0);

data drs01_riskset;
  set r01_risksetdata_1_1 ;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='01';
data drs02_riskset;
  set r01_risksetdata_2_1;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='02';
data drs03_riskset;
  set r01_risksetdata_3_1;
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
run;
proc print data = HR;
   title 'DRS-matched riskset - HR';
run;
title;

proc datasets library=work;
 delete r01_risksetdata_1_1 r01_risksetdata_2_1  r01_risksetdata_3_1;
run;

*binary outcome - DRS-matched, site-stratified logistic regression;
proc sql;
  create table summary01 as
  select unique site,  exposure as AGB, count(exposure) as count  
  from drs01_forriskset
  group by site, exposure;

  create table summary02 as
  select unique site,  exposure as AGB, count(exposure) as count 
  from drs02_forriskset
  group by site, exposure;

  create table summary03 as
  select unique site,  exposure as AGB, count(exposure) as count 
  from drs03_forriskset
  group by site, exposure;
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

proc sort data = drs_risksetbin1;
  by site agb;
proc sort data = summary;
  by site agb;
data drs_risksetbin2;
  merge drs_risksetbin1 summary;
  by site agb;
  case1=case; if case=. then case1=0;
run;
 
proc logistic data = drs_risksetbin2;
    strata site ;
    model case1/count=agb;
	exact agb / estimate=both;
	title 'DRS-matched riskset -OR';
run;
title;

*********************************
*COMPARISON2-summary-table 
*********************************;
proc sql;
  create table drssur_matchST as
  select unique site, agb, sum(event_itt) as event, count(CustomPatid) as total, sum(followuptime_itt) as pt, log(calculated pt) as lnpt
  from &out
  group by site, agb;
quit;

*binary outcomes - DRS-matched, site-stratified exact conditional logistic regression; 
proc logistic data=drssur_matchST exactonly;
    strata site;
	model event/total=agb;
	exact agb / estimate=both;
	title 'DRS-matched summary table -OR';
run;

*survival outcome - DRS-matched, site-stratified conditional poisson regression;
proc genmod data=drssur_matchST;
    strata site;
	model event=agb  /offset=lnpt dist=poisson link=log;
	exact agb  / estimate=both;
    title 'DRS-matched summary table -HR';
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
 title 'DRS-matched fixed-effect meta-analysis';
run;
title;

*both fixed-effect and random-effects meta-analysis;
%metaanal(data=&infile, beta=estimate, se_or_var=s, se=stderr, studylab=site, loglinear=t, coeff=t, printcoeff=t,
          explabel=AGB, outcomelabel=outcome, pooltype=both);
%mend;


*call macros to perform analysis on BMI loss <=5% outocme ;
*using DRS estimated from binary outcome model;
%drs_match(in1=drs01_lt05, in2=drsbin01_lt05match, outc="LT05", fu='1Y', outfile=drsbin01_lt05matched, dp=dp1, site='01');
%drs_match(in1=drs02_lt05, in2=drsbin02_lt05match, outc="LT05", fu='1Y', outfile=drsbin02_lt05matched, dp=dp2, site='02');
%drs_match(in1=drs03_lt05, in2=drsbin03_lt05match, outc="LT05", fu='1Y', outfile=drsbin03_lt05matched, dp=dp3, site='03');
%analysis(in1=drsbin01_lt05matched, 
          in2=drsbin02_lt05matched, 
          in3=drsbin03_lt05matched,          out=drsbin_lt05matchedshare);
%meta(infile=hrparmest_drsm,  var=parameter,outfile=ma_lt05, para=HR_drsm);
%meta(infile=condparmest_drsm, var=variable,outfile=ma_lt05, para=OR_drsm);

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;
*using DRS estimated from survival outcome model;
%drs_match(in1=drs01_lt05, in2=drssur01_lt05match, outc='LT05', outfile=drssur01_lt05matched, dp=dp1, site='01');
%drs_match(in1=drs02_lt05, in2=drssur02_lt05match, outc='LT05', outfile=drssur02_lt05matched, dp=dp2, site='02');
%drs_match(in1=drs03_lt05, in2=drssur03_lt05match, outc='LT05', outfile=drssur03_lt05matched, dp=dp3, site='03');
%analysis(in1=drssur01_lt05matched, 
          in2=drssur02_lt05matched, 
          in3=drssur03_lt05matched,            out=drssur_lt05matchedshare);
%meta(infile=hrparmest_drsm,   var=parameter, outfile=ma_lt05, para=HR_drsm);
%meta(infile=condparmest_drsm, var=variable,  outfile=ma_lt05, para=OR_drsm);
