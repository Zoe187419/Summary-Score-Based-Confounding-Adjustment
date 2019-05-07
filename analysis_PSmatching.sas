***************************************************************************************************
* NAME: analysis_PSmatching.sas
* PURPOSE: Matched analysis 
*          Confounder summary score: Propensity Score
*          Confounding adjustment methods: Matching
*          Data-sharing approaches: reference: pooled individual-level data
                                    comparison 1: risk-set data-sharing 
*                                   comparison 2: summary-table data-sharing
*                                   comparison 3: meta-analysis data-sharing
*          PS matching method: Nearest neighbor matching
*          Outcomes: binary and survival
*
* CREATED BY: XLI
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org 
*
* INPUT DATA: data sets with estimated propensity scores
*
* OUTPUT: Results from privacy-protecting analytical methods where confounding is adjusted through 
*         PS matching; 
*
*
**************************************************************************************************/;

dm 'log;clear;output;clear;odsresults;select all;clear;';

%include "myproject/Program/ms_nearestneighbormatch.sas"; *macro for PS matching;
%include "myproject/Program/ms_createriskset_test1.sas";  *macro to create risk-set data;
%include "myproject/Program/metaanal.sas";                *macro for meta-analysis - available at 
                                                           https://www.hsph.harvard.edu/donna-spiegelman/software/metaanal/;

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

*1. estimate PS within each data site - see program analysis_PSest.sas;
*2. do PS matching;
%macro psmatch(in=, out=);
data pstest;
  set &in;
  exposure=agb;
  patid=custompatid;
run;
%ms_NearestNeighborMatch(infile=pstest,
						 outfile=&out,  
						 matchvars=patid,
						 psvar=ps,
						 matchratio=1,
						 withreplacement=N,
						 caliper=0.2,
						 digits=0.0000000000000001);
%mend;
%psmatch(in=final.ps01, out=ps01_match);
%psmatch(in=final.ps02, out=ps02_match);
%psmatch(in=final.ps03, out=ps03_match);

*merge back so every outcome group has PS;
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

%macro ps_match(in1=, in2=, outc=, fu=, outfile=, dp=, site=);
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
*PSmatched survival analysis;
ods output HazardRatios=HR_psm_&dp ParameterEstimates=HRparmest_psm_&dp;
proc phreg data=&outfile;  
   model followuptime_itt*event_itt(0)= agb;
   hazardratio agb;
run;
data HRparmest_psm_&dp;
  set HRparmest_psm_&dp;
  site=&site;
run;

*For binary outcomes;
*PSmatched logistic regression;
*with proc logistic;
ods output ParameterEstimates=condparmest_psm_&dp;
proc logistic data=&outfile;
	model event_itt (event="1")=agb;
run;
data condparmest_psm_&dp;
  set condparmest_psm_&dp;
  site=&site;
  if variable="Intercept" then delete;
run;

*getting exact estimate;
ods output ExactParmEst=exactparmest_psm_&dp;
proc logistic data=&outfile exactonly;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
run;
data exactparmest_psm_&dp;
  set exactparmest_psm_&dp;
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
	output out=mh_psm_&dp(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
run;
%mend;

*3. carry out analysis with each privacy-protecting analytical method;
%macro analysis();
*survival outcomes;
data hrparmest_psm;
  set hrparmest_psm_dp1  
      hrparmest_psm_dp2
      hrparmest_psm_dp3;
proc print data = hrparmest_psm;
run;

data hr_psm;
  set hr_psm_dp1  
      hr_psm_dp2
      hr_psm_dp3;
proc print data = hr_psm;
  title 'site-specific PS-matched HR';
run;
title;

*Binary outcomes;
*conditional;
data condparmest_psm;
  set condparmest_psm_dp1  
      condparmest_psm_dp2
      condparmest_psm_dp3;
  exp=exp(estimate);
proc print data = condparmest_psm;
  title 'site-specific PS-matched conditinal OR';
run;

*exact estimate;
data exactparmest_psm;
  set exactparmest_psm_dp1  
      exactparmest_psm_dp2
      exactparmest_psm_dp3;
  exp=exp(estimate);
proc print data = exactparmest_psm;
  title 'site-specific PS-matched exact OR';
run;

*MH estimate;
data mh_psm;
  set mh_psm_dp1  
      mh_psm_dp2
      mh_psm_dp3;
proc print data = mh_psm;
  title 'site-specific PS-matched MH OR';
run;
title;

*********************************
*REFERENCE ANALYSIS 
*********************************;
*the result is used as the benchmark to evaluate the performance of the combinations 
*of data-sharing approaches and confonding adjustment methods;
*pooled data across sites;
data ps01_matchedshare;
   set ps01_matched;
   keep custompatid group site proc agb event_itt followuptime_itt ps set_num;
data ps02_matchedshare;
   set ps02_matched;
   keep custompatid group site proc agb event_itt followuptime_itt ps set_num;
data ps03_matchedshare;
   set ps03_matched;
   keep custompatid group site proc agb event_itt followuptime_itt ps set_num;
run;

*stack them together;
data ps_matchedshare;
  set ps01_matchedshare
      ps02_matchedshare
	  ps03_matchedshare;
run;

*survival outcomes -PS-matched, site-stratified Cox PH;
proc phreg data=ps_matchedshare;
   model followuptime_itt*event_itt(0)=agb;
   strata  site;
   hazardratio agb;
   title 'pooled PS-matched HR';
run;

*binary outcomes -PS-matched, site-stratified logistic;
*using the exact conditional logistic-can't use strata statement when using weight statement;
proc logistic data=ps_matchedshare exactonly;
    strata site;
	model event_itt (event="1")=agb;
	exact agb / estimate=both;
	title 'pooled PS-matched OR';
run;

*using the MH estimate;
proc sort data=ps_matchedshare;
   by site descending AGB descending event_Itt;
proc freq data= ps_matchedshare;
	table site*agb*event_Itt / nocol cmh;
	output out=mh_pool(keep=_mhor_ l_mhor u_mhor _mhrrc1_ l_mhrrc1 u_mhrrc1) cmh ;
	title 'PS-matched MH OR';
run;
title;

*Privacy-protecting data-sharing approaches;

*********************************
*COMPARISON 1-risk-set 
*********************************;
data ps01_forriskset;
  set ps01_matchedshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=1;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps01_forriskset,var=subgroup,COVARNUM=0);

data ps02_forriskset;
  set ps02_matchedshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=2;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps02_forriskset,var=subgroup,COVARNUM=0);

data ps03_forriskset;
  set ps03_matchedshare;
  patid=custompatid;
  exposure=agb;
  subgroup="Allpts";
run;
%let RUNID=r01;
%let COMP=3;
%let LOOK=1;
%let CurrentCat=Overall;
%ms_CreateRiskSet(file=ps03_forriskset,var=subgroup,COVARNUM=0);

data ps01_riskset;
  set r01_risksetdata_1_1 ;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='01';
data ps02_riskset;
  set r01_risksetdata_2_1;
  lnodds=log(exposureprobability/(1-exposureprobability));
  site='02';
data ps03_riskset;
  set r01_risksetdata_3_1;
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
data HR;
  set OR;
  where variable="Intercept";
  HR =exp(estimate);
  HR_ll=exp(estimate -1.96*stderr);
  HR_ul=exp(estimate +1.96*stderr);
proc print data = HR;
   title 'PS-matched riskset - HR';
run;
title;

proc datasets library=work;
 delete r01_risksetdata_1_1 r01_risksetdata_2_1  r01_risksetdata_3_1;
run;

*binary outcome - PS-matched, site-stratified logistic regression;
proc sql;
  create table summary01 as
  select unique site, exposure as AGB, count(exposure) as count 
  from ps01_forriskset
  group by site, exposure;

  create table summary02 as
  select unique site, exposure as AGB, count(exposure) as count 
  from ps02_forriskset
  group by site, exposure;

  create table summary03 as
  select unique site, exposure as AGB, count(exposure) as count
  from ps03_forriskset
  group by site, exposure;
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

proc sort data = ps_risksetbin1;
  by site agb;
proc sort data = summary;
  by site agb;
data ps_risksetbin2;
  merge ps_risksetbin1 summary;
  by site agb;
  case1=case; if case=. then case1=0;
run;

proc logistic data = ps_risksetbin2;
    strata site ;
    model case1/count=agb;
	exact agb / estimate=both;
	title 'PS-matched riskset -OR';
run;
title;


*********************************
*COMPARISON2-summary-table 
*********************************;
proc sql;
  create table ps_matchST as
  select unique site, agb, sum(event_itt) as event, count(CustomPatid) as total, sum(followuptime_itt) as pt, log(calculated pt) as lnpt
  from ps_matchedshare
  group by site, agb;
quit;

proc freq data = ps_matchedshare;
   tables site*agb*event_Itt/cmh;
run;

*binary outcome - PS-matched, site-stratified exact conditional logistic regression; 
proc logistic data=ps_matchST exactonly;
    strata site;
	model event/total=agb;
	exact agb / estimate=both;
	title 'PS-matched summary table -OR';
run;

*survival outcome - PS-matched, site-stratified conditional poisson regression;
proc genmod data=ps_matchST;
    strata site;
	model event=agb /offset=lnpt dist=poisson link=log;
	exact agb / estimate=both;
    title 'PS-matched summary table -HR';
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
 title 'PS-matched fixed-effect meta-analysis';
run;
title;

*both fixed-effect and random-effects meta-analysis;
%metaanal(data=&infile, beta=estimate, se_or_var=s, se=stderr, studylab=site, loglinear=t, coeff=t, printcoeff=t,
          explabel=AGB, outcomelabel=outcome, pooltype=both);
%mend;

dm 'log;clear;output;clear;odsresults;select all;clear;';
title;

*call macros to perform analysis on BMI loss <=5% outocme ;
%ps_match(in1=ps01all, in2=ps01_match, outc='LT05', fu='1Y', outfile=ps01_matched, dp=dp1, site='01');
%ps_match(in1=ps02all, in2=ps02_match, outc="LT05", fu='1Y', outfile=ps02_matched, dp=dp2, site="02");
%ps_match(in1=ps03all, in2=ps03_match, outc="LT05", fu='1Y', outfile=ps03_matched, dp=dp3, site="03");
%analysis();
%meta(infile=hrparmest_psm, var=parameter, outfile=ma_LT05, para=HR_psm);
%meta(infile=condparmest_psm, var=variable, outfile=ma_LT05, para=OR_psm);

