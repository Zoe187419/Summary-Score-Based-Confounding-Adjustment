***************************************************************************************************
* NAME: analysis_dataprep.sas
* PROGRAM PURPOSE: data prepreparation, data exploration, and baseline characteristics
* CREATED BY: XLI
* CONTACT INFO: xiaojuan_li@harvardpilgrim.org
*
* INPUT DATA: individual-level data sets from participating data partners
* 
* OUTPUT: Descriptive statistics of baseline characteristics for Table 1;
*         Data sets with baseline characteristics for confounder summary score estimation; 
*
* NOTE: The program was built for a method study where the analytical core had access to individual-
*       level data sets from all participating data partners in this 3-site distributed network study.
*       In a typical distributed data network study, the analytical core will not have access to these 
*       individual-level data sets, and this step will be carried out at each site of 
*       data partner. Some parts of the program will not be applicable.  
*       
**************************************************************************************************/;
dm 'log;clear;output;clear;odsresults;select all;clear;';

libname dp1 "myproject/data/DP01" access=readonly; 
libname dp2 "myproject/data/DP02" access=readonly; 
libname dp3 "myproject/data/DP03" access=readonly; 
libname final "myproject/data/Analysis";

*get overall information from each site;
data dp1;
  set dp1.r01_t2_cida;
  where level="000";
  site="01";
  keep group level npts episodes dispensings daysupp amtsupp eps_wEvents All_Events tte site;
data dp2;
  set dp2.r01_t2_cida;
  where level="000";
  site="02";
  keep group level npts episodes dispensings daysupp amtsupp eps_wEvents All_Events tte site;
data dp3;
  set dp3.r01_t2_cida;
  where level="000";
  site="03";
  keep group level npts episodes dispensings daysupp amtsupp eps_wEvents All_Events tte site;
run;

*stack the datasets together;
data all; 
  set dp1 dp2 dp3;
run;

*transpose the dataset from long to wide;
proc sort data = all;
   by group site;
proc transpose data = all out=all_site1 prefix=Nptsite;
    by group;
	id site;
	var Npts;
proc transpose data = all out=all_site2 prefix=Episite;
    by group;
	id site;
	var Episodes;
proc transpose data = all out=all_site3 prefix=EpswEsite;
    by group;
	id site;
	var Eps_wEvents;
proc transpose data = all out=all_site4 prefix=AllEsite;
    by group;
	id site;
	var All_Events;
proc transpose data = all out=all_site5 prefix=ttesite;
    by group;
	id site;
	var tte;
run;

*to create an overview of the number of pts at risk, number of events at each site; 
data all_site;
  merge all_site1 (drop=_name_) 
        all_site2 (drop=_name_)  
        all_site3 (drop=_name_) 
        all_site4 (drop=_name_) 
        all_site5 (drop=_name_);
  by group;
run;

proc export data = all_site 
	     outfile = "myproject/data/Analysis/Bariatric_overall_site.csv"
	        dbms = csv
	     replace;
run;

%macro cov(dt=, dp=, prefix= );
	data count;
	   set &dt;
	run;

	proc transpose data=count out=countl;
	   by group;
	run;
	data countl1y;
	   set countl;
       *abstract procedure;
 	   proc=scan(group, 2, '_'); 
	   *abstract outcome;
		   if find(group,'REOPER','i') ge 1 OR find(group, 'REHOSP', 'i') ge 1 or find(group, 'DEATH', 'i') ge 1
	   then outc=scan(group, 3, '_');
	   else outc=scan(group, 4, '_');
	   *abstract FU length;
	   fu=scan(group, -1, '_');
	   if fu="1Y" and outc="LT05" then output; 
	run;

	proc sort data = countl1y;
	   by _name_ _label_;
	proc transpose data = countl1y out=countw1y_&dp prefix=&prefix;
	   by _name_ _label_;
	   id proc;
	   var col1;
	run;
%mend;
%cov(dt=dp1.r01_baseline_1, dp=dp1, prefix=dp1_);
%cov(dt=dp2.r01_baseline_1, dp=dp2, prefix=dp2_);
%cov(dt=dp3.r01_baseline_1, dp=dp3, prefix=dp3_);

*merge the three sites together;
data countw1y;
   merge countw1y_dp1 countw1y_dp2 countw1y_dp3;
   by _name_;
run;
proc export data = countw1y 
	 outfile="myproject/data/Analysis/Bariatric_baseline_site.csv"
	 dbms=csv
	 replace;
run;

*To add additional baseline information;
%macro bsln(dt=, dp=, site=);
data bsln;
  set &dt;
  proc=scan(group, 2, '_');  *abstract procedure;

  if find(group,'REOPER','i') ge 1 OR find(group, 'REHOSP', 'i') ge 1 or find(group, 'DEATH', 'i') ge 1
  then outc=scan(group, 3, '_');  *abstract outcome;
  else outc=scan(group, 4, '_');

  fu=scan(group, -1, '_');        *abstract FU length;
  keep custompatid group year race age pcombined_score_num icombined_score_num diastolic systolic los losbfproc
       proc outc fu baselinebmi;
run; 

*To fill the missing bslnbmi with nonmissing values;
data bslnbmi;
  set &dt; 
  proc=scan(group, 2, '_');
       if find(group,'REOPER','i') ge 1 OR find(group, 'REHOSP', 'i') ge 1 or find(group, 'DEATH', 'i') ge 1
  then outc=scan(group, 3, '_');
  else outc=scan(group, 4, '_');
  fu=scan(group, -1, '_');
  keep custompatid group year baselinebmi proc outc fu;
run; 

*Because some patients have eligible episodes in 2 different years, thus 2 different baseline bmis;
proc sort data = bslnbmi;
   by year custompatid descending baselinebmi;
data bslnbmia;
   set bslnbmi;
   by year custompatid;
   bslnbmi=baselinebmi;
   if first.custompatid;
run; 
proc freq data=bslnbmia;
   tables fu;
run; 

proc sort data = bsln;
  by custompatid year;
proc sort data = bslnbmia;
  by custompatid year;
data bsln_a;
  merge bsln bslnbmia (keep=custompatid year bslnbmi);
  by custompatid year;
run;

data bsln_&dp;
  set bsln_a;
  
  site="&site";
       if 18<=age<45 then agecat=1;
  else if 45<=age<65 then agecat=2;
  else if 65<=age    then agecat=3;

       if     bslnbmi<35 then bslnbmicat=1;
  else if 35<=bslnbmi<40 then bslnbmicat=2;
  else if 40<=bslnbmi<50 then bslnbmicat=3;
  else if 50<=bslnbmi<60 then bslnbmicat=4;
  else if 60<=bslnbmi    then bslnbmicat=5;
run;
%mend;
%bsln(dt=dp1.r01_mstr, dp=dp1, site=01);
%bsln(dt=dp2.r01_mstr, dp=dp2, site=02);
%bsln(dt=dp3.r01_mstr, dp=dp3, site=03);

*Baseline characteristics by site;
%macro bsln1(dt=, dp=, dt1=, dt2=);
*distribution of continuous variables;
proc means data =&dt n nmiss mean std median q1 q3 min max nolabels;
   where outc="LT05" and fu='1Y';
   class group;
   var age bslnbmi los losbfproc pcombined_score_num icombined_score_num systolic diastolic ;
run;

*distribution of categorical variables;
ods rtf file="Table1-cat-&dp..rtf" style=journal bodytitle;
proc tabulate data =&dt;
  where outc="LT05" and fu='1Y';
  class group year race agecat bslnbmicat;
  table agecat race bslnbmicat year, group*(N*f=comma9.0 COLPCTN='%'*f=9.1)/rts=32;
run;
ods rtf close;

*create the final baseline dataset;
proc sort data=&dt;
   by custompatid group year;
proc sort data=&dt1 out=ps_&dp;
   by custompatid group year;
data &dt2;
   merge &dt ps_&dp (drop=age race);
   by custompatid group year;
run;

%mend;
%bsln1(dt=bsln_dp1, dp=dp1, dt1=dp1.r01_ads_psmatch_1, dt2=bsln_dp1f);
%bsln1(dt=bsln_dp2, dp=dp2, dt1=dp2.r01_ads_psmatch_1, dt2=bsln_dp2f);
%bsln1(dt=bsln_dp3, dp=dp3, dt1=dp3.r01_ads_psmatch_1, dt2=bsln_dp3f);


*Baseline characteristics in pooled dataset;
data final.bsln_f;*merge three datasets together;
  set bsln_dp1f 
      bsln_dp2f 
      bsln_dp3f;
  where proc ^= "SG";
  Female=(sex="F");
     AGB=(proc="AGB"); 
run;

proc freq data = final.bsln_f;
  where outc="LT05" and fu='1Y';
  tables site*proc*female/list;
run;

*distribution of continuous variables;
proc means data =final.bsln_f n nmiss mean std median q1 q3 min max nolabels;
   where outc="LT05" and fu='1Y';
   class group;
   var age bslnbmi los losbfproc
       pcombined_score_num icombined_score_num systolic diastolic 
       numav numed numip numis numoa numGeneric numClass numRX ;
run;

*distribution of categorical variables;
ods rtf file="Table1-cat-pooled.rtf" style=journal bodytitle;
proc tabulate data =final.bsln_f ;
  where outc="LT05" and fu='1Y';
  class group year race agecat bslnbmicat;
  table agecat race bslnbmicat year, group*(N*f=comma9.0 COLPCTN='%'*f=9.1)/rts=32;
run;
ods rtf close;



