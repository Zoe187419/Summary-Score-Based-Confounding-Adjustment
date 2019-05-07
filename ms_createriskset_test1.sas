****************************************************************************************************
*                                           PROGRAM OVERVIEW
****************************************************************************************************
*
* PROGRAM: ms_createriskset_test1.sas  
*
* Created (mm/dd/yyyy): 03/09/2015
* Last modified: 11/09/2015
* Version: 1.1
*
*--------------------------------------------------------------------------------------------------
* PURPOSE:
*   This macro will compute three aggregated datasets to run the sequential analysis on returned data from the DPs.  
*     These datasets will not contain any patient level data and their use will be equivalent to 
*     running the same analyses on patient level data.
*   
*  Program inputs:    
*    -Dataset containing the patient level data
*                                     
*  Program outputs:       
*	 -MSOC.&RUNID._RISKSETDATA_&COMP._&LOOK.
*    -MSOC.&RUNID._RISKDIFFDATA_&COMP._&LOOK.
*    -MSOC.&RUNID._SURVIVALDATA_&COMP._&LOOK.
*                                                                                                  
*  PARAMETERS:  
*  	-file		= Dataset containing the patient level data
*   -var		= The variable name which defines the subgroup
*   -COVARNUM	= The associated covarnum for the subgroup
*            
*  Programming Notes:                                                                                
*                                                                           
*
*
*--------------------------------------------------------------------------------------------------
* CONTACT INFO: 
*  xiaojuan_li@harvardpilgrim.org
*
*--------------------------------------------------------------------------------------------------
*  CHANGE LOG: 
*
*   Version   Date       Initials      Comment (reference external documentation when available)
*   -------   --------   --------   ---------------------------------------------------------------
*   1.1       11/09/15   EM          QCI-52 - adapting the matched KM cruves sur that the curve is
*                                    censored when the riskset becomes non-informative.
*
* 
***************************************************************************************************;
%MACRO ISDATA(dataset=);

%PUT =====> MACRO CALLED: ms_macros v1.0 => ISDATA;

	%GLOBAL NOBS;
	%let NOBS=0;
	%if %sysfunc(exist(&dataset.))=1 and %LENGTH(&dataset.) ne 0 %then %do;
		data _null_;
		dsid=open("&dataset.");
		call symputx("NOBS",attrn(dsid,"NLOBS"));
		run;
	%end;	
%PUT &NOBS.;

%put NOTE: ********END OF MACRO: ms_macros v1.0 => ISDATA ********;

%MEND ISDATA; 

%macro ms_CreateRiskSet(file=,var=,COVARNUM=);

  %put =====> MACRO CALLED: ms_CreateRiskSet v1.1;

  data RS_file0;
  set &file.;
  *where missing(&varmiss.)=0 and missing(PatId)=0;
  where missing(PatId)=0;
  Exp=0;UnExp=0;EVExp=0;EVUnExp=0;FUTimeExp=0;FUTimeUnExp=0;
  if exposure=1 then do;
    Exp=1;
    FUTimeExp=Followuptime_ITT;
    if Event_ITT then EVExp=1;
  end;
  if exposure=0 then do;
    UnExp=1;
    FUTimeUnExp=Followuptime_ITT;
    if Event_ITT then EVUnExp=1;
  end;

  if "&var." in ("AllPts","percentile") then Type="AllPts ";
  else Type="Matched";

  if "&var." ne "percentile" then Percentile="0";

  keep &var. Type Followuptime_ITT Exp UnExp EVExp EVUnExp FUTimeExp FUTimeUnExp Percentile;
  run;

  %ISDATA(dataset=RS_file0);
  %IF %EVAL(&NOBS.>0) %THEN %DO;

    /******************************/
    /* aggregate data for riskSet */
    /******************************/

    *Create one record per Followuptime (collapse);
    proc means data=RS_file0 nway noprint;
    var Exp UnExp EVExp EVUnExp;
    class &var. Type Followuptime_ITT;
    output out=RS_file(drop=_:) sum= FU_Exp FU_UnExp FU_EVExp FU_EVUnExp;
    run;

    *Total Number of Pts per &var.;
    proc means data=RS_file nway noprint;
    var FU_Exp FU_UnExp;
    class &var.;
    output out=RS_file2(drop=_:) sum=Tot_Exp Tot_UnExp;
    run;

    *Creating RiskSetID and StratumName;
    data RS_file3;
    merge RS_file
          RS_file2;
    by &var.;

    if _N_=1 then RiskSetID=0;

    *initialize total count in new variables (retain);
    if first.&var. then do;
      *this represents the total number of patients at the first FU time;
      totExposed=Tot_Exp;
      totUnexposed=Tot_UnExp;
    end;
      
    *ExposureProbability for up to this FU time;
    ExposureProbability=totExposed/(totExposed+totUnexposed);     

    format StratumName $40.;
    StratumName="&CurrentCat.";

	if "&var."="AllPts" then &var.="";

    *output events;
    do i=1 to FU_EVExp;
     Case_Exposure=1;
     RiskSetID=RiskSetID+1;
     output;
    end;

    do i=1 to FU_EVUnExp;
     Case_Exposure=0;
     RiskSetID=RiskSetID+1;
     output;
    end;

    *Deduct patients for this FU time (after output);
    totExposed=totExposed-FU_Exp;
    totUnexposed=totUnexposed-FU_UnExp;

    rename &var.= MatchId;
    retain RiskSetID totExposed totUnexposed;
    keep &var. Type Followuptime_ITT ExposureProbability Case_Exposure RiskSetID StratumName;      
    run;

    /**************************************/
    /* aggregate data for risk difference */
    /**************************************/

    *Compute sums by stratum and exposure;
    proc means data=RS_file0 nway noprint;
    var Exp UnExp EVExp EVUnExp FUTimeExp FUTimeUnExp;
    class &var.;*either All Pts or &matchid.;
    output out=RD_file(drop=_:) sum=;
    run;

    *one row per stratum;
		data stratatotals;
			set RD_file;
			N1sq = FUTimeExp**2;
			N0sq = FUTimeUnExp**2;
			numerator = (N1sq*N0sq);
			denominator =  (EVExp*N0Sq + EVUnExp*N1sq) ;
			weight = 0;
			if denominator > 0 then do;
                weight = numerator/denominator ; *1.;
			    diff = (EVExp/FUTimeExp) - (EVUnExp/FUTimeUnExp); *2.;
			    weighted_diff = weight*diff; *3.;
			end;
%if  "&var."="percentile"  %then %do;
      PercentileValue=input(percentile,best.);
%end;
%else %do;
      PercentileValue=0;
%end;
		run; 

    *create data for Risk Difference analysis (trivial for AllPts);
%if  "&var."="percentile"  %then %do;
    proc means data=stratatotals /*nway*/ noprint;
    var Exp UnExp EVExp EVUnExp FUTimeExp FUTimeUnExp
        weight weighted_diff;
class PercentileValue;
    output out=RD_file(drop=_:) sum=;
    run;



%end;
%else %do;
        proc means data=stratatotals nway noprint;
    var Exp UnExp EVExp EVUnExp FUTimeExp FUTimeUnExp
        weight weighted_diff;

    output out=RD_file(drop=_:) sum=;
    run;
%end;

data RD_file;
set RD_file;
if PercentileValue=. then PercentileValue=0;
run;

  /***********************************/
  /* aggregate data for Survival Data */
  /***********************************/

%if  "&var."="percentile"  %then %do;
    %let PercentileValue=PercentileValue;
    data RS_file0;
    set RS_file0;
    PercentileValue=input(percentile,best.);
    output;
    PercentileValue=0;  * 0 means all patients (all percentile patients, regardless of percentile values);
    output;
    run;
%end;
%else %do;
    %let PercentileValue=;
    data RS_file0;
    set RS_file0;
    PercentileValue=0;  * 0 means all patients (all percentile patients, regardless of percentile values);
    run;
%end;

%if %eval("&var." NE "AllPts" & "&var." NE "percentile") %then %do;
    *find maximum followuptime in each arm for each matchid;
    proc means data=RS_file0 nway noprint;
    var FUTimeExp FUTimeUnExp;
    class &var.;
    output out=MAX_FU_TIME(drop=_:) max=maxFUTimeExp maxFUTimeUnExp;
    run;

    data MAX_FU_TIME;
    set MAX_FU_TIME;
    stopFU=min(maxFUTimeExp,maxFUTimeUnExp);
    keep &var. stopFU;
    run;

    proc sql noprint undo_policy=none;
    create table RS_file0 as
    select dat.*,
           fu.stopFU
    from RS_file0 as dat left join
         MAX_FU_TIME as fu
    on dat.&var. = fu.&var.;
    quit; 

    data RS_file0;
    set RS_file0;

    if FUTimeExp > stopfu then do;
      FUTimeExp = stopfu;
      Followuptime_ITT= stopfu;
      EVExp=0;
    end;
    if FUTimeUNExp > stopfu then do;
      FUTimeUNExp = stopfu;
      Followuptime_ITT= stopfu;
      EVUnExp=0;
    end;
    run;
%end;


    *Create one record per Followuptime (collapse);
    proc means data=RS_file0 nway noprint;
    var Exp UnExp EVExp EVUnExp;
	class PercentileValue Followuptime_ITT;
    output out=KM_file(drop=_:) sum= FU_Exp FU_UnExp FU_EVExp FU_EVUnExp;
    run;

    *Total Number of Pts at time 0;
    proc means data=KM_file nway noprint;
    var FU_Exp FU_UnExp;
	class PercentileValue; 
    output out=KM_file2(drop=_:) sum=Tot_Exp Tot_UnExp;
    run;
    
    data KM_file3;
    set KM_file2(in=a)
        KM_file;
    by PercentileValue;
	if first.PercentileValue then do;
      *Set remaining pts;
      NExp=Tot_Exp;
      NUnExp=Tot_UnExp;
    end;
    else do;
	%if %eval("&var." NE "AllPts" & "&var." NE "percentile") %then %do;
	if NExp ne 0 and NUnExp ne 0 then output;
      NExp=NExp-FU_Exp;
      NUnExp=NUnExp-FU_UnExp;
    end; 
	*Deduct patients for this FU time (after output);
    rename  Followuptime_ITT=FollowUpDay
            FU_EVExp=EVExp
            FU_EVUnExp=EVUnExp;
    retain NExp NUnExp;
	%end;
	%if %eval("&var." = "AllPts" or "&var." = "percentile") %then %do;
	output;
      NExp=NExp-FU_Exp;
      NUnExp=NUnExp-FU_UnExp;
    end; 
	*Deduct patients for this FU time (after output);
    rename  Followuptime_ITT=FollowUpDay
            FU_EVExp=EVExp
            FU_EVUnExp=EVUnExp;
    retain NExp NUnExp;
	%end;
    keep NExp NUnExp FU_EVExp FU_EVUnExp Followuptime_ITT PercentileValue;
    run;

  %END;
  %ELSE %DO;
       	%put WARNING: no patients for covarnum &covarnum. and stratum &Stratvar. = &CurrentCat.;
  %END;

  *Programming note: The handling of empty datasets bewteen RS_FILE3 (or _RISKSETDATA_) and  RD_file/KM_file3 (or _RISKDIFFDATA_/_SURVIVALDATA_)
   are different as per the specifications.  The latter is the only one including a "MatchId" variable and it will be set to 
   missing in the case of AllPts and will be populated in the case of a valid "Matched" analysis.  When there are no data 
   in the matched analysis, then an empty line is generated with MatchID set to missing.  Because there are no MatchId 
   in the _RISKDIFFDATA_/_SURVIVALDATA_, the variable type=AllPts or Matched is necessary to distinguish between the two level of analysis;

  *Manage empty riskset dataset;
  %ISDATA(dataset=RS_file3);
  %IF %EVAL(&NOBS.=0) %THEN %DO;

    data  RS_StratumName;
    format StratumName $40.;
    StratumName="&CurrentCat.";
    run;
    
    data RS_file3;
    set RS_file0(obs=0 keep=&var. Type Followuptime_ITT) /*to copy &var. formats*/
        RS_StratumName;
    format ExposureProbability Case_Exposure RiskSetID best.;
    call missing(ExposureProbability,Case_Exposure,RiskSetID);
    rename &var.= MatchId;
    keep &var. Type Followuptime_ITT ExposureProbability Case_Exposure RiskSetID StratumName;      
    run;
  %END;

  *Manage empty riskDiff dataset;
  %ISDATA(dataset=RD_file);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    data RD_file;
    Exp=.;UnExp=.;EVExp=.;EVUnExp=.;FUTimeExp=.;FUTimeUnExp=.;weight=.; weighted_diff=.;PercentileValue=.;
    output;
    run;
  %END;

  *Manage empty survival dataset;
  %ISDATA(dataset=KM_file3);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    data KM_file3;
    NExp=.; NUnExp=.; EVExp=.; EVUnExp=.;FollowUpDay=.;PercentileValue=.;
    output;
    run;
  %END;

  /***********************************/
  /*appending subgroup to create one */
  /* final dataset per comp/look     */
  /***********************************/

  %ISDATA(dataset=&RUNID._RISKSETDATA_&COMP._&LOOK.);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    *first run for this combination; 
  	*Apply formating to variables. Must rename some of them to be able to reduce their length without w.a.r.n.i.n.g.;
    data &RUNID._RISKSETDATA_&COMP._&LOOK.;
    format covarnum 4. Type $20. StratumName $50. _MatchId $20. _Case_Exposure 1. _ExposureProbability 6.5 _Followuptime_ITT 6.;	
	length covarnum 4 _Case_Exposure 3 _ExposureProbability 6 _Followuptime_ITT 6;
	drop MatchId Case_Exposure ExposureProbability Followuptime_ITT;
	rename _MatchId=MatchId;
	rename _Case_Exposure=Case_Exposure;
	rename _ExposureProbability=ExposureProbability;
	rename _Followuptime_ITT=Followuptime_ITT;
    set RS_file3;
	_MatchId=MatchId;
	_Case_Exposure=Case_Exposure;
	_ExposureProbability=ExposureProbability;
	_Followuptime_ITT=Followuptime_ITT;
    covarnum=&covarnum.;
    label _MatchId="MatchId of the riskset &var.";
	label CovarNum="CovarNum";
	label RiskSetID="RiskSetID";
	label Type="Type";
	label StratumName="StratumName";
	label _Case_Exposure="Case_Exposure";
	label _ExposureProbability="ExposureProbability";
	label _Followuptime_ITT="Followuptime_ITT";

	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;
  %ELSE %DO;
  	*Apply formatting to variables;
  	data RS_file3_Fmt;
	format covarnum 4. Type $20. StratumName $50. _MatchId $20. _Case_Exposure 1. _ExposureProbability 6.5 _Followuptime_ITT 6.;	
	length covarnum 4 _Case_Exposure 3 _ExposureProbability 6 _Followuptime_ITT 6;
	drop MatchId Case_Exposure ExposureProbability Followuptime_ITT;
	rename _MatchId=MatchId;
	rename _Case_Exposure=Case_Exposure;
	rename _ExposureProbability=ExposureProbability;
	rename _Followuptime_ITT=Followuptime_ITT;
	set RS_file3;
	_MatchId=MatchId;
	_Case_Exposure=Case_Exposure;
	_ExposureProbability=ExposureProbability;
	_Followuptime_ITT=Followuptime_ITT;
    covarnum=&covarnum.;
	run;

    data &RUNID._RISKSETDATA_&COMP._&LOOK.;	
    set &RUNID._RISKSETDATA_&COMP._&LOOK.
        RS_file3_Fmt;
	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;

  %ISDATA(dataset=&RUNID._RISKDIFFDATA_&COMP._&LOOK.);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    *first time for this combination; 
    data &RUNID._RISKDIFFDATA_&COMP._&LOOK.;
    set RD_file;
    format covarnum 4. StratumName $40. Type $20.;
    covarnum=&covarnum.;
    StratumName="&CurrentCat.";
    if "&var." in ("AllPts","percentile")  then Type="AllPts ";
    else type="Matched";

	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;
  %ELSE %DO;
    data &RUNID._RISKDIFFDATA_&COMP._&LOOK.;
    set &RUNID._RISKDIFFDATA_&COMP._&LOOK.
        RD_file(in=a);
    if a then do;
      covarnum=&covarnum.;
      StratumName="&CurrentCat.";
    if "&var." in ("AllPts","percentile")  then Type="AllPts ";
      else type="Matched";
    end;
	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;

  %ISDATA(dataset=&RUNID._SURVIVALDATA_&COMP._&LOOK.);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    *first time for this combination; 
    data &RUNID._SURVIVALDATA_&COMP._&LOOK.;
    set KM_file3;
    format covarnum 4. StratumName $40. Type $20.;
    covarnum=&covarnum.;
    StratumName="&CurrentCat.";
    if "&var." in ("AllPts","percentile")  then Type="AllPts ";
    else type="Matched";

	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;
  %ELSE %DO;
    data &RUNID._SURVIVALDATA_&COMP._&LOOK.;
    set &RUNID._SURVIVALDATA_&COMP._&LOOK.
        KM_file3(in=a);
    if a then do;
      covarnum=&covarnum.;
      StratumName="&CurrentCat.";
    if "&var." in ("AllPts","percentile")  then Type="AllPts ";
      else type="Matched";
    end;
	*Relabel Type for unconditional since cannot stratify by matchid (&var.) in code above;
	if index(StratumName, "Unconditional") > 0 then Type="Matched";
    run;
  %END;

  proc datasets library=work nolist nowarn;
  delete RS_: RD_: KM_: stratatotals; 
  quit;

  %put NOTE: ********END OF MACRO: ms_CreateRiskSet v1.1********;

%mend ms_CreateRiskSet;

