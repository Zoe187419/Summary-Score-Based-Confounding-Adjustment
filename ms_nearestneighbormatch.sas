
****************************************************************************************************
*                                           PROGRAM OVERVIEW
****************************************************************************************************
*
* PROGRAM: ms_nearestneighbormatch.sas  
*
* Created (mm/dd/yyyy): 01/26/2015
* Last modified: 03/09/2016
* Version: 1.1
*
*--------------------------------------------------------------------------------------------------
* PURPOSE:
*  This program will perform propensity score nearest neighbor optimal matching between two groups
*
*  Program inputs:                                                                                   
*  -dataset containing the propensity scores  
*
* 
* PARAMETERS:
*  InFile=  The name of the dataset containing the propensity scores.
*  OutFile= The name of the dataset containing the matches and their controls.
*  MatchVars= The combination of variables that defines the records match level.
*  PSVar= The name of the variable containing the propensity scores.
*  MatchRatio= The maximum number of controls for one match.
*  WithReplacement= Matching with or without replacement.
*  Caliper= The maximum allowable distance between match and control propensity scores.
*
*  Program outputs:                                                                                                                                                            
*	-dataset with the matches and their controls 
*
*  Programming Notes:                                                                                
*
*
*--------------------------------------------------------------------------------------------------
* CONTACT INFO: 
*  xiaojuan_li@harvardpilgrim.org
*
*--------------------------------------------------------------------------------------------------
*  CHANGE LOG: 
*
*   Version   Date       Initials	   Comment (reference external documentation when available)
*   -------   --------   --------   ---------------------------------------------------------------
*   1.1       03/09/16   AH         Added Digits as a macro parameter instead of caluclating it 
*                                   from Caliper.
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

%macro ms_NearestNeighborMatch(InFile=,
                               OutFile=,
                               MatchVars=,
                               PSVar=,
                               MatchRatio=2,
                               WithReplacement=N,
                               Caliper=1.1,
                               Digits=0.0000000000000001);

    %put =====> MACRO CALLED: ms_NearestNeighborMatch v1.1;


  %macro findInCtrl(dir,key);
      rc=ictrl.setcur(key:&key.);
 *    dist=0;
      do while(CTRL_YN=0 and rc=0 /*and dist<MinDist*/); 
        rc=ictrl.&dir.();     *Skip Treatment placholders - will always start at a trt placeholder;
         if abs(Ctrl_PS - trt_PS) >= MinDist then leave; 
*        if floor(abs(Ctrl_PS - trt_PS)/&digits.)*&digits. >= MinDist then leave; 
      end;
  %mend;

  %macro findInTrt(dir,key);
      rc=itrt.setcur(key:&key.);
*      dist=0;
      do while(TRT_YN=0 and rc=0 /*and dist <MinDist*/);  
      rc=itrt.&dir.();     *Skip Control placholders  - will always start at a ctrl placeholder;
      if abs(Ctrl_PS - trt_PS) >= MinDist then leave; 
*     if floor(abs(Ctrl_PS - trt_PS)/&digits.)*&digits. >= MinDist then leave; 
      end;
  %mend;

  %macro AddToHeap(YesNoVar);
    dist = floor(abs(Ctrl_PS - trt_PS)/&digits.)*&digits.; 
    if &YesNoVar.=1 and dist < minDist then do;
        heap_Ctrl_PS=Ctrl_PS;
        heap_trt_PS=trt_PS;
        Rand=ranuni(123455);
        N=N+1;
        heap.add();  *output this pair of neighbors to the heap space;
    end;
  %mend;

  %macro DeleteFromHeap;
    *bring the hash iterator to a null pointer thereby
     allowing for the clear() method on the hash object itself to be succesful.;
    rc0=iheap.first();
    rc0=iheap.prev();
    rc0=heap.remove(); 
  %mend;

  %macro DeleteFrom(type1,type2,key);  *one key only;
    key0=&key.;
    *bring the hash iterator to a null pointer thereby
     allowing for the clear() method on the hash object itself to be succesful.;
    rc0=&type1..first();
    rc0=&type1..prev();
    rc0=&type2..remove(key: key0); 
  %mend;

  *Sort raw data;
  proc sort data=&InFile. out=_raw;
  by &PSVar. exposure;
  run;

  %LET match_Error=;

%do i=1 %to &MatchRatio.;

  data treated(rename=CTRL_N=TRT_CTRL_N)
       controls(rename=TRT_N=CTRL_TRT_N);
  set _raw;

  TRT_YN=exposure=1;
  CTRL_YN=exposure=0;

  lastExposure=lag(Exposure);

  if _N_=1 then do;
    TRT_N=1;
    CTRL_N=1;
  end;  
  else do;

    if exposure=1 then TRT_N=TRT_N+1;
    else if exposure=0 and lastExposure=1 then TRT_N=TRT_N+1;

    if exposure=0 then CTRL_N=CTRL_N+1;
    else if exposure=1 and lastExposure=0 then CTRL_N=CTRL_N+1;

  end;

  retain TRT_N CTRL_N;
  run;

* post process placeholders for direct access;

  data treated;
  set treated;
  by TRT_N;
  if _N_=1 then Set_num=0;
  if first.TRT_N;

  if TRT_YN=1 then Set_num=Set_num+1;


*  ForSort=ranuni(1250);
  rename &PSVar.=TRT_PS;
  retain Set_num;
  keep &MatchVars. &PSVar. TRT_N TRT_CTRL_N TRT_YN Set_num /*ForSort*/;
  run;

  data Controls;
  set Controls;
  by CTRL_N;
  if first.CTRL_N;
*  ForSort=ranuni(1251);
  rename &PSVar.=CTRL_PS;
  keep &MatchVars. &PSVar. CTRL_N CTRL_TRT_N CTRL_YN /*ForSort*/ exposure;
  run;

  

  data _NULL_;

  *Reading the "Augmented" Controls space into memory;
  if _N_= 1 then do;
      if 0 then set controls;
      declare hash ctrl(hashexp:20, dataset: "Controls", ordered: 'a');  *Note:  here speed increase with ordered - Ctrl_N sequential so no worries;
      declare hiter ictrl('ctrl');
      ctrl.defineKey('Ctrl_N');
      ctrl.defineData('Ctrl_N','Ctrl_PS','Ctrl_TRT_N','CTRL_YN','PatId','exposure');
      ctrl.defineDone();
      call missing(Ctrl_N,Ctrl_PS,Ctrl_TRT_N,CTRL_YN);
  end;

  *The Match space - Memory for Storing Final Selected Matches;
  if _N_= 1 then do;
     declare hash mat();
      mat.defineKey('Trt_N','MatchNumber'); 
      mat.defineData('Trt_N','Ctrl_N','MatchNumber');
      mat.defineDone();
      call missing(Trt_N,Ctrl_N,N,MatchNumber);
  end;



  *The Treated space - reading all Treated into Memory;
  if _N_= 1 then do;
      declare hash trt(hashexp:20, dataset: "Treated", ordered: 'a');  *Note:  here speed increase with ordered - Ctrl_N sequential so no worries;
      declare hiter itrt('trt');
      trt.defineKey('TRT_N');
      trt.defineData('TRT_N','TRT_PS','TRT_CTRL_N','TRT_YN');
      trt.defineDone();
      call missing(TRT_N,TRT_PS,TRT_CTRL_N,TRT_YN);
  end;

  *The Heap space - temporary space to store potential matches or matched pairs;
  if _N_= 1 then do;
      declare hash heap(ordered: 'a');
      heap.defineKey('Dist','Rand','N');                    *n is to ensure uniqueness; 
      heap.defineData('Dist','Rand','N','Trt_N','Ctrl_N','heap_trt_PS','heap_Ctrl_PS');  *rand is for ties;
      heap.defineDone();
      call missing(Dist,Ctrl_N,N,Rand,heap_trt_PS,heap_Ctrl_PS);
  end;

  *Variable initialization;
  N=0;
  MatchNumber=&i.;
  MinDist=&caliper.;
  MatchRatio=&MatchRatio.;

  *Fill the heap with neighbor pairs;
  rc1=itrt.first();
  do while(rc1=0);
    if TRT_YN=1 then do;
      %findInCTRL(prev,TRT_CTRL_N);  * means the control number in control for this treatment : TRT_CTRL_N;
      %AddToHeap(CTRL_YN);
      %findInCTRL(next,TRT_CTRL_N);
      %AddToHeap(CTRL_YN);
    end;
    rc1=itrt.next();
  end;

  *poll -Top-Down or sorted FILO);
  declare hiter iheap('heap');
  rc2=iheap.first();
  do while(rc2=0);

    *Are both trt_N and ctrl_N still available to create a pair?;
    IsTrtStillAvail=trt.check(key:TRT_N)=0;
    IsCtrlStillAvail=ctrl.check(key:CTRL_N)=0;
    dist1=.;
    dist2=.;


		*dist = abs(heap_trt_PS - heap_Ctrl_PS); 
*		dist = round(abs(heap_trt_PS - heap_Ctrl_PS),&digits.); 

    if IsTrtStillAvail=1 and IsCtrlStillAvail=1 /*and (dist < mindist)*/ then do;  *caliper was aldready checked by AddToHeap;
      rc0=mat.add();
      %DeleteFrom(itrt,trt,TRT_N);  *one key only;
      %DeleteFrom(ictrl,ctrl,Ctrl_N);  *one key only;
      %DeleteFromHeap;
    end;
    else if IsTrtStillAvail=1 and IsCtrlStillAvail=0 then do;  *caliper was aldready checked by AddToHeap;;

      *REFERENCE PS FOR DISTANCE CALCULATIONS;
      DeletedPS=heap_Ctrl_PS;

      *Control not available for matching;
      %DeleteFromHeap;

      rc=trt.find(key: TRT_N);  * load TRT_CTRL_N and TRT_PS into active memory;

      *Look left;
      %findInCTRL(prev,TRT_CTRL_N);
      if CTRL_YN=1 then do;
        dist1=abs(Ctrl_PS - DeletedPS);
        Ctrl_PS1=Ctrl_PS;
        Ctrl_N1=Ctrl_N;
      end;

      *Look Right;
      %findInCTRL(next,TRT_CTRL_N);
      if CTRL_YN=1 then do;
        dist2=abs(Ctrl_PS - DeletedPS);
        Ctrl_PS2=Ctrl_PS;
        Ctrl_N2=Ctrl_N;
      end;

      *which is better;
      if dist1 ne . and dist2 eq . then do;
*        dist=dist1;
        Ctrl_PS=Ctrl_PS1;
        Ctrl_N=Ctrl_N1;
      end;
      else if dist1 eq . and dist2 ne . then do;
*        dist=dist2;
        Ctrl_PS=Ctrl_PS2;
        Ctrl_N=Ctrl_N2;
      end;
      else if . < dist1 <= dist2 then do;
*        dist=dist1;
        Ctrl_PS=Ctrl_PS1;
        Ctrl_N=Ctrl_N1;
      end;
      else if . < dist2 < dist1 then do;
*        dist=dist2;
        Ctrl_PS=Ctrl_PS2;
        Ctrl_N=Ctrl_N2;
      end;

    if dist1 ne . or dist2 ne . then do;
      CTRL_YN=1;
      %AddToHeap(CTRL_YN);
    end; 
    end;
    else if IsTrtStillAvail=0 and IsCtrlStillAvail=1 then do;  *caliper was aldready checked by AddToHeap;;

      *REFERENCE PS FOR DISTANCE CALCULATIONS;
      DeletedPS=heap_trt_PS;

      *treated not available for matching;
      %DeleteFromHeap;

      rc=ctrl.find(key: CTRL_N);  * load CTRL_TRT_N and CTRL_PS into active memory;

      *Look left;
      %findInTRT(prev,CTRL_TRT_N);
      if TRT_YN=1 then do;
        dist1=abs(Trt_PS - DeletedPS);
        Trt_PS1 =Trt_PS; 
        trt_N1=trt_N;
      end;

      *Look Right;
      %findInTRT(next,CTRL_TRT_N);
      if TRT_YN=1 then do;
        dist2=abs(Trt_PS - DeletedPS);
        Trt_PS2 =Trt_PS; 
        trt_N2=Trt_N;
      end;

      *which is better;
      if dist1 ne . and dist2 eq . then do;
*        dist=dist1;
        trt_PS=trt_PS1;
        trt_N=trt_N1;
      end;
      else if dist1 eq . and dist2 ne . then do;
*        dist=dist2;
        trt_PS=trt_PS2;
        trt_N=trt_N2;
      end;
      else if . < dist1 <= dist2 then do;
*        dist=dist1;
        trt_PS=trt_PS1;
        trt_N=trt_N1;
      end;
      else if . < dist2 < dist1 then do;
*        dist=dist2;
        trt_PS=trt_PS2;
        trt_N=trt_N2;
      end;

      if dist1 ne . or dist2 ne . then do;
        TRT_YN=1;
        %AddToHeap(TRT_YN);
      end;
    end;
    else do;
      %DeleteFromHeap;
    end;

    *poll;
    rc2 = iheap.first();          
  end;
  /*
  itrt.delete();
  trt.delete();
  iheap.delete();
  heap.delete();
*/
  mat.output(dataset:"res1"); 
  %if %eval(&MatchRatio. >1  and &i. <&MatchRatio.) %then %do;
    ctrl.output(dataset:"ctrl"); 
  %end;
  run;

*Reassinging PatIds and PS values;
  *Assign Matchvars of treated;
  proc sql noprint undo_policy=none;
  create table res2 as
  select base.*,
         match.ctrl_N,
         match.MatchNumber
  from Treated as base,
       res1 as match
  where match.Trt_N=base.Trt_N;
  quit;

  *Assign Matchvars of controls;
  proc sql noprint undo_policy=none;
  create table _outfile(drop=Trt_N Ctrl_N) as
  select match.*,
         ctrl.patid as ctrl_Patid,
         ctrl.Ctrl_PS,
         abs(Ctrl_PS-Trt_PS) as match_distance
  from res2 as match,
       Controls as ctrl
  where ctrl.ctrl_N=match.ctrl_N
  order by Patid, MatchNumber;
  quit;


  %if &i.=1 %then %do;*overwrite;
   data &outfile.;
   set _outfile;
   run;
  %end;
  %else %do;
   data &outfile.;
   set &outfile. _outfile;
   by Patid MatchNumber;
   run;
  %end;


  %if %eval(&MatchRatio. >1  and &i. <&MatchRatio.) %then %do;
    *Assign Matchvars of controls;

  	data _null_;
  	dsid=open("ctrl");
  	call symputx("NOBSCTRL",attrn(dsid,"NLOBS"));
  	run;
    %if &NOBSCTRL.=0 %then %do;
      %let i=&MatchRatio.+1;
    %end;
    %else %do;
      data _raw;
      set _raw(where=(exposure=1))
          ctrl(where=(exposure=0) rename=Ctrl_PS=&PSVar.);
      by &PSVar. exposure;
      *if b then exposure=0;
      keep PatId &PSVar. exposure;
      run;
    %end;  
  %end;

*Cleaning up in case of multiple runs;
    proc datasets library=work nowarn nolist;
    delete _outfile;
    *delete Treated Controls res1 res2;
    quit;
%end;

*if outfile is empty, no matches were found;
  %ISDATA(dataset=&outfile.);
  %IF %EVAL(&NOBS.=0) %THEN %DO;
    %let match_Error=1;
  %END;

    %put NOTE: ********END OF MACRO: ms_NearestNeighborMatch v1.1********;


%mend;
