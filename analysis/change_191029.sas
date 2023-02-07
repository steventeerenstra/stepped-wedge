* note: home > ward > nurse > measurement time, but most nurses only 1 measurement , so home > ward;
* note: each measurement time has several moments observed, but in this analysis these are analysed separately;

* hard coding of three levels, make better later in sw_bin_infer;
* using " " for making labels in the subgroup statement, else macro vars not resolved
label &outcome="&outcome_lab | subgroup &counter : &subgroup_name";
* using the wording: "analysis with interaction term: &interaction";

* to do still in continuous version: no output if model not converged;

* questions: how comes that the observed versus predicted plots on cluster level look so good?;
* I would expect that the shape must be same over time, and only the height changes as only a random effect for cluster;

/* to automate writing contrast and estimate statements: use x=" "; xx = cat(repeat(' 0',5), x);
to make a series of _0_0_0 etc; and then substr(xx,4, 1)=1; to replace a 0 into a 1
*/

**** MACRO DEFINITIONS ****;

%macro rm_design(ds=ds, cluster=cluster,trt=trt, period=period);
* general check on repeated measures design;
title4 "treatments applied per cluster-period";
title4 "check whether each cluster-period has one treatment?";
proc tabulate data=&ds;class &cluster &period &trt;format &period 4.0;
table &cluster*&trt, &period*N*f=4. /rts=15 ;
run;
%mend;


* continuous outcome";
* descriptives;
%macro sw_cont_descrip(ds=ds,outcome=, cluster=cluster, period=period, trt=trt,
outcome_min=-1000, outcome_max=+1000,period_min=1,period_max=14,dsdescrip=);

*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;
title3 "description   ****&outcome: &outcome_label *****";
title4 "all values";
proc means data=&ds n nmiss min max p5 p95 Q1 median Q3 mean std; var &outcome; run;
proc sgpanel data=&ds;
panelby &cluster / onepanel;
histogram &outcome;
run;
title4 "data size per cluster-period";
proc tabulate data=&ds; class &cluster &period; var &outcome;
table &cluster, &period="&outcome: N per cluster-period"*&outcome=' '*(N=' ')*f=5.0 / rts=10;  run;
title4 "mean per cluster-period";
proc tabulate data=&ds; class &cluster &period; var &outcome;
table &cluster, &period="&outcome: mean per cluster-period"*&outcome=' '*(mean=' ')*f=5.2 / rts=10;  run;

title3 "analysis    ****&outcome: &outcome_label ****";
title4 "time series plot per cluster";
proc means data=&ds noprint; id &trt; class &cluster &period;
var &outcome; output out=_ClusPeriodMeans mean(&outcome)=mean;run;
proc sgpanel data=_clusPeriodmeans;
panelby &cluster / onepanel; *or columns=1 to get only one column;
series x=&period y=mean;
scatter x=&period y=mean /group=trt markerattrs=(size=9) ;
rowaxis integer min=&outcome_min max=&outcome_max alternate;
colaxis min=&period_min max=&period_max alternate integer;
run;
* to do: time series analysis per cluster using standard techniques?*;
title4 "disregarding possible time trend:";
title5 "shift in distribution by intervention";
proc means data=&ds n nmiss min max p5 p95 Q1 median Q3 mean std; class &trt;var &outcome; run;
data _shift; set &ds;
if &trt=1 then intervention_&outcome=&outcome;
if &trt=0 then control_&outcome=&outcome;
run;
* shift over all clusters;
proc sgplot data=_shift;
histogram intervention_&outcome /transparency=0.3;
histogram control_&outcome /transparency=0.6;
run;
* shift by cluster;
proc sgpanel data=_shift;
panelby &cluster / onepanel;
histogram intervention_&outcome /transparency=0.3;
histogram control_&outcome /transparency=0.6;
run;
title5 "shift in means, medians, boxplots";
proc sgpanel data=&ds;
panelby &cluster / onepanel;
*vbar &trt / response=&outcome stat=mean;
hbox &outcome / category=&trt;
run;
*show descriptives and if need save them***;
proc means data=&ds n nmiss mean std median Q1 Q3; var &outcome; class &trt;
%IF &dsdescrip ne %THEN %DO;
output out=_descrip n=n nmiss=nmiss mean=mean std=std median=median Q1=Q1 Q3=Q3 ;
%END;
run;
%IF &dsdescrip ne %THEN %DO;
data _descrip; length outcome $30; length outcome_label $ 100; length stat $ 10;
set _descrip; outcome="&outcome";outcome_label="&outcome_label";
array aux[7] n nmiss mean std median Q1 Q3;
do i=1 to 7;
stat=vname( aux[i]); value=aux[i];output;
end;
drop  n nmiss mean std median Q1 Q3;
if _type_=1;
keep outcome outcome_label &trt stat value;
run;
proc append force base=&dsdescrip data=_descrip;run;
%END;
%mend sw_cont_descrip;

%macro sw_cont_infer(ds=ds,outcome=, cluster=cluster, subcluster=subcluster, period=period, trt=trt,
period_min=1,period_max=14,covars_cat=, covars=, interaction=,
dseff=, random_statement=%str(random intercept /subject=&cluster;random intercept /subject=&subcluster(&cluster); ),
repeated_statement=%str( ), estimate_statement=%str( ), contrast_statement=%str( )
);
* to do: ITS analyse by cluster;
* to do: Hooper-Girling, exponential decay model;
* to do: different size for observed vs predicted? GTL scatter plot has an option for that;

*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;
title3 "analysis   **** &outcome: &outcome_label ****";
title4 "modeling categorical time trend, random effect cluster (Hussey & Hughes)";
* if dseff asked then save effects for &trt and possibly estimate from estimate statement;
%IF &dseff ne %THEN
%DO;
ods output SolutionF=_solutions;
    %IF &estimate_statement ne %THEN ods output Estimates=_estimates;;
    %IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
    * if interaction asked then also make interaction term;
    %local interaction_term; %let interaction_term=;
    %IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;
****** Hussey and Hughes model ****;
proc mixed data=&ds;
class &cluster &subcluster &period &covars_cat ;
model &outcome= &trt &period &covars &interaction_term/ solution cl outpred=_respred ;
&random_statement;;
&repeated_statement;;
&estimate_statement;;
&contrast_statement;;
run;
* process the effects;
%IF &dseff ne %THEN
%DO;
ods output close;
data _solutions;
    length effect $50;length outcome $30;length outcome_label $100;length model $100;
    set _solutions;
    if index(effect,"&trt")>0;
    outcome="&outcome";outcome_label="&outcome_label";
    * get the name of the level of the interaction variable from the format;
    %IF &interaction ne %THEN effect=catx(' : ',effect,&interaction, vvalue(&interaction));;
    model="&trt &period &covars &interaction_term // &random_statement &repeated_statement";
    rename effect=label;
    keep outcome outcome_label effect estimate stderr probt lower upper model ;
    run;
    proc append force base=&dseff data=_solutions ;run;
%IF &estimate_statement ne %THEN %DO;
    data _estimates;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _estimates;
    outcome="&outcome";outcome_label="&outcome_label";
    model="&trt &period &covars &interaction_term // &random_statement &repeated_statement";
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
    proc append force base=&dseff data=_estimates;run;
    %END;
%IF &contrast_statement ne %THEN %DO;
    data _contrasts;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _contrasts;
    outcome="&outcome";outcome_label="&outcome_label";
    model="&trt &period &covars &interaction_term // &random_statement &repeated_statement";
    probt=probf; estimate=.;stderr=.;lower=.;upper=.;
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
    proc append force base=&dseff data=_contrasts;run;
    %END;
%END;
*** end Hussey & Hughes model**;

title6 "residual plot: all values";
proc sgplot data=_respred; scatter x=pred y=resid ; refline 0/axis=y;run;
title6 "residual plot: by cluster";
proc sgpanel data=_respred; panelby &cluster /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
title6 "residual plot: by period";
proc sgpanel data=_respred; panelby &period /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
/* too much detail
title6 "residual plot: by cluster-period";
proc sgpanel data=_respred; panelby &cluster &period ; scatter x=pred y=resid; refline 0 / axis=y;run;
*/
* fit, residuals and observed by cluster;
proc means data=_respred noprint; class &cluster &period; id &trt;
output out=_pred1 mean(pred &outcome )=mean_pred mean_obs ;run;
data _pred1; set _pred1; mean_resid=mean_obs-mean_pred;run;
title6 "observed vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
series x=&period y=mean_pred;
scatter x=&period y=mean_obs;
rowaxis integer;
run;
title6 "residual vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot of *all* clusters";
proc sgplot data=_pred1;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
%mend sw_cont_infer;

%macro sw_cont(ds=ds,outcome=, cluster=cluster, subcluster=subcluster, period=period, trt=trt,
outcome_min=-1000, outcome_max=+1000,period_min=1,period_max=14,
covars_cat=, covars=, interaction=,
dseff=, dsdescrip=,random_statement=%str(random intercept /subject=&cluster; ),
repeated_statement=%str( ), estimate_statement=%str( ), contrast_statement=%str( )
);
* descriptives;
%sw_cont_descrip(ds=&ds,outcome=&outcome, cluster=&cluster, period=&period, trt=&trt,
outcome_min=&outcome_min, outcome_max=&outcome_max,period_min=&period_min,period_max=&period_max,
dsdescrip=&dsdescrip
);
*inference;
%sw_cont_infer(ds=&ds,outcome=&outcome, cluster=&cluster, subcluster=&subcluster, period=&period, trt=&trt,
period_min=&period_min,period_max=&period_max,covars_cat=&covars_cat, covars=&covars, interaction=&interaction,
dseff=&dseff, random_statement=&random_statement,
repeated_statement=&repeated_statement, estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%mend sw_cont;

******;
******;
******;
%macro sw_bin_descrip(ds=ds, outcome=,cluster=cluster, period=period, trt=trt,period_min=1, period_max=14,
dsdescrip=
);
** assumes variable is {0,1} coded;

data &ds; set &ds;
pct_&outcome=100*&outcome;* for calculating percentages via a means-operator;
dummy=1;* for having a categorical variable with only one value;
run;
*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;
title3 "description   **** &outcome: &outcome_label ****";
title4 "all values";
proc freq data=&ds; table &outcome /missing; run;
title4 "percentage &outcome=1 (assuming 0,1 coding)";
proc format; value dummyft 1=' ';run;
proc sgpanel data=&ds;format dummy dummyft.;
panelby &cluster / onepanel;
* need dummy variable to have a category variable for vbar statement;
vbar dummy /  stat=mean response=pct_&outcome;
colaxis label=' ';
run;
title4 "data size per cluster-period";
proc tabulate data=&ds; class &cluster &period; var &outcome;
table &cluster, &period="&outcome: N per cluster-period"*&outcome=' '*(N=' ')*f=5.0 / rts=10;  run;
title4 "percent per cluster-period";
proc tabulate data=&ds; class &cluster &period; var pct_&outcome;
table &cluster, &period="&outcome: percent per cluster-period"*pct_&outcome=' '*(mean=' ')*f=5.2 / rts=10;run;


title3 "analysis   **** &outcome: &outcome_label ****";
title4 "time series plot per cluster";
proc means data=&ds noprint; id &trt; class &cluster &period;
var pct_&outcome; output out=_ClusPeriodMeans mean(pct_&outcome)=percent;run;
proc sgpanel data=_clusPeriodmeans;
panelby &cluster / onepanel; *or columns=1 to get only one column;
series x=&period y=percent;
scatter x=&period y=percent /group=trt markerattrs=(size=9) ;
rowaxis min=0 max=100 alternate;
colaxis min=&period_min max=&period_max alternate integer;
run;
* to do: time series analysis per cluster using standard techniques?*;
title4 "disregarding possible time trend:";
title5 "shift in percentage by intervention";
proc sgpanel data=&ds;
panelby &cluster / onepanel;
vbar &trt / response=pct_&outcome stat=mean;
run;
title6 "descriptives with missing data";
proc freq data=&ds; table &outcome*&trt / missing norow nopercent nocol; run;
title6 "descriptives with only non-missing data";
proc freq data=&ds; table &outcome*&trt / norow nopercent; run;
* save percentages (n, ntot, nmissing) for non-missing data (n, ntot, perc, nmiss per treatment group);
* if asked;
%IF &dsdescrip ne %THEN %DO;
proc means data=&ds noprint; class &trt; var &outcome; 
output out=_descrip mean=mean sum=sum n=n nmiss=nmiss;
run;

data _descrip; length stat $ 10;length outcome $30; length outcome_label $ 100;
set _descrip(drop=_type_ _freq_ where=(&trt ne .));
outcome="&outcome";outcome_label="&outcome_label";
*dataset that is set contains one record for each treatment;
stat="perc"; value=100*mean; output;
stat="n";value=sum;output;
stat="ntot";value=n;output;
stat="nmiss";value=nmiss;output;
keep outcome outcome_label &trt stat value;
run;
proc append force base=&dsdescrip data=_descrip;run;
%END;
%mend sw_bin_descrip;


%macro sw_bin_infer(ds=ds, outcome=,cluster=cluster, subcluster=subcluster, period=period, trt=trt,period_min=1, period_max=14,
random_statement=%str( ),repeated_statement=%str( ),
estimate_statement=%str( ), 
contrast_statement=%str( ),dseff=,covars_cat=, covars=, interaction=,
refcat=first
);
** reference category for outcome is &refcat
** works if variable is {0,1} coded;
** error is printed if &outcome is 0 in both arms of is 1 in both arms, but this is handled ;

data &ds; set &ds;
pct_&outcome=100*&outcome;* for calculating percentages via a means-operator;
run;
*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;

title3 "analysis   **** &outcome: &outcome_label ****";
title4 "percent point difference (linear mixed model)";
* if dseff asked then save effects for &trt and possibly estimate from estimate statement;
%IF &dseff ne %THEN
%DO;
ods output SolutionF=_solutions;
    %IF &estimate_statement ne %THEN ods output Estimates=_estimates;;
    %IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
    * if interaction asked then also make interaction term;
    %local interaction_term; %let interaction_term=;
    %IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;


****** % difference: Hussey and Hughes model ****;
title5 "LMM: model &outcome =&trt  &period &covars &interaction_term; class &cluster &period &covars_cat";
ods output ConvergenceStatus=_converged;*update convergence status;
proc mixed data=&ds;
class &cluster &subcluster &period &covars_cat ;
model pct_&outcome= &trt &period &covars &interaction_term/ solution cl outpred=_respred ;
random intercept / subject=&cluster;
random intercept / subject=&subcluster(&cluster);
&repeated_statement;;
&estimate_statement;;
&contrast_statement;;
run;
ods output close; * close the ods for at least the convergence checking;

* depending on convergence status, status=0 means convergenced, estimates, plots etc are set missing;
* first assume that the model did not convergence;
%local nonconverged; %let nonconverged=1;
* update that after running the model, provided the model was able to run;
data _null_;  call symput('exist', exist("_converged") ); run; 
%IF &exist %THEN %DO; data _null_; set _converged; call symput('nonconverged', status);run;%END;

* process the effects;
%IF &dseff ne %THEN
%DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _solutions; length label $50;length outcome $30;length outcome_label $100;length model $100;
	outcome="&outcome"; outcome_label="&outcome_label";label="no &trt (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;
	model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";run;
	%END;
	%IF &nonconverged eq 0 %THEN %DO;
	data _solutions;
    length effect $50;length outcome $30;length outcome_label $100;length model $100;
    set _solutions;
    if index(effect,"&trt")>0;
    outcome="&outcome";outcome_label="&outcome_label";
    * get the name of the level of the interaction variable from the format;
    	%IF &interaction ne %THEN effect=catx(' : ',effect,&interaction, vvalue(&interaction));;
    model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";
    rename effect=label;
    keep outcome outcome_label effect estimate stderr probt lower upper model ;
    run;
	%END;
    proc append force base=&dseff data=_solutions ;run;
%END;

%IF &estimate_statement ne %THEN 
%DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _estimates;length label $50;length outcome $30;length outcome_label $100; length model $100;
	outcome="&outcome"; outcome_label="&outcome_label";label="no estimate (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;
	model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";
	run;
	%END;
	%IF &nonconverged eq 0 %THEN %DO;
	data _estimates;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _estimates;
    outcome="&outcome";outcome_label="&outcome_label";
    model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
    %END;
    proc append force base=&dseff data=_estimates;run;
%END;

%IF &contrast_statement ne %THEN 
%DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _contrasts;length label $50;length outcome $30;length outcome_label $100; length model $100;
	outcome="&outcome"; outcome_label="&outcome_label";label="no contrast (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;
	model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";
	run;
	%END;
	%IF &nonconverged eq 0 %THEN %DO;
	data _contrasts;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _contrasts;
    outcome="&outcome";outcome_label="&outcome_label";
    model="LMM: &trt &period &covars &interaction_term // &random_statement &repeated_statement";
    probt=probf; estimate=.;stderr=.;lower=.;upper=.;
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
	%END;
    proc append force base=&dseff data=_contrasts;run;
%END;


* only if model converged, provide plots;
%IF &nonconverged eq 0 %THEN %DO;
title6 "residual plot: all values";
proc sgplot data=_respred; scatter x=pred y=resid ; refline 0/axis=y;run;
title6 "residual plot: by cluster";
proc sgpanel data=_respred; panelby &cluster /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
title6 "residual plot: by period";
proc sgpanel data=_respred; panelby &period /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;

* calculate residuals, observed and fitted for cluster x time averages;
proc means data=_respred noprint; class &cluster &period; id &trt;
output out=_pred1 mean(pred pct_&outcome)=mean_pred mean_obs;run;
data _pred1; set _pred1; mean_resid=mean_obs-mean_pred;run;
title6 "observed vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
series x=&period y=mean_pred;
scatter x=&period y=mean_obs;
rowaxis integer;
run;
title6 "residual vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot of *all* clusters";
proc sgplot data=_pred1;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
%END;
*** end % difference: Hussey & Hughes model**;

title4 "** odds ratio ** (generalized linear mixed model)";
* first remove the convergence status from the previous model;
proc delete data=_converged; run;

%IF &dseff ne %THEN
%DO;
ods output ParameterEstimates=_solutions;
    %IF &estimate_statement ne %THEN ods output Estimates=_estimates;;
    %IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
    * if interaction asked then also make interaction term;
    %local interaction_term; %let interaction_term=;
    %IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;
*model;
title5 "GLMM: model &outcome(reference=&refcat)=&trt  &period &covars &interaction_term; class &cluster &period &covars_cat";
ods output ConvergenceStatus=_converged;*check on convergence status;
proc glimmix data=&ds;
class &cluster &subcluster &period &covars_cat;
model &outcome(reference=first)=&trt  &period &covars &interaction_term
/solution cl distribution=binary oddsratio; *&period &covars oddsratio;
random intercept / subject=&cluster;
random intercept / subject=&subcluster(&cluster);
&estimate_statement;;
&contrast_statement;;
output out=_respred pred=lp;
run;
ods output close; * close at least convergence checking;

* depending on convergence status, status=0 means convergenced, estimates, plots etc are set missing;
* assume first the model did not converge;
%let nonconverged=1;
* update that after running the model, provided the model was able to run;
data _null_;  call symput('exist', exist("_converged") ); run; 
%IF &exist %THEN %DO; data _null_; set _converged; call symput('nonconverged', status);run;%END;

* process the effects;
%IF &dseff ne %THEN
%DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _solutions;length label $50;length outcome $30;length outcome_label $100;length model $100;
	outcome="&outcome"; outcome_label="log(odds) &outcome_label";label="no &trt (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;
	model="GLMM: &trt &period &covars &interaction_term | &random_statement ";run;
	%END;
	%IF &nonconverged eq 0 %THEN %DO; 
	data _solutions;
    length effect $50;length outcome $30;length outcome_label $100;length model $100;
    set _solutions;
    if index(effect,"&trt")>0;
    outcome="&outcome";outcome_label="log(odds) &outcome_label";
    * get the name of the level of the interaction variable from the format;
    	%IF &interaction ne %THEN effect=catx(' : ',effect,&interaction, vvalue(&interaction));;
    model="GLMM: &trt &period &covars &interaction_term | &random_statement ";
    rename effect=label;
    keep outcome outcome_label effect estimate stderr probt lower upper model ;
    run;
	%END;
    proc append force base=&dseff data=_solutions ;run;
%END;

%IF &estimate_statement ne %THEN 
%DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _estimates;length label $50;length outcome $30;length outcome_label $100; length model $100;
	outcome="&outcome"; outcome_label="log(odds) &outcome_label";label="no estimate (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;
	model="GLMM: &trt &period &covars &interaction_term | &random_statement ";run;
	%END;  
	%IF &nonconverged eq 0 %THEN %DO; 
	data _estimates;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _estimates;
    outcome="&outcome";outcome_label="log(odds) &outcome_label";
    model="GLMM: &trt &period &covars &interaction_term | &random_statement ";
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
	%END;
    proc append force base=&dseff data=_estimates;run;
%END;

%IF &contrast_statement ne %THEN %DO;
	%IF &nonconverged ne 0 %THEN %DO; 
	data _contrasts;length label $50;length outcome $30;length outcome_label $100; length model $100;
	outcome="&outcome"; outcome_label="log(odds) &outcome_label";label="no contrast (model not converged)"; 
	estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
	model="GLMM: &trt &period &covars &interaction_term | &random_statement ";run;
	%END;
	%IF &nonconverged eq 0 %THEN %DO;
	data _contrasts;
    length label $50;length outcome $30;length outcome_label $100; length model $100;
    set _contrasts;
    outcome="&outcome";outcome_label="log(odds) &outcome_label";
    model="GLMM: &trt &period &covars &interaction_term // &random_statement &repeated_statement";
    probt=probf; estimate=.;stderr=.;lower=.;upper=.;
    keep outcome outcome_label label estimate stderr probt lower upper model;
    run;
	%END;
    proc append force base=&dseff data=_contrasts;run;
%END;
*** end odds ratio**;

* only if model converged, provide plots;
%IF &nonconverged eq 0 %THEN %DO;
title6 "observed vs predicted plot: by cluster";
* calculate residuals, observed and fitted for cluster x time averages;
data _respred; set _respred; pred_percent=100*exp(lp)/(1+exp(lp)); run;
proc means data=_respred noprint; class &cluster &period;
output out=_pred1 mean(pred_percent pct_&outcome)=mean_pred mean_obs;run;
data _pred1; set _pred1; mean_resid=mean_obs-mean_pred;run;
title6 "observed vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
series x=&period y=mean_pred;
scatter x=&period y=mean_obs;
rowaxis integer;
run;
title6 "residual vs predicted plot of clusters: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot of *all* clusters";
proc sgplot data=_pred1;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
%END;
%mend sw_bin_infer;

%macro sw_bin(ds=ds, outcome=,cluster=cluster, subcluster=subcluster, period=period, trt=trt,period_min=1, period_max=14,
random_statement=%str(random intercept /subject=&cluster; ),repeated_statement=%str( ),
estimate_statement=%str(estimate 'ctl' intercept 1 period  0 0 0 0 0 0 0 0 0 0 0 0 0 1 / cl; estimate 'exp' intercept 1 trt 1 period  0 0 0 0 0 0 0 0 0 0 0 0 0 1 / cl),
contrast_statement=%str( ),
dseff=,dsdescrip=, covars_cat=, covars=, interaction=
);
* descrip;
%sw_bin_descrip(ds=&ds, outcome=&outcome,cluster=&cluster, period=&period, trt=&trt,period_min=&period_min,
period_max=&period_max, dsdescrip=&dsdescrip
);
* inference;
%sw_bin_infer(ds=&ds, outcome=&outcome,cluster=&cluster, subcluster=&subcluster, period=&period, trt=&trt,period_min=&period_min,
period_max=&period_max, random_statement=&random_statement,repeated_statement=&repeated_statement,
estimate_statement=&estimate_statement,contrast_statement=&contrast_statement,
dseff=&dseff,covars_cat=&covars_cat, covars=&covars,
interaction=&interaction
);
%mend sw_bin;




* subgroup analyses for binary;
%macro sw_subgr_bin(ds=ds, cat_min=0, cat_max=3, outcome=dlr,
interaction=admission_type, estimate_statement=, contrast_statement=, dsdescrip=table4_bin, dseff=table4_infer);
* requires numeric variable for subgroup and possibly with a value label;
%DO counter=&cat_min  %TO &cat_max;
data _null; set &ds;where &interaction=&counter;
call symput('outcome_lab', vlabel(&outcome)); *variable label;
call symput('subgroup_name',vvalue(&interaction));* value label;
run;
data ds1; set &ds;where &interaction=&counter;
label &outcome="&outcome_lab | subgroup &counter : &subgroup_name";run;
%sw_bin_descrip(ds=ds1,outcome=&outcome, dsdescrip=&dsdescrip);
%END;
title2 "analysis with interaction term: &interaction";
%sw_bin_infer(ds=&ds,outcome=&outcome, covars_cat=&interaction, covars=&interaction, interaction=&interaction,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement,
dseff=&dseff);
%mend;
* per moment versies van een variabele maken;
%macro makebymoment(var=);
%do i=1 %to 5;
* for moment &i make the variable that takes the value of that variable only at moment &i;
if M&i=1 then &var._M&i= &var;
else &var._M&i =.;
%end;
%mend makebymoment;



********* READ 	DATA *******;
options pagesize=60;
libname dir ".";
* get formats;
options fmtsearch= (dir.formats);
* show variables in dataset;
title6 "variables in ruwe dataset";
proc contents data=dir.change;run;

* show format, look up name of format in the formatslibrary;
/* 1‚Helpende; 2‚Leerling verpleegkundige, 3‚Overig, 4‚Serviceassistent,
5‚Verpleegkundige, 6‚Verzorgende */ 
proc format lib=dir.formats fmtlib;
   select BEROEPS; 
run;
/* show format id_afde*/
proc format lib=dir.formats fmtlib;
  select ID_AFDE;
run;



proc format; 
value beroep 1="verpleegkundige+verzorgende"
                          2="helpende+serviceassistent"
						  3="lerende"
;
value lerende 0="niet-leerling"
			  1="leerling"	
;
value type 1="revalidatie"
           2="psy.geriatrisch"
		   3="somatisch"
;
value num 0="0=nee"
		  1="1=ja"
;
run; 

   
* get dataset;
* define variables;
data ds; set dir.change; 
* moment van handhygiene;
M1=Voorpati__ntcontact; label M1="M1: voor patientcontact";
M2=Vooraseptischehandeling;label M2="M2: voor aseptische beh.";
M3=Nacontactlichaamsvloeistoffen;label M3="M3: na contact lich.vloeistof";
M4=Napati__ntcontact;label M4="M4: na patientcontact";
M5=Nacontactomgevingpati__nt;label M5="M5: na contact omgev.pat.";
* levels;
cluster=id_verpleeghuis;
subcluster=id_afdeling;
period = meting;
* intervention;
trt= (meting >= tijdstipinterventie); 
********;
*0-1 format instead of "ja/nee";
format HH_uitgev_gecorr_asep num.;
format HH_uitgevoerd_M1 num.;format HH_uitgevoerd_M2 num.;
format HH_uitgevoerd_M3 num.;format HH_uitgevoerd_M4 num.;
format HH_uitgevoerd_M5 num.;
* derived outcomes;
* over alle metingen heen;
if HH_uitgev_gecorr_asep=1 and handenwassen=1 then HH_asep_wassen=1; else  HH_asep_wassen=0;
if HH_uitgev_gecorr_asep=1 and Handdesinfectie=1 then HH_asep_desinfec=1; else  HH_asep_desinfec=0;
if HH_uitgev_gecorr_asep=1 and Handschoenen=1 then HH_asep_handschoen=1; else  HH_asep_handschoen=0;
*	handenschoen gecombineerd met ofwel handenwasen ofwel handdesinfectie ;
if HH_uitgev_gecorr_asep=1 and Handschoenen=1 and handenwassen=1 then HH_asep_hschoen_wassen=1; else  HH_asep_hschoen_wassen=0;
if HH_uitgev_gecorr_asep=1 and Handschoenen=1 and Handdesinfectie=1 then HH_asep_hschoen_des=1; else  HH_asep_hschoen_des=0;
* per meting
*     n.b. geteld alleen binnen de observaties waarvoor dat hygiene moment geindiceerd was;
%makebymoment(var=HH_asep_desinfec);
%makebymoment(var=HH_asep_wassen);
%makebymoment(var=HH_asep_handschoen);
%makebymoment(var=HH_asep_hschoen_wassen);
%makebymoment(var=HH_asep_hschoen_des);
* gecorrigeerde HH_uitgev_gecorr_asep alleen op M2  ;
if M2=1 then HH_uitgev_gecorr_asep_M2=HH_uitgev_gecorr_asep; else HH_uitgev_gecorr_asep_M2=.;
*********;
* subgroepen;
* lerende vs rest;
if beroepssubgroep_num ne . then lerende= (beroepssubgroep_num=2);
label lerende='leerling verpleegkundige';
format lerende lerende.;
* verpleegkundige+verzorgende, helpende+serviceassistent, lerende (overige weg);
if beroepssubgroep_num ne 3 then do;
if beroepssubgroep_num in (5,6) then beroep=1;
else if beroepssubgroep_num in (1,4) then beroep=2;
else if beroepssubgroep_num in (2) then beroep=3;
end;
format beroep beroep.;
label beroep='opleiding (1=vpk/vz,2=hlp/srv,3=ll)';
* type afdeling;
	* "VA1", "VA2", "VA3","VD2", "VE", "VH": revalidation;
if id_afdeling in (1, 2,3,8,9,12) then type_afd=1;
	* "VB","VC1","VF", "VG", "VI3", "VJ1", "VK", "VM1", "VM3": pg;
else if id_afdeling in (4,5,10,11,15, 17, 19,21,23) then type_afd=2;
	* "VC2","VD1", "VI1", "VJ2","VL", "VN1", "VN2", "VN3": som;
else if id_afdeling in (6,7,13,18,20,25,26,27) then type_afd=3;
format type_afd type.;
label type_afd='type afdeling (1=rev,2=pg,3=som)';
run;








***** DATA CHECKS ***********;
title "data checks";

proc freq data=ds; table id_verpleeghuis meting tijdstipinterventie/ missing;run;
proc sort data=ds; by tijdstipinterventie; run;
proc freq data=ds; by tijdstipinterventie; table id_verpleeghuis/missing;run;
proc freq data=ds; by tijdstipinterventie; table id_verpleeghuis*trt*meting/missing nocol;run;

title2 "types verpleegafdelingen";
proc freq data=ds; table id_afdeling*type_afd/missing;run;
title2 "Check maken beroepsgroepen";
proc freq data=ds; table beroepssubgroep_num *lerende/missing; table beroepssubgroep_num*beroep/missing;run;

title2 "check aangemaakte hygiene variabelen: mogelijke patronen";
proc tabulate data=ds; class HH_uitgev_gecorr_asep handenwassen HH_asep_wassen; var id_verpleeghuis;
table HH_uitgev_gecorr_asep*handenwassen*HH_asep_wassen, id_verpleeghuis*(n nmiss);
run; 
proc tabulate data=ds; class HH_uitgev_gecorr_asep handschoenen HH_asep_handschoen; var id_verpleeghuis;
table HH_uitgev_gecorr_asep*handschoenen* HH_asep_handschoen, id_verpleeghuis*(n nmiss);
run; 
proc tabulate data=ds; class HH_uitgev_gecorr_asep handdesinfectie HH_asep_desinfec; var id_verpleeghuis;
table HH_uitgev_gecorr_asep*handdesinfectie* HH_asep_desinfec, id_verpleeghuis*(n nmiss);
run; 
proc tabulate data=ds; class HH_uitgev_gecorr_asep handschoenen handenwassen HH_asep_hschoen_wassen; var id_verpleeghuis;
table HH_uitgev_gecorr_asep*handschoenen*handenwassen*HH_asep_hschoen_wassen, id_verpleeghuis*(n nmiss);
run; 
proc tabulate data=ds; class HH_uitgev_gecorr_asep handschoenen handdesinfectie HH_asep_hschoen_des; var id_verpleeghuis;
table HH_uitgev_gecorr_asep*handschoenen*handdesinfectie*HH_asep_hschoen_des, id_verpleeghuis*(n nmiss);
run; 

*check by moment versions of the variables;
%macro checkmoment(var=);
%DO i=1 %TO 5;
title6 "variable &var._M&i: &var at moment &i";
proc tabulate data=ds missing ; class M&i &var &var._M&i; var id_verpleeghuis;
table M&i*&var*&var._M&i, id_verpleeghuis*(n nmiss );run;
%END;
%MEND;

* hh_uitgev_gecorr_asep
%checkmoment(var=hh_uitgev_gecorr_asep);

* op M2 komt alleen de waarde 0 voor;
%checkmoment(var=hh_asep_wassen);
proc freq data=ds; by M2; where M2=1; table HH_uitgev_gecorr_asep* handenwassen /missing;run; 
/* 
title "VA1 VA2 are wards within home VA etc";
proc freq data=ds; table id_verpleeghuis*id_afdeling /missing;run;
*/




*** STATISTICAL ANALYSIS ***;

**********************************************************;
************* ITT ANALYSES *******************************;
**********************************************************;
ods pdf style=statistical file="ITT_CHANGE.pdf" compress=9;
title1 "ITT analyse (CHANGE trial)";
%let pathname =%sysget(SAS_EXECFILEPATH);
footnote "&pathname";*place and name of program used;

title2 "variables in dataset";
proc contents data=ds;run;

title2 "all data by verpleeghuis";
%rm_design(ds=ds, cluster=id_verpleeghuis);
title2 "all data by unit";
%rm_design(ds=ds, cluster=id_afdeling);

title2 "Across all moments";
%sw_bin(outcome=HH_uitgev_gecorr_asep, dseff=table2, dsdescrip=table2_descrip,period_min=0, period_max=4);
* vraag reviewer over confounding: %sw_bin(outcome=HH_uitgev_gecorr_asep, dseff=table2, dsdescrip=table2_descrip,period_min=0, period_max=4,covars=beroep,covars_cat=beroep);
%sw_bin(outcome=hh_asep_wassen, dseff=table2, dsdescrip=table2_descrip, period_min=0, period_max=4);

%sw_bin(outcome=HH_asep_desinfec, dseff=table2, dsdescrip=table2_descrip, period_min=0, period_max=4);

%sw_bin(outcome=HH_asep_handschoen, dseff=table2, dsdescrip=table2_descrip, period_min=0, period_max=4);

	options orientation=landscape;
	ods rtf style=minimal file="Table2_overallmoments.doc";
	title3 "effects on binary outcomes";
	data table2x; set table2; 
	*add the exponentiated estimates and 95%-CI for log transformed variables;
	if index(outcome_label,'log')>0 then do; 
	e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
	if index(label, 'trt') >0 ; * only the treatment estimates; 
	run;
	proc print data=table2x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; run;

	title3 "descriptives binary outcomes";
	proc tabulate data=table2_descrip order=data; class outcome outcome_label stat trt;var value;
	table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
	ods rtf close;



*** per moment**;
%macro bymoment(var=);
%do i=1 %to 5;
data _moment; set ds; if &var._M&i ne .;run; *only the observations for which this moment is applicable;
%sw_bin(ds=_moment,outcome=&var._M&i, dseff=table3, dsdescrip=table3_descrip, period_min=0, period_max=4);
%end;
%mend bymoment;
title2 "Handhygiene uitgevoerd per moment";
%bymoment(var=HH_uitgevoerd);
title2 "***APART:  hh gecorrigeeerd op moment 2***";
%sw_bin(outcome=HH_uitgev_gecorr_asep_M2, dseff=table2, dsdescrip=table2_descrip,period_min=0, period_max=4);

title2 "Handwassen per moment";
%bymoment(var=HH_asep_wassen);


title2 "Handdesinfectie per moment";
%bymoment(var=HH_asep_desinfec);

title2 "Handschoenen per moment";
%bymoment(var=HH_asep_handschoen);

* ook voor de combinaties;
title2 "Handschoenen en wassen per moment";
%bymoment(var=HH_asep_hschoen_wassen);
title2 "Handschoenen en desinfectie per moment";
%bymoment(var=HH_asep_hschoen_des);

	options orientation=landscape;
	ods rtf style=minimal file="table3_permoment.doc";
	title2 "effects on binary outcomes";
	data table3x; set table3; 
	*add the exponentiated estimates and 95%-CI for log transformed variables;
	if index(outcome_label,'log')>0 then do; 
	e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
	if index(label, 'trt')>0 ; * only the treatment estimates; 
	run;
	proc print data=table3x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; run;

	title2 "descriptives binary outcomes";
	proc tabulate data=table3_descrip order=data; class outcome outcome_label stat trt;var value;
	table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
	ods rtf close;




title2 "descriptief per huis: pooled ctl vs pooled interventie periodes";

%macro ctl_exp(data=ds, by=id_verpleeghuis, trt=trt, outcome=);
* make percentages based on {0,1} coded data;
data &data; set &data;
pct_&outcome=100*&outcome;
%DO i=1 %TO 5; 
pct_&outcome._M&i=100*&outcome._M&i;
%END;
run;
* overall...; 
data _null_; set &data; call symput('varlabel', vlabel(&outcome));run;
title5 "variable (pct_&outcome): &varlabel over all moments";
proc tabulate data=&data;class &by &trt;
var pct_&outcome;
table &by, &trt*(pct_&outcome=' ')*(mean='percent')*f=5.2 &trt*(pct_&outcome=' ')*(N Nmiss)/ rts=10;
run;
* for each moment ...; 
%DO i=1 %TO 5;
data _null_; set &data; call symput('moment',vlabel(M&i)); run;
title5 "variable (pct_&outcome._M&i): &varlabel at moment M&i";
title6 "M&i= &moment";
proc tabulate data=&data;class &by &trt;
var pct_&outcome._M&i;
table &by, &trt*(pct_&outcome._M&i=' ')*(mean='percent')*f=5.2 &trt*(pct_&outcome._M&i=' ')*(N Nmiss)/ rts=10;
run;
%END;
%mend;

	options orientation=landscape;
	ods rtf style=minimal file="table descriptief per huis ctl vs interventie.doc";
title2 "Handhygiene uitgevoerd gecorr asep: controle vs interventie periode per huis";
%ctl_exp(outcome=HH_uitgev_gecorr_asep);
/* check
proc freq data=ds; table &by*HH_uitgevoerd_M1*&trt /missing;run; 
*/
title2 "Handhygiene uitgevoerd : controle vs interventie periode per huis";
%ctl_exp(outcome=HH_uitgevoerd);

title2 "Handhygiene uitgevoerd asep + wassen : controle vs interventie periode per huis";
%ctl_exp(outcome=HH_asep_wassen);
title2 "Handhygiene uitgevoerd asep + desinfectie : controle vs interventie periode per huis";
%ctl_exp(outcome=HH_asep_desinfec);
title2 "Handhygiene uitgevoerd asep + handschoenen: controle vs interventie periode per huis";
%ctl_exp(outcome=HH_asep_handschoen);
	ods rtf close;



title2 "subgroepen";

title2 "subgroep beroep: verpl/verz, help/serv, lerende (overige niet meegenomen)";
* hh asep gecorrigeerd;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=HH_uitgev_gecorr_asep,
interaction=beroep, 
estimate_statement=%str(estimate '1:vpk/vz' trt 1 trt*beroep 0 0 1/cl; estimate '2: hlp/srv' trt 1 trt*beroep 1 0 0/cl ; 
estimate '3:ll' trt 1 trt*beroep 0 1 0/cl; ), contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep wassen;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=hh_asep_wassen,
interaction=beroep, 
estimate_statement=%str(estimate '1:vpk/vz' trt 1 trt*beroep 0 0 1/cl; estimate '2: hlp/srv' trt 1 trt*beroep 1 0 0/cl ; 
estimate '3:ll' trt 1 trt*beroep 0 1 0/cl; ), contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep desinfectie;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=HH_asep_desinfec,
interaction=beroep, 
estimate_statement=%str(estimate '1:vpk/vz' trt 1 trt*beroep 0 0 1/cl; estimate '2: hlp/srv' trt 1 trt*beroep 1 0 0/cl ; 
estimate '3:ll' trt 1 trt*beroep 0 1 0/cl; ), contrast_statement=, dsdescrip=table4_descrip, dseff=table4);


title2 "subgroep: lerende vs niet-lerende";
%let outcome=HH_uitgev_gecorr_asep;
%sw_subgr_bin(ds=ds, cat_min=0, cat_max=1, outcome=&outcome,
interaction=lerende, 
estimate_statement=%str(estimate 'lerende=0=nee' trt 1 trt*lerende 0 1/cl; estimate 'lerende=1=ja' trt 1 trt*lerende 1 0/cl;), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep wassen;
%let outcome=hh_asep_wassen;
%sw_subgr_bin(ds=ds, cat_min=0, cat_max=1, outcome=&outcome,
interaction=lerende, 
estimate_statement=%str(estimate 'lerende=0=nee' trt 1 trt*lerende 0 1/cl; estimate 'lerende=1=ja' trt 1 trt*lerende 1 0/cl;), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep desinfectie;
%let outcome=HH_asep_desinfec;
%sw_subgr_bin(ds=ds, cat_min=0, cat_max=1, outcome=&outcome,
interaction=lerende, 
estimate_statement=%str(estimate 'lerende=0=nee' trt 1 trt*lerende 0 1/cl; estimate 'lerende=1=ja' trt 1 trt*lerende 1 0/cl;), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);

title2 "subgroep type afdeling: revalidatie/pg/somatisch"; 
* asep;
%let outcome=HH_uitgev_gecorr_asep;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=&outcome,
interaction=type_afd, 
estimate_statement=%str(estimate 'revalidatie' trt 1 trt*type_afd 0 1 0 / cl; estimate 'pg' trt 1 trt*type_afd 1 0 0 / cl; 
estimate 'somatisch' trt 1 trt*type_afd 0 0 1/cl; ), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep wassen;
%let outcome=hh_asep_wassen;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=&outcome,
interaction=type_afd, 
estimate_statement=%str(estimate 'revalidatie' trt 1 trt*type_afd 0 1 0 / cl; estimate 'pg' trt 1 trt*type_afd 1 0 0 / cl; 
estimate 'somatisch' trt 1 trt*type_afd 0 0 1/cl; ), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);
* asep desinfectie;
%let outcome=HH_asep_desinfec;
%sw_subgr_bin(ds=ds, cat_min=1, cat_max=3, outcome=&outcome,
interaction=type_afd, 
estimate_statement=%str(estimate 'revalidatie' trt 1 trt*type_afd 0 1 0 / cl; estimate 'pg' trt 1 trt*type_afd 1 0 0 / cl; 
estimate 'somatisch' trt 1 trt*type_afd 0 0 1/cl; ), 
contrast_statement=, dsdescrip=table4_descrip, dseff=table4);


	options orientation=landscape; 
	ods rtf style=minimal file="Table4_Subgroups.doc";
	title "Subgroups: Efficacy ";
	data table4x; set table4; 
	*add the exponentiated estimates and 95%-CI for log transformed variables;
	if index(outcome_label,'log')>0 then do; 
	e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
	if index(label, 'trt')=0; * not the rows with 'trt' in it; 
	run;
	proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; run;

	title "Subgroups: descriptives ";
	proc tabulate data=table4_descrip order=data; class outcome outcome_label stat trt;var value;
	table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
	ods rtf close;


title1 "reden observeren van invloed op trial resultaten?";
proc freq data=ds; table RedenObserverenbekend/missing;run;
title2 " nee -> weten=0, ja -> weten=1, overig -> weten=.";
data ds1; set ds;
if compress(RedenObserverenbekend)="nee" then weten=0;
else if compress(RedenObserverenbekend)="ja" then weten=1;
else weten=.;
run;

proc freq data=ds1; table weten;run;
title2 "verbetering in effect door weten geobserveerd te zijn";
proc mixed data=ds1; class cluster subcluster period;
model HH_uitgev_gecorr_asep= trt period weten weten*trt / solution;
random intercept / subject=cluster;
random intercept / subject=subcluster(cluster);
estimate 'effect in weten=0' trt 1 /cl;
estimate 'effect in weten=1' trt 1 trt*weten 1 /cl;
run;

ods pdf close;


ods rtf file="verzoek_Anja_20200403.doc";

data ds1; set ds; pct_HH_uitgev_gecorr_asep=100*HH_uitgev_gecorr_asep;run;
proc tabulate data=ds1;
class id_afdeling trt;
var pct_HH_uitgev_gecorr_asep;
table id_afdeling, trt*pct_HH_uitgev_gecorr_asep*(mean n nmiss);run;
ods rtf close;
