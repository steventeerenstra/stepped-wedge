* remove: period_min and period_max in the infer macro's: determine them from the dataset?;
* add comment variable  in the output of the proc mixed/glimmix, so that PP/ITT, adjusted etc 
      can be put in the outcome table; 
* add option to turn off residual-pred and obs-pred plots at individual level or cluster level; 
* ordering in ITS plots on order when intervention and a reference line for that?;
* do the data handling to adapt the estimates table for making a forest much more ;
* add a use of fit-method in proc glimmix part;
*     separately from the macro that gives the estimates;
/* to automate writing contrast and estimate statements: use x=" "; xx = cat(repeat(' 0',5), x);
to make a series of _0_0_0 etc; and then substr(xx,4, 1)=1; to replace a 0 into a 1
*/

**** MACRO DEFINITIONS ****;

%macro rm_clusperiod(ds=ds, cluster=center,period=period);
* rm= repeated measures;
* which clusters and periods present?;
title6 "clusters in the design";
proc freq data=&ds; table &cluster/ missing nocum ; run;
title6 "periods in the design";
proc tabulate data=&ds; class &period; table &period*N;run;
%mend;

%macro rm_design(ds=ds, cluster=center,trt=trt, period=period);
* general check on repeated measures design;
title6 "treatments applied per cluster-period";
title6 "check whether each cluster-period has one treatment?";
proc tabulate data=&ds;class &cluster &period &trt;format &period 4.0; 
table &cluster*&trt, &period*N*f=4. /rts=15 ;
run;
%mend;


* continuous outcome";
* descriptives;
%macro sw_cont_descrip(ds=ds,outcome=, cluster=center, period=period, trt=trt,
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
title5 "shift in distribution by intervention: over all clusters";
proc sgplot data=_shift;
histogram intervention_&outcome /transparency=0.3; *fillattrs=(color=lightgreen);
histogram control_&outcome /transparency=0.6; *fillattrs=(color=lightred);
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

%macro sw_cont_infer(ds=ds,outcome=, cluster=center, period=period, trt=trt,
period_min=1,period_max=14,covars_cat=, covars=, interaction=,
dseff=, random_statement=%str(random intercept /subject=&cluster; ), ddfm=satterthwaite,
repeated_statement=%str( ), estimate_statement=%str( ), contrast_statement=%str( ), alpha=0.05
);

* to do: ITS analyse by cluster;
* to do: Hooper-Girling, exponential decay model;
* to do: check and describe which variables are in dseff: 
        effect, outcome, outcome_label, effect, model, estimate, stderr, lower, upper;
* to do: different size for observed vs predicted? GTL scatter plot has an option for that;
* to do: reorder flow of calculation: if non-converged set to . everything, if not then fill in estimates;

*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;

* get for this outcome (for calculation Searle CI for the ICC DOI: 10.1002/sim.1330);
* average cluster size, weighted average cluster size, and total number of non-missing observations;
ods exclude all;
ods output table=_clustersizes;
proc tabulate data=&ds ;class &cluster; var &outcome; table &cluster*&outcome*(N) ;run; 
ods output close;
ods exclude none;
data _clustersizes; set _clustersizes end=last; retain n_clus 0 sum_terms 0 sum_squares 0;
n_clus+1; sum_terms=sum_terms+&outcome._N; sum_squares=sum_squares+&outcome._N**2;
if last then do; n0=(1/(n_clus-1))*( sum_terms - (1/sum_terms)*sum_squares);
	call symputx('clussize_0', n0,'L');* L for local, G for global macro variable;
	call symputx('nclus', n_clus,'L');
	call symputx('nobs', sum_terms, 'L');* total number of observations;
	call symputx('clussize_av',sum_terms/n_clus, 'L');
	end;
run;

title3 "analysis   **** &outcome: &outcome_label ****";
title4 "modeling categorical time trend, random effect cluster (Hussey & Hughes)";
* if dseff asked then save effects for &trt and possibly estimate from estimate statement; 
%IF &dseff ne %THEN 
%DO;
ods output ConvergenceStatus=_converged;
ods output SolutionF=_solutions;
ods output CovParms=_covparms;
	%IF &estimate_statement ne %THEN ods output Estimates=_estimates;;	
	%IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
	* if interaction asked then also make interaction term;
	%local interaction_term; %let interaction_term=;
	%IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;
****** Hussey and Hughes model ****;
proc mixed data=&ds;
class &cluster &period &covars_cat ;
model &outcome= &trt &period &covars &interaction_term/ solution cl outpred=_respred ddfm=&ddfm; 
&random_statement;;
&repeated_statement;;
&estimate_statement;;
&contrast_statement;;
run;
* process the effects;
%IF &dseff ne %THEN
%DO; 
ods output close;
* depending on convergence status, status=0 means convergenced, estimates etc are set missing;
data _null_; set _converged; call symput('nonconverged', status);run;

	* get trt solution in one row;
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
	* calculate and save intracluster correlations (assuming Hussey &Hussey model);
	data _covparms(keep=ICC ICC_lower ICC_upper var_cluster var_resid); set _covparms; 
	retain var_cluster var_resid;
	if covparm="Intercept" then var_cluster=estimate;
	if covparm="Residual" then do;*this is the last record;
        var_resid=estimate;
		ICC=var_cluster/(var_resid+var_cluster); 
		* confidence interval using Searle method, doi: 10.1002/sim.1330;
		* alternative would be Fisher transformation;
		F_U=quantile('F', 1-&alpha/2,&nclus-1,&nobs-&nclus); F_L=quantile('F',&alpha/2,&nclus-1,&nobs-&nclus);
		F=1+&clussize_0*var_cluster/var_resid;
		ICC_lower=(F/F_U -1)/(&clussize_0 + F/F_U -1);
			ICC_lower=max(ICC_lower,0); *negative ICCs are truncated at 0;
		ICC_upper=(F/F_L-1)/(&clussize_0 + F/F_L -1);
		*ICC_F= (F-1)/(F+&clussize_0 -1);* this is how ICC is calculated from F-value;
		output;
		end;
	run;
	data _solutions_icc; merge _solutions _covparms;run;

		%IF &nonconverged ne 0 %THEN %DO; 
		data _solutions_icc;set _solutions_icc; estimate=.;stderr=.;probt=.;lower=.;upper=.;icc=.;icc_lower=.; icc_upper=.;var_cluster=.; var_resid=.; run;
		%END;
	proc append force base=&dseff data=_solutions_icc ;run;
%IF &estimate_statement ne %THEN %DO;
	data _estimates; 
	length label $50;length outcome $30;length outcome_label $100; length model $100;
	set _estimates; 
	outcome="&outcome";outcome_label="&outcome_label";
    model="&trt &period &covars &interaction_term // &random_statement &repeated_statement";
	keep outcome outcome_label label estimate stderr probt lower upper model; 
	run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _estimates;set _estimates; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
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
		%IF &nonconverged ne 0 %THEN %DO; 
		data _contrasts;set _contrasts; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
	proc append force base=&dseff data=_contrasts;run;
	%END;
%END;
*** end Hussey & Hughes model**;

title6 "residual vs predicted plot at subject level";
proc sgplot data=_respred; scatter x=pred y=resid ; refline 0/axis=y;run;
title6 "residual vs predicted plot at subject level: by cluster";
proc sgpanel data=_respred; panelby &cluster /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
title6 "residual vs predicted plot of at subject level: by period";
proc sgpanel data=_respred; panelby &period /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
/* too much detail
title6 "residual plot: by cluster-period";
proc sgpanel data=_respred; panelby &cluster &period ; scatter x=pred y=resid; refline 0 / axis=y;run;
*/
* fit: residuals and observed values at cluster level; 
proc means data=_respred noprint; class &cluster &period; id &trt;
output out=_pred1 mean(pred &outcome )=mean_pred mean_obs ;run;
data _pred1; set _pred1; mean_resid=mean_obs-mean_pred;run;
title6 "observed and predicted plot >averaged at cluster level<: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
series x=&period y=mean_pred;
scatter x=&period y=mean_obs;
rowaxis integer;
run;
title6 "residual vs predicted plot >averaged at cluster level<";
proc sgplot data=_pred1; 
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot >averaged at cluster level<: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot >averaged at cluster level<: by period";
proc sgpanel data=_pred1; panelby &period /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
%mend sw_cont_infer;


%macro sw_cont(ds=ds,outcome=, cluster=center, period=period, trt=trt,
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
%sw_cont_infer(ds=&ds,outcome=&outcome, cluster=&cluster, period=&period, trt=&trt,
period_min=&period_min,period_max=&period_max,covars_cat=&covars_cat, covars=&covars, interaction=&interaction,
dseff=&dseff, random_statement=&random_statement,
repeated_statement=&repeated_statement, estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%mend sw_cont;

******; 
******; 
******; 
%macro sw_bin_descrip(ds=ds, outcome=,cluster=center, period=period, trt=trt,period_min=1, period_max=14,
dsdescrip=
); 
** assumes variable is {0,1} coded;

data &ds; set &ds; 
percent_&outcome=100*&outcome;* for calculating percentages via a means-operator;
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
vbar dummy /  stat=mean response=percent_&outcome;
colaxis label=' ';
run;
title4 "data size per cluster-period";
proc tabulate data=&ds; class &cluster &period; var &outcome;
table &cluster, &period="&outcome: N per cluster-period"*&outcome=' '*(N=' ')*f=5.0 / rts=10;  run;
title4 "percent per cluster-period";
proc tabulate data=&ds; class &cluster &period; var percent_&outcome;
table &cluster, &period="&outcome: percent per cluster-period"*percent_&outcome=' '*(mean=' ')*f=5.2 / rts=10;run;


title3 "analysis   **** &outcome: &outcome_label ****";
title4 "time series plot per cluster";
proc means data=&ds noprint; id &trt; class &cluster &period;
var percent_&outcome; output out=_ClusPeriodMeans mean(percent_&outcome)=percent;run;
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
vbar &trt / response=percent_&outcome stat=mean;
run;
* show descriptives and save if asked ;
%IF &dsdescrip ne %THEN ods output CrossTabFreqs=_descrip;; 
proc freq data=&ds; table &outcome*&trt / missing norow nopercent; run;
%IF &dsdescrip ne %THEN %DO;
ods output close;
data _descrip; length stat $ 10;length outcome $30; length outcome_label $ 100; 
set _descrip; 
outcome="&outcome";outcome_label="&outcome_label";
if &trt ne . and &outcome=. then do;
stat="ntot"; value=frequency; output;
end;
if &outcome=1 and &trt ne . then do; 
stat="n"; value=frequency; output;
stat="perc"; value=colpercent; output;
end;
keep outcome outcome_label &trt stat value;
run;
proc append force base=&dsdescrip data=_descrip;run;
%END;
%mend sw_bin_descrip;

%macro sw_bin_infer(ds=ds, outcome=,cluster=center, period=period, trt=trt,period_min=1, period_max=14,
random_statement=%str(random intercept /subject=&cluster; ),ddfm=satterthwaite, repeated_statement=%str( ),
estimate_statement=%str(estimate 'ctl' intercept 1 period  0 0 0 0 0 0 0 0 0 0 0 0 0 1 / cl; estimate 'exp' intercept 1 trt 1 period  0 0 0 0 0 0 0 0 0 0 0 0 0 1 / cl),
contrast_statement=%str( ),dseff=,covars_cat=, covars=, interaction=
); 
** assumes variable is {0,1} coded;

data &ds; set &ds; 
percent_&outcome=100*&outcome;* for calculating percentages via a means-operator;
dummy=1;* for having a categorical variable with only one value;
run;
*get outcome label;
data _null_; set &ds; call symput('outcome_label', vlabel(&outcome));run;

title3 "analysis   **** &outcome: &outcome_label ****";
title4 "modeling categorical time trend, random effect cluster (Hussey & Hughes)";
title5 "** % difference ** (linear mixed model)";
* if dseff asked then save effects for &trt and possibly estimate from estimate statement; 
%IF &dseff ne %THEN 
%DO;
ods output SolutionF=_solutions;
ods output convergencestatus=_converged;
	%IF &estimate_statement ne %THEN ods output Estimates=_estimates;;	
	%IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
	* if interaction asked then also make interaction term;
	%local interaction_term; %let interaction_term=;
	%IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;
****** % difference: Hussey and Hughes model ****;
proc mixed data=&ds;
class &cluster &period &covars_cat ;
model percent_&outcome= &trt &period &covars &interaction_term/ solution cl outpred=_respred ddfm=&ddfm ; 
&random_statement;;
&repeated_statement;;
&estimate_statement;;
&contrast_statement;;
run;
* process the effects;
%IF &dseff ne %THEN
%DO; 
ods output close;
data _null_; set _converged; call symput('nonconverged', status);run;

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
	* calculate and save intracluster correlations (assuming Hussey &Hussey model);
	data _covparms(keep=ICC ICC_lower ICC_upper var_cluster var_resid); set _covparms; retain var_cluster var_resid;
	if covparm="Intercept" then var_cluster=estimate;
	if covparm="Residual" then do;*this is the last record;
        var_resid=estimate;
		ICC=var_cluster/(var_resid+var_cluster); 
		* confidence interval using Searle method, doi: 10.1002/sim.1330;
		* alternative would be Fisher transformation;
		F_U=quantile('F', 1-&alpha/2,&nclus-1,&nobs-&nclus); F_L=quantile('F',&alpha/2,&nclus-1,&nobs-&nclus);
		F=1+&clussize_0*var_cluster/var_resid;
		ICC_lower=(F/F_U -1)/(&clussize_0 + F/F_U -1);
			ICC_lower=max(ICC_lower,0); *negative ICCs are truncated at 0;
		ICC_upper=(F/F_L-1)/(&clussize_0 + F/F_L -1);
		*ICC_F= (F-1)/(F+&clussize_0 -1);* this is how ICC is calculated from F-value;
		output;
		end;
	run;
	data _solutions_icc; merge _solutions _covparms;run;

		%IF &nonconverged ne 0 %THEN %DO; 
		data _solutions_icc;set _solutions_icc; estimate=.;stderr=.;probt=.;lower=.;upper=.;icc=.;icc_lower=.; icc_upper=.;var_cluster=.; var_resid=.; run;
		%END;
	proc append force base=&dseff data=_solutions_icc ;run;
%IF &estimate_statement ne %THEN %DO;
	data _estimates; 
	length label $50;length outcome $30;length outcome_label $100; length model $100;
	set _estimates; 
	outcome="&outcome";outcome_label="&outcome_label";
    model="LMM: &trt &period &covars &interaction_term | &random_statement &repeated_statement";
	keep outcome outcome_label label estimate stderr probt lower upper model; 
	run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _estimates;set _estimates; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
	proc append force base=&dseff data=_estimates;run;
	%END;
%IF &contrast_statement ne %THEN %DO;
	data _contrasts;
	length label $50;length outcome $30;length outcome_label $100; length model $100;
	set _contrasts;
	outcome="&outcome";outcome_label="&outcome_label";
    model="LMM: &trt &period &covars &interaction_term // &random_statement &repeated_statement";
	probt=probf; estimate=.;stderr=.;lower=.;upper=.;
	keep outcome outcome_label label estimate stderr probt lower upper model; 
	run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _contrasts; set _contrasts; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
	proc append force base=&dseff data=_contrasts;run;
	%END;	
%END;
*** end % difference: Hussey & Hughes model**;
title6 "residual vs predicted plot at subject level";
proc sgplot data=_respred; scatter x=pred y=resid ; refline 0/axis=y;run;
title6 "residual vs predicted plot at subject level: by cluster";
proc sgpanel data=_respred; panelby &cluster /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
title6 "residual vs predicted plot of at subject level: by period";
proc sgpanel data=_respred; panelby &period /onepanel; scatter x=pred y=resid; refline 0 / axis=y;run;
/* too much detail
title6 "residual plot: by cluster-period";
proc sgpanel data=_respred; panelby &cluster &period ; scatter x=pred y=resid; refline 0 / axis=y;run;
*/
* fit: residuals and observed values at cluster level; 
proc means data=_respred noprint; class &cluster &period; id &trt;
output out=_pred1 mean(pred &outcome )=mean_pred mean_obs ;run;
data _pred1; set _pred1; mean_resid=mean_obs-mean_pred;run;
title6 "observed and predicted plot >averaged at cluster level<: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
series x=&period y=mean_pred;
scatter x=&period y=mean_obs;
rowaxis integer;
run;
title6 "residual vs predicted plot >averaged at cluster level<";
proc sgplot data=_pred1; 
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot >averaged at cluster level<: by cluster";
proc sgpanel data=_pred1; panelby &cluster /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;
title6 "residual vs predicted plot >averaged at cluster level<: by period";
proc sgpanel data=_pred1; panelby &period /onepanel;
scatter x=mean_pred y=mean_resid; refline 0 / axis=y;
run;


title5 "** odds ratio ** (generalized linear mixed model)";
%IF &dseff ne %THEN 
%DO;
ods output ParameterEstimates=_solutions;
ods output Covparms=_covparms;
ods output ConvergenceStatus=_converged;
	%IF &estimate_statement ne %THEN ods output Estimates=_estimates;;	
	%IF &contrast_statement ne %THEN ods output Contrasts=_contrasts;;
%END;
	* if interaction asked then also make interaction term;
	%local interaction_term; %let interaction_term=;
	%IF &interaction ne %THEN %DO;%let interaction_term=&trt*&interaction; %END;
*model;
proc glimmix data=&ds; 
class &cluster &period &covars_cat;
model &outcome(reference='0')=&trt  &period &covars &interaction_term
/solution cl distribution=binary oddsratio; *&period &covars oddsratio;
&random_statement;;
&estimate_statement;;
&contrast_statement;;
output out=_respred pred=lp;
run;
* process the effects;
%IF &dseff ne %THEN
%DO; 
ods output close;
* depending on convergence status, status=0 means convergenced, estimates etc are set missing;
data _null_; set _converged; call symput('nonconverged', status);run;

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
	* calculate and save intracluster correlations (assuming Hussey &Hussey model);
	* icc on logistic scale;
	/*  this code has still to be augmented and checked;
	    data _covparms(keep=ICC cov_cluster); set _covparms; 
		if parameter="&cluster" then do; 
		cov_cluster=estimate;icc=cov_cluster/(cov_cluster + constant("Pi")*constant("Pi")/3); 
		output;end;
		run;
	data _solutions_icc; merge _solutions _covparms;run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _solutions;set _solutions; estimate=.;stderr=.;probt=.;lower=.;upper=.; icc=.; cov_cluster=.;run;
		%END;
	proc append force base=&dseff data=_solutions_icc ;run;
	*/
	proc appende force base=&dseff data=_solutions;run;
%IF &estimate_statement ne %THEN %DO;
	data _estimates; 
	length label $50;length outcome $30;length outcome_label $100; length model $100;
	set _estimates; 
	outcome="&outcome";outcome_label="log(odds) &outcome_label";
    model="GLMM: &trt &period &covars &interaction_term | &random_statement ";
	keep outcome outcome_label label estimate stderr probt lower upper model; 
	run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _estimates;set _estimates; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
	proc append force base=&dseff data=_estimates;run;
	%END;	
%IF &contrast_statement ne %THEN %DO;
	data _contrasts;
	length label $50;length outcome $30;length outcome_label $100; length model $100;
	set _contrasts;
	outcome="&outcome";outcome_label="log(odds) &outcome_label";
    model="GLMM: &trt &period &covars &interaction_term // &random_statement &repeated_statement";
	probt=probf; estimate=.;stderr=.;lower=.;upper=.;
	keep outcome outcome_label label estimate stderr probt lower upper model; 
	run;
		%IF &nonconverged ne 0 %THEN %DO; 
		data _contrasts;set _contrasts; estimate=.;stderr=.;probt=.;lower=.;upper=.;run;
		%END;
	proc append force base=&dseff data=_contrasts;run;
	%END;
%END;
*** end odds ratio: Hussey & Hughes model**;
title6 "observed vs predicted plot: by cluster";
;
* calculate residuals, observed and fitted for cluster x time averages;
data _respred; set _respred; pred_percent=100*exp(lp)/(1+exp(lp)); run;
proc means data=_respred noprint; class &cluster &period;
output out=_pred1 mean(pred_percent percent_&outcome)=mean_pred mean_obs;run;
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
%mend sw_bin_infer;

%macro sw_bin(ds=ds, outcome=,cluster=center, period=period, trt=trt,period_min=1, period_max=14,
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
%sw_bin_infer(ds=&ds, outcome=&outcome,cluster=&cluster, period=&period, trt=&trt,period_min=&period_min, 
period_max=&period_max, random_statement=&random_statement,repeated_statement=&repeated_statement,
estimate_statement=&estimate_statement,contrast_statement=&contrast_statement,
dseff=&dseff,covars_cat=&covars_cat, covars=&covars, 
interaction=&interaction
); 
%mend sw_bin;
*options nosymbolgen nomprint nomlogic;
*options symbolgen mprint mlogic;

********* READ 	DATA *******;
options pagesize=60;
libname data ".\data";
* get formats;
options fmtsearch= (data.formats);
* show variables in dataset;
proc contents data=data.effampart;run;

* show coding of selected formats;
*libname library 'SAS-library';
proc format library=data.formats fmtlib;
   select control pt_sex;* name of the format without the .;
run;

/*** needed:
1.	Tevredenheid (CQI_general_satisfaction)
2.	Angst (HADS_ANXIETY)
3.	Depressie (HADS_DEP)
4.	PTSS (IES_R_MEAN)
***/

* process dataset;
data ds0; set data.effampart;
trt=control_intervention; *0=control, 1=intervention;
tevredenheid=cqi_general_satisfaction; label tevredenheid="tevredenheid (CQI)";
angst=hads_anxiety;label angst="angst (HADS)";
depressie=hads_dep; label depressie="depressie (HADS)";
ptss=ies_r_mean; label ptss="PTSS (IES-R)";
log_tevr_revers=log(11-tevredenheid); label log_tevr_revers="tevredenheid (log[11- ])";
log_ptss=log(ptss+1); label log_ptss="PTSS (log+1)";
log_depressie=log(depressie+1); ;label log_depressie="depressie (log+1)";
log_angst=log(angst); label log_angst="angst (log)";
* covariaat die scheef is;
log_los=log(pt_los_icu);
run;

* final ITT dataset;
data itt; set ds0; run;
 
data pp; set ds0; where subanalyse_pp=1;run;


***************************************************************************************************;
%let date=230206;

ods pdf style=statistical file="EFFAMPART_st&date._ITT.pdf" compress=9; * style=journal kan ook;

***** program used  **************;
title " ";footnote "";
%let pathname =%sysget(SAS_EXECFILEPATH);
proc odstext;
p "The program used for these analysis is: &pathname";*place and name of program used;
run;
ods startpage=now;

title1 "EFFAMPART";
***** DATA CHECKS ***********;
title2 "hoe ziet design eruit?";
footnote1 "Isala en ETZ hebben de meeste metingen; periods >= 15 hebben weinig metingen";
%rm_clusperiod(ds=ds0,cluster=center,period=period);


title3 "aantallen toegewezen behandelingen in elke center-period (inclusief trainingsperiodes)";
footnote1 "metingen in beide condities: MMC in periode 6,Bernhoven in periode 11"; 
footnote2 "alleen Isala heeft vrijwel alle periodes metingen";
footnote3 "ETZ heeft alleen controle metingen";
%rm_design(ds=ds0, cluster=center, period=period,trt=trt);


title2 "hoe ziet verdeling uitkomst variabelen eruit";
footnote "relatie tussen ptss & depressi, ptts & angst, ptts & depressies rechts scheve verdelingen behalve tevredenheid";
proc sgscatter data=ds0;
matrix tevredenheid angst depressie ptss
       / diagonal=(histogram);
run;

title2 "verdeling tevredenheid angst depressie ptss over alle metingen heen"; 
footnote1 "Minimale en maximale waarden";
footnote2 "Tevredenheid: 0 (Zeer ontevreden -10 Zeer tevreden)";
footnote3 "HADS (Angst/Depressie): 0 – 21; 8 of hoger is indicatief voor angst/depressie";
footnote4 "IESR: 0-88, van de gemiddelde score is 1.6 of hoger indicatief voor ptss";
proc means data=ds0 n nmiss min q1 median q3 max mean std; 
var tevredenheid angst depressie ptss;run;


footnote "getransformeerde versies meer symmetrisch verdeeld";
footnote1 "nb: log tevredenheid is andersom gecodeerd (hoger is slechter)";
footnote2 "maar de verdeling wordt er niet symmetrischer van";
proc sgscatter data=ds0;
matrix log_tevr_revers log_angst log_depressie log_ptss
       / diagonal=(histogram);
run;



*** STATISTICAL ANALYSIS ***;

title " ";footnote "";

**********************************************************;
************* ITT ANALYSES *******************************;
**********************************************************;
* shortcut name for the dataset;
data ds; set itt;run;
%let table_eff=eff_ITT;
%let table_descrip=descrip_ITT_cont;


******* missing value analysis ***************************;
title2 "Missing Data analysis: welke combinaties van missingness komen voor?";
footnote "missing data vooral op ptts+angst+depressie terwijl tevredenheid wel bekend";
ods select MissPattern;
proc mi data=ds nimpute=0;
var tevredenheid angst depressie ptss;
run;


title2 "ITT analyses per variable";
footnote "";
* tevredenheid;
proc odstext; 
p "tevredenheid op de oorspronkelijke schaal: " /style={textindent=100}; 
p "ondanks scheve verdeling is de fit redelijk (en niet beter op loggetransformeerde schaal)";
p "niet-relevante verslechtering -0.6 maar wel p< 0.05";
p "fit op logschaal is niet wezenlijk beter, dus nemen we de getallen op oorspronkelijke schaal";
run;
ods startpage=no;
%sw_cont(outcome=tevredenheid,outcome_min=0,outcome_max=15, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
%sw_cont(outcome=log_tevr_revers,outcome_min=0,outcome_max=15, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

*angst;
proc odstext;
p "angst : "/style ={textindent=100};
p "goede fit zowel op oorspronkelijk als op log-schaal";
run;
ods startpage=no;
%sw_cont(outcome=angst, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
%sw_cont(outcome=log_angst, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

* depressie;
proc odstext;
p "depressie: "/style={textindent=100};
p "goede fit op logschaal, lijkt ook beter dan op oorspronkelijke schaal (residuen beter rondom 0 gespreid, vooral voor MMC";
run;
ods startpage=no;
%sw_cont(outcome=depressie, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
%sw_cont(outcome=log_depressie, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;


* ptss;
proc odstext;
p "ptss: "/style={textindent=100};
p "goede fit op logschaal, lijkt wat beter dan op oorspronkelijke schaal (residuen beter rondom 0 gespreid, vooral voor de spreiding omhoog en omlaag";
run;
ods startpage=no;
%sw_cont(outcome=ptss, outcome_min=0, outcome_max=70, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
%sw_cont(outcome=log_ptss, outcome_min=0, outcome_max=70, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;


*** summary ********;
title2 "ITT effects for *tevredenheid* in termen van *verschil* intervention minus control";
footnote "";
data table_itt_eff;set eff_ITT; where outcome="tevredenheid";run;
proc print data=table_itt_eff noobs; var outcome_label estimate lower upper probt icc icc_lower icc_upper;run;


title2 "ITT effects in termen of ratio of intervention divided by control";
footnote "ICC is ook op log schaal";
data table_itt_eff;set eff_ITT; where find(outcome,"log","i")> 0 and find(outcome,"log_tevr_revers","i") <= 0;
ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);run;
proc print data=table_itt_eff noobs; var outcome_label ratio ratio_ll ratio_ul probt icc icc_lower icc_upper;run;
footnote "";

ods pdf close;




**************************************************************;
*********** ITT met covariaat correctie Pt_LOS_ICU *******************;
**************************************************************;
data ds; set itt;run;
%let table_eff=eff_ITT_LOScorr;
%let table_descrip=descrip_ITT_LOScorr;


ods pdf style=statistical file="EFFAMPART_st&date._ITT_adj.pdf" compress=9; * style=journal kan ook;
title1 "EFFAMPART";

title2 "verdeling LoS in ITT";
proc means data=ds nmiss min max mean median ; var Pt_LOS_ICU;run;
proc sgplot data=ds; histogram pt_los_icu;run;
title2 "verdeling LoS loggetransformeerd";
proc sgplot data=ds; histogram log_los;run;

title2 "ITT analyses per variable, correctie voor log(LoS)";
proc odstext; 
p "gezien dat de natuurlijke log(LoS) normaal verdeeld is";
p "zullen we voor de loggetransformeerde LoS corrigeren";run;
ods startpage=yes;

proc odstext;
p "Tevredenheid met correctie: niet-relevante verslechtering (p<0.05)";
run;
ods startpage=no;
%sw_cont(outcome=tevredenheid,covars=log_los,outcome_min=0,outcome_max=15, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

proc odstext;
p "angst, depressie en ptss op logschaal, omdat de primaire analyse op logschaal is";
run;
*angst (log);
%sw_cont(outcome=log_angst, covars=log_los,outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

* depressie (log);
%sw_cont(outcome=log_depressie, covars=log_los,outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

* ptss (log);
%sw_cont(outcome=log_ptss, covars=log_los,outcome_min=0, outcome_max=70, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;


*** summary ***;

title2 "ITT effects for *tevredenheid* in termen van *verschil* intervention minus control";
footnote " ";
data table_itt_eff;set eff_ITT_LOScorr; where outcome="tevredenheid";run;
proc print data=table_itt_eff noobs; var outcome_label estimate lower upper probt icc icc_lower icc_upper;run;


title2 "ITT effects in termen of ratio of intervention divided by control";
footnote "ICC is ook op log schaal";
data table_itt_eff;set eff_ITT_LOScorr; where find(outcome,"log","i")> 0 and find(outcome,"log_tevr_revers","i") <= 0;
ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);run;
proc print data=table_itt_eff noobs; var outcome_label ratio ratio_ll ratio_ul probt icc icc_lower icc_upper;run;
footnote " ";



ods pdf close;



**********************************************************;
************* PP ANALYSES *******************************;
**********************************************************;
* shortcut name for the dataset;
data ds; set pp;run;
%let table_eff=eff_pp;
%let table_descrip=descrip_pp_cont;

ods pdf style=statistical file="EFFAMPART_st&date._PP.pdf" compress=9; * style=journal kan ook;
******* missing value analysis ***************************;
title2 "PP: Missing Data analysis: welke combinaties van missingness komen voor?";
footnote "missing data vooral op ptts+angst+depressie terwijl tevredenheid wel bekend";
ods select MissPattern;
proc mi data=ds nimpute=0;
var tevredenheid angst depressie ptss;
run;


title2 "PP analyses per variable";
footnote " ";


proc odstext;
p "in de PP analyse (zonder correctie voor LoS) is het effect voor tevredenheid hetzelfde";
run;
ods startpage=no;
* tevredenheid;
%sw_cont(outcome=tevredenheid,outcome_min=0,outcome_max=15, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

proc odstext; 
p "omdat de ITT analyse voor angst, depressie en ptts op de loggetransformeerde data is, is dat in de PP dat ook: " /style={textindent=100}; 
run;

*angst;
%sw_cont(outcome=log_angst, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

* depressie;
%sw_cont(outcome=log_depressie, outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;


* ptss;
%sw_cont(outcome=log_ptss, outcome_min=0, outcome_max=70, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;



*** summary***;
title2 "PP effects for *tevredenheid* in termen van *verschil* intervention minus control";
footnote " ";
data table_pp_eff;set eff_pp; where outcome="tevredenheid";run;
proc print data=table_pp_eff noobs; var outcome_label estimate lower upper probt icc icc_lower icc_upper;run;


title2 "PP effects in termen of ratio of intervention divided by control";
footnote "ICC is op log schaal";
data table_pp_eff;set eff_pp; where find(outcome,"log","i")> 0 and find(outcome,"log_tevr_revers","i") <= 0;
ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);run;
proc print data=table_pp_eff noobs; var outcome_label ratio ratio_ll ratio_ul probt icc icc_lower icc_upper;run;
footnote " ";

ods pdf close;





********** PP analyse met correctie voor LoS **************************;
data ds; set pp;run;
%let table_eff=eff_pp_corr;
%let table_descrip=descrip_pp_corr_cont;

ods pdf style=statistical file="EFFAMPART_st&date._PP_corr.pdf" compress=9; * style=journal kan ook;

proc odstext;
p "in de PP analyse (met correctie voor LoS) vergelijkbare effecten als in ITT";
run;
ods startpage=no;
* tevredenheid;
%sw_cont(outcome=tevredenheid,covars=log_los, outcome_min=0,outcome_max=15, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

*angst;
%sw_cont(outcome=log_angst, covars=log_los,outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

* depressie;
%sw_cont(outcome=log_depressie, covars=log_los,outcome_min=0, outcome_max=20,period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;


* ptss;
%sw_cont(outcome=log_ptss, covars=log_los,outcome_min=0, outcome_max=70, period_min=1, period_max=18,dseff=&table_eff, dsdescrip=&table_descrip);
ods startpage=yes;

title2 "PP effects *with* adjustment for log_LoS in terms of ratio of intervention divided by control";
data table_pp_eff_corr;set eff_pp_corr; where find(outcome,"log","i")> 0;
ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);run;
proc print data=table_pp_eff_corr noobs; var outcome_label ratio ratio_ll ratio_ul probt;run;

*** summary***;
title2 "PP effects for *tevredenheid* in termen van *verschil* intervention minus control";
footnote " ";
data table_pp_eff_corr;set eff_pp_corr; where outcome="tevredenheid";run;
proc print data=table_pp_eff_corr noobs; var outcome_label estimate lower upper probt icc icc_lower icc_upper;run;


title2 "PP effects in termen of ratio of intervention divided by control";
footnote "ICC is op log schaal";
data table_pp_eff_corr;set eff_pp_corr;; where find(outcome,"log","i")> 0 and find(outcome,"log_tevr_revers","i") <= 0;
ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);run;
proc print data=table_pp_eff_corr noobs; var outcome_label ratio ratio_ll ratio_ul probt icc icc_lower icc_upper;run;
footnote " ";



ods pdf close;


***********************************************;
********** word table with summaries **********;
title "EFFAMPART";
title "summary table";

* effect estimates;
data eff_itt; length descrip $ 100; set eff_itt; descrip="ITT";run;
data eff_itt_LOScorr; length descrip $ 100; set eff_itt_LoScorr; descrip="ITT, adjusted for log(LoS)";run;
data eff_pp; length descrip $ 100; set eff_pp; descrip="PP";run;
data eff_pp_corr; length descrip $ 100; set eff_pp_corr; descrip="PP, adjusted for log(LoS)";run;

data summary_effects; set eff_itt eff_itt_LOScorr eff_pp eff_pp_corr;
if find(outcome,"log", 'i') > 0 then do; 
	ratio=exp(estimate); ratio_ll=exp(lower); ratio_ul=exp(upper);
	end;
run;

* descriptions only for ITT and PP, both unadjusted;
data descrip_ITT_cont; length descrip $ 100; set descrip_ITT_cont; descrip="ITT";run;
data descrip_pp_cont;; length descrip $ 100; set descrip_pp_cont;; descrip="PP";run;
data summary_descrip; set descrip_ITT_cont descrip_pp_cont;run;


ods rtf file="summary_&date..rtf" style=minimal;

proc odstext;
p "descriptive statistics by population and by treatment";
run;
ods startpage=no;
proc tabulate data=summary_descrip order=data;
by descrip; where trt ne . ; *by population, only taking trt=0 and =1;
class outcome_label stat trt; var value; 
*by taking n instead of sum, you can check that only one value is taken, as should; 
table outcome_label,trt*(stat="")*(value="")*(sum="");
run;

ods startpage=yes;
ods startpage=yes;

proc odstext;
p "effect tevredenheid is een verschil op de oorspronkelijke schaal";
p "de ICC is op de oorspronkelijke schaal";
run;
ods startpage=no;
proc sort data=summary_effects; by outcome;run;
proc print data=summary_effects noobs;where outcome in ("tevredenheid"); by outcome;
  var descrip estimate lower upper probt icc icc_lower icc_upper; run;

ods startpage=yes;
proc odstext;
p "effect (angst, depressie, ptss) is een ratio op de oorspronkelijke schaal";
p "onderliggende data is log-getransformeerd";
p "de ICC is op de log schaal";
run;
ods startpage=no;
proc print data=summary_effects noobs;
  by outcome;
  where find(outcome, "log",'i')> 0 and outcome ne "log_tevr_revers"; 
  var descrip ratio ratio_ll ratio_ul probt icc icc_lower icc_upper; run;

ods rtf close;












/* weglaten tot nu toe


* delirium medication days: prevention;
title1 "prevention";
data ds1; set ds; where dlr_prev=1;label med_d='medication days (prevention)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (prevention)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table3, dsdescrip=table3_cont);

* delirium medication days: treatment;
title1 "treatment";
data ds1; set ds; where dlr_trt=1;label med_d='medication days (treatment)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (treatment)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table3, dsdescrip=table3_cont);

* delirium medication days: prevention and treatment;
title1 "prevention and treatment";
data ds1; set ds; where dlr_prevtrt=1;label med_d='medication days (treat+prev)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (treat+prev)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table3, dsdescrip=table3_cont);

title " ";
*delirium incidence;
%sw_bin(outcome=dlr, dseff=table3, dsdescrip=table3_bin); 

*mechanical ventilation days;
%sw_cont(outcome=vent_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
	* log;
%sw_cont(outcome=vent_d_log, outcome_min=0, outcome_max=4, dseff=table3, dsdescrip=table3_cont);

* re intubation ;
%sw_bin(outcome=retub, dseff=table3, dsdescrip=table3_bin); 

*readmission;
%sw_bin(outcome=readm,dseff=table3, dsdescrip=table3_bin); 

* unplanned removal devices;
%sw_bin(outcome=remov, dseff=table3, dsdescrip=table3_bin); 

* physical restraints;
%sw_bin(outcome=fixed, dseff=table3, dsdescrip=table3_bin); 

* physical restraints duration;
%sw_cont(outcome=fix_d, dseff=table3, dsdescrip=table3_cont);
	*log;
%sw_cont(outcome=fix_d_log, outcome_min=0, outcome_max=4, dseff=table3, dsdescrip=table3_cont); 

*LOS ICU;
%sw_cont(outcome=los_icu, outcome_min=1, outcome_max=130, dseff=table3, dsdescrip=table3_cont);
	*log;
%sw_cont(outcome=los_icu_log, outcome_min=0, outcome_max=5, dseff=table3, dsdescrip=table3_cont);

*LOS hosp;
%sw_cont(outcome=los_hosp, outcome_min=1, outcome_max=150, dseff=table3, dsdescrip=table3_cont);
	*log;
%sw_cont(outcome=los_hosp_log, outcome_min=0, outcome_max=5,dseff=table3, dsdescrip=table3_cont);

*D28;
%sw_bin(outcome=d28, dseff=table3, dsdescrip=table3_bin); 

*D90;
%sw_bin(outcome=d90, dseff=table3, dsdescrip=table3_bin); 

ods pdf close;


ods rtf style=minimal file="Table 3 ITT with trainingperiods.doc";
title "ITT: Efficacy Inference Table 3 with training periods";

data table3x; set table3; 
*add the exponentiated estimates and 95%-CI for log transformed variables;
if index(outcome_label,'log')>0 then do; 
e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
if label="trt"; * only the treatment estimates; 
run;
proc print data=table3x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; run;


title "ITT: Efficacy descriptives Table 3: continuous outcomes";
proc tabulate data=table3_cont order=data;
class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;

title "ITT: Efficacy descriptives Table 3: binary outcomes";
proc tabulate data=table3_bin order=data; class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
ods rtf close;


**********************************************************;
************* SUBGROUPS ANALYSES *************************;
**********************************************************;

title "Subgroup analyses Table 4"; 
%macro sw_subgr_cont(ds=ds, cat_min=0, cat_max=3, outcome=dcf,outcome_min=0, outcome_max=28,
interaction=admission_type, estimate_statement=, contrast_statement=, dsdescrip=table4_cont, dseff=table4_infer);
* requires numeric variable for subgroup and possibly with a value label;
%DO counter=&cat_min  %TO &cat_max;
data _null; set &ds;where &interaction=&counter;  
call symput('outcome_lab', vlabel(&outcome));
call symput('subgroup_name',vvalue(&interaction));
run;
data ds1; set &ds;where &interaction=&counter; 
label &outcome='&outcome_lab | subgroup &counter : &subgroup_name';run;
%sw_cont_descrip(ds=ds1,outcome=&outcome, outcome_min=&outcome_min, outcome_max=&outcome_max,dsdescrip=&dsdescrip);
%END;
title2 "interaction analysis";
%sw_cont_infer(ds=&ds,outcome=&outcome, covars_cat=&interaction, covars=&interaction, interaction=&interaction,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement,
dseff=&dseff);
%mend;


%macro sw_subgr_bin(ds=ds, cat_min=0, cat_max=3, outcome=dlr,
interaction=admission_type, estimate_statement=, contrast_statement=, dsdescrip=table4_bin, dseff=table4_infer);
* requires numeric variable for subgroup and possibly with a value label;
%DO counter=&cat_min  %TO &cat_max;
data _null; set &ds;where &interaction=&counter;  
call symput('outcome_lab', vlabel(&outcome)); *variable label;
call symput('subgroup_name',vvalue(&interaction));* value label;
run;
data ds1; set &ds;where &interaction=&counter; 
label &outcome='&outcome_lab | subgroup &counter : &subgroup_name';run;
%sw_bin_descrip(ds=ds1,outcome=&outcome, dsdescrip=&dsdescrip);
%END;
title2 "interaction analysis";
%sw_bin_infer(ds=&ds,outcome=&outcome, covars_cat=&interaction, covars=&interaction, interaction=&interaction,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement,
dseff=&dseff);
%mend;



ods pdf style=statistical compress=9 file="table4_adm_predeliric.pdf";
title2 "admission group";
%let interaction=admission_type;
%let estimate_statement=%str(estimate 'chirurgisch' trt 1 trt*admission_type 1 0 0 / cl;
estimate 'medisch(intern)' trt 1 trt*admission_type 0 1 0 / cl;
estimate 'trauma' trt 1 trt*admission_type 0 0 1 / cl;
);
%let contrast_statement=%str(contrast ' p-value interaction admission_type' trt*admission_type 1 0 -1, trt*admission_type 0 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=2, outcome=dcf, outcome_min=0, outcome_max=28, 
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=2, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=2, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=2, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=2, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);

title2 "pre-deliric"; 
proc format; 
*upper value not included, lower is included;
value dlr_cat    0-<  0.38 = '0'
				0.38-< 0.42 = '1'
				0.42-< 0.48 = '2'
				0.48-1.00  = '3'
;
run;

proc format;
value dlr_risk 	0='E-PRE-DELIRC: 0-<38%'
				1='E-PRE-DELIRC: 38-<42%'
				2='E-PRE-DELIRC: 42-<48%'
				3='E-PRE-DELIRC: 48-100%'
;
run;

data ds; set ds; 
dlr_risk=input(put(E_PRE_DELIRIC_Calc,dlr_cat.),best12.);
format dlr_risk dlr_risk.;
run;

%let interaction=dlr_risk;
%let estimate_statement=%str(estimate 'dlr risk: 0-<38%' trt 1 trt*dlr_risk 1 0 0 0/ cl;
estimate 'dlr risk: 38-<42%' trt 1 trt*dlr_risk 0 1 0 0 / cl;
estimate 'dlr risk: 42-<48%' trt 1 trt*dlr_risk 0 0 1 0 / cl;
estimate 'dlr risk: 48-100%' trt 1 trt*dlr_risk 0 0 0 1 / cl;
);
%let contrast_statement=%str(contrast ' p-valueinteraction dlr_risk' trt*dlr_risk 1 0 0 -1, trt*dlr_risk  0 1 0 -1, trt*dlr_risk 0 0 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dcf, outcome_min=0, outcome_max=28, 
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);

ods pdf close;


**;
ods pdf compress=9 style=statistical file="table4_apachLoS.pdf";
title2 "apache"; 
proc format; 
*upper value not included, lower is included;
value apache_cat    0-<  60 = '0'
				60 -< 80 = '1'
				80-< 100 = '2'
				100-high  = '3'
;
run;

proc format;
value apache_risk 	0='0=APACHE: 0-< 60'
					1='1=APACHE: 60-< 80'
					2='2=APACHE: 80-< 100'
					3='3=APACHE: 100 ->'
;
run;

data ds; set ds; 
apache_risk=input(put(APACHE_IV,apache_cat.),best12.);
format apache_risk apache_risk.;
run;


%let interaction=apache_risk;
%let estimate_statement=%str(estimate 'apache_risk:  0-<60' trt 1 trt*apache_risk 1 0 0 0/ cl;
estimate 'apache_risk:  60-<80' trt 1 trt*apache_risk 0 1 0 0 / cl;
estimate 'apache_risk:  80-100' trt 1 trt*apache_risk 0 0 1 0 / cl;
estimate 'apache_risk: 100-' trt 1 trt*apache_risk 0 0 0 1/ cl;
);
%let contrast_statement=%str(contrast ' p-value interaction apache_risk' trt*apache_risk 1 0 0 -1, trt*apache_risk 0 1 0 -1, trt*apache_risk 0 0 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dcf, outcome_min=0, outcome_max=28, 
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);


title2 "LOS_ICU risk"; 
proc format; 
*upper value not included, lower is included;
value los_cat    low-<  3 = '0'
				 3 -< 6 = '1'
				 6-< 12 = '2'
				 >12  = '3'
;
run;

proc format;
value los_risk 		0='0=LOS ICU:  0-< 3'
					1='1=LOS ICU:  3-< 6'
					2='2=LOS ICU:  6-< 12'
					3='3=LOS ICU: 12->'
;
run;

data ds; set ds; 
losicu_risk=input(put(los_icu,los_cat.),best12.);
format losicu_risk los_risk.;
run;


%let interaction=losicu_risk;
%let estimate_statement=%str(estimate 'losicu_risk:  0-3' trt 1 trt*losicu_risk 1 0 0 0/ cl;
estimate 'losicu_risk:  3-6' trt 1 trt*losicu_risk 0 1 0 0 / cl;
estimate 'losicu_risk:  6-12' trt 1 trt*losicu_risk 0 0 1 0 / cl;
estimate 'losicu_risk: 12-' trt 1 trt*losicu_risk 0 0 0 1/ cl;
);
%let contrast_statement=%str(contrast ' p-value interaction losicu_risk*' trt*losicu_risk 1 0 0 -1, trt*losicu_risk 0 1 0 -1, trt*losicu_risk 0 0 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dcf, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=3, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=3, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement);

ods pdf close;



 
ods pdf style=statistical file="Table4_surv.pdf";
title2 "D28 survival"; 
proc format;* format name must end with character;
value d28f	 	1='dead at day 28'
				0='alive at day 28'
;
run;


data ds; set ds;
label d28='dead at day 28';
format d28 d28f.;
run;
%let interaction=d28;
%let estimate_statement=%str(estimate 'd28 alive' trt 1 trt*d28 1 0 / cl;
estimate 'd28 dead' trt 1 trt*d28 0 1/ cl;
);
%let contrast_statement=%str(contrast ' p-value interaction died28' trt*d28 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dcf, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
*%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=1, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
); * no survival outcome when subgroups based on survival;

title2 "D90 survival"; 
proc format;* format name must end in character;
value d90f	 	1='dead at day 90'
					0='alive at day 90'
;
run;

data ds; set ds;
label d90='dead at day 90'; 
format d90 d90f.;
run;
%let interaction=d90;
%let estimate_statement=%str(estimate 'd90 alive' trt 1 trt*d90 1 0 / cl;
estimate 'd90 dead' trt 1 trt*d90 0 1/ cl;
);
%let contrast_statement=%str(contrast ' p-value interaction died90' trt*d90 1 -1;
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dcf, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr_d_log, outcome_min=0, outcome_max=3,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
%sw_subgr_cont(interaction=&interaction,cat_min=0, cat_max=1, outcome=dlr_d, outcome_min=0, outcome_max=28,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);
*%sw_subgr_bin(interaction=&interaction,cat_min=0, cat_max=1, outcome=d90,
estimate_statement=&estimate_statement, contrast_statement=&contrast_statement
);  * no survival outcome when subgroups based on survival;
* end;
ods pdf close;

	ods rtf style=minimal file="Table4_Subgroups.doc";
	title "Subgroups: Efficacy Inference Table 4";
data table4x; set table4_infer; 
*add the exponentiated estimates and 95%-CI for log transformed variables;
if index(outcome_label,'log')>0 then do; 
e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
if index(label, 'trt')=0; * not the rows with 'trt' in it; 
run;

* in the order of Table 4;
proc sort data=table4x; by label;
title2 "admission type";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'chirurgisch') > 0 or index(label,'medisch') > 0 or index(label,'trauma') > 0 or index(label, 'admission_type')> 0)
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;
title2 "delirium risk";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'dlr') > 0 )
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;
title2 "Apache score";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'apache') > 0 )
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;
title2 "LOS icu";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'losicu') > 0 )
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;
title2 "d28 survival status";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'d28') > 0 )
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;
title2 "d90 survival status";
proc print data=table4x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; 
where (index(label,'d90') > 0 )
and outcome_label ne "delirium incidence" and outcome_label ne "delirium (days)" and outcome_label ne "d90";
run;

	title "Subgroups: Efficacy descriptives Table 4: continuous outcomes";
proc tabulate data=table4_cont order=data; class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
	title "Subgroups: Efficacy descriptives Table 4: binary outcomes";
proc tabulate data=table4_bin order=data; class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
	ods rtf close;




**********************************************************;
************* PP ANALYSES *************************;
**********************************************************;
ods pdf compress=9 style=statistical file="Table5_PP_analyses.doc";

title1 "PP: per_protocol_analysis=1 and trainings periods removed"; 
data ds; set itt; 
if Per_protocol_analysis=1;;
*exclude training periods;
if center=1 and period=5 then delete;
if center=2 and period=2 then delete;
if center=3 and period=13 then delete;
if (center=4 or center=6) and period=7 then delete;
if center=5 and period=11 then delete;
if center=7 and period=9 then delete;
if center=8 and period=12 then delete;
if center=9 and period=10 then delete;
if center=10 and period=3 then delete;
run; 

title2 "design with treatment assignment";
%rm_design(ds=ds);

title2 " ";
* delirium-coma free;
%sw_cont(outcome=dcf,outcome_min=0,outcome_max=28, dseff=table5, dsdescrip=table5_cont);

*delirium days;
%sw_cont(outcome=dlr_d, outcome_min=0, outcome_max=28,dseff=table5, dsdescrip=table5_cont);
	*log;
%sw_cont(outcome=dlr_d_log, outcome_min=0, outcome_max=4,dseff=table5, dsdescrip=table5_cont);

* coma days;
%sw_cont(outcome=coma_d, outcome_min=0, outcome_max=28, dseff=table5, dsdescrip=table5_cont);
	* log;
%sw_cont(outcome=coma_d_log, outcome_min=0, outcome_max=4, dseff=table5, dsdescrip=table5_cont);

* sedation days; * increase!;
%sw_cont(outcome=sed_d, outcome_min=0, outcome_max=28, dseff=table5, dsdescrip=table5_cont);
	*log;
%sw_cont(outcome=sed_d_log, outcome_min=0, outcome_max=4, dseff=table5, dsdescrip=table5_cont);

* delirium medication days: prevention;
title1 "prevention";
data ds1; set ds; where dlr_prev=1;label med_d='medication days (prevention)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (prevention)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table5, dsdescrip=table5_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table5, dsdescrip=table5_cont);

* delirium medication days: treatment;
title1 "treatment";
data ds1; set ds; where dlr_trt=1;label med_d='medication days (treatment)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (treatment)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table5, dsdescrip=table5_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table5, dsdescrip=table5_cont);

* delirium medication days: prevention and treatment;
title1 "prevention and treatment";
data ds1; set ds; where dlr_prevtrt=1;label med_d='medication days (treat+prev)';
med_d_log=log(med_d+1);label med_d_log='log medication days (+1) (treat+prev)';run;
%sw_cont(ds=ds1, outcome=med_d, outcome_min=0, outcome_max=28,dseff=table5, dsdescrip=table5_cont);
%sw_cont(ds=ds1, outcome=med_d_log, outcome_min=0, outcome_max=4,dseff=table5, dsdescrip=table5_cont);

title " ";
*delirium incidence;
%sw_bin(outcome=dlr, dseff=table5, dsdescrip=table5_bin); 

*mechanical ventilation days;
%sw_cont(outcome=vent_d, outcome_min=0, outcome_max=28,dseff=table5, dsdescrip=table5_cont);
	* log;
%sw_cont(outcome=vent_d_log, outcome_min=0, outcome_max=4, dseff=table5, dsdescrip=table5_cont);

* re intubation ;
%sw_bin(outcome=retub, dseff=table5, dsdescrip=table5_bin); 

*readmission;
%sw_bin(outcome=readm,dseff=table5, dsdescrip=table5_bin); 

* unplanned removal devices;
%sw_bin(outcome=remov, dseff=table5, dsdescrip=table5_bin); 

* physical restraints;
%sw_bin(outcome=fixed, dseff=table5, dsdescrip=table5_bin); 

* physical restraints duration;
%sw_cont(outcome=fix_d, dseff=table5, dsdescrip=table5_cont);
	*log;
%sw_cont(outcome=fix_d_log, outcome_min=0, outcome_max=4, dseff=table5, dsdescrip=table5_cont); 

*LOS ICU;
%sw_cont(outcome=los_icu, outcome_min=1, outcome_max=130, dseff=table5, dsdescrip=table5_cont);
	*log;
%sw_cont(outcome=los_icu_log, outcome_min=0, outcome_max=5, dseff=table5, dsdescrip=table5_cont);

*LOS hosp;
%sw_cont(outcome=los_hosp, outcome_min=1, outcome_max=150, dseff=table5, dsdescrip=table5_cont);
	*log;
%sw_cont(outcome=los_hosp_log, outcome_min=0, outcome_max=5,dseff=table5, dsdescrip=table5_cont);

*D28;
%sw_bin(outcome=d28, dseff=table5, dsdescrip=table5_bin); 

*D90;
%sw_bin(outcome=d90, dseff=table5, dsdescrip=table5_bin); 


* end Table 5;
ods pdf close;


ods rtf style=minimal file="Table 5.doc";
title "PP: Efficacy Inference Table 5";

data table5x; set table5; 
*add the exponentiated estimates and 95%-CI for log transformed variables;
if index(outcome_label,'log')>0 then do; 
e_est=exp(estimate);e_lower=exp(lower); e_upper=exp(upper);end;
if label="trt"; * only the treatment estimates; 
run;
proc print data=table5x noobs; var outcome outcome_label label estimate lower upper probt e_est e_lower e_upper; run;


title "ITT: Efficacy descriptives Table 5: continuous outcomes";
proc tabulate data=table5_cont order=data;
class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;

title "ITT: Efficacy descriptives Table 5: binary outcomes";
proc tabulate data=table5_bin order=data; class outcome outcome_label stat trt;var value;
table outcome*outcome_label, trt*(stat=' ')*(value=' ')*(mean=' ')*f=7.2;run;
ods rtf close;
 


**********************************************************;
************* Forest plot ********************************;
**********************************************************;

* use the itt set again;
data ds; set itt;run;


title "Subgroup analyses Figure"; 

%macro forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=, estimate_statement=%str( ),  dsforest=forest);
* note the output is appended to &dsforest;

*run the interaction analysis with appropriate estimate statements to get the estimates in the subgroups;
%sw_cont_infer(ds=&ds,outcome=&outcome, cluster=&cluster, period=&period, trt=&trt,
period_min=&period_min,period_max=&period_max,covars_cat=&interaction, covars=&interaction, 
interaction=&interaction, dseff=_subgroups, random_statement=%str(random intercept /subject=&cluster; ),
repeated_statement=%str( ), estimate_statement=&estimate_statement, contrast_statement=%str( )
);

%IF &dsforest ne %THEN %DO;
	* if subgroups are requested;
	%IF &interaction ne %THEN %DO;
		* get the variable label of the interaction variable from the master dataset;
		data _null_; set ds(obs=1 keep=&interaction); call symput('label', vlabel(&interaction));run;
		* first record to state variable determining the subgroups;
		data _firstrow; length label $50;length outcome $30;length outcome_label $100; length model $100;
		level=1; outcome="&outcome";label="&label";estimate=.;lower=.; upper=.;stderr=.;probt=.;run;
		* lines with subgroup values (in the label of the estimate statement;
		data _otherrows; set _subgroups;level=2;
		if index(label, "&trt")=0; * keep only the estimate statement records;
		run;
 		* apppend datasets;
		data _allrows; set _firstrow _otherrows;run;
		proc append force base=&dsforest data=_allrows;run;
		%END;
	%ELSE %DO;
	* if no subgroups are requested (e.g. for overall unadjusted estimate) pick from estimate statement;
	* and put at level=1 for the table;
		data _otherrows; set _subgroups;level=1;
		if index(label, "&trt")=0; * keep only the estimate statement records;
		run;
 		* apppend datasets;
		data _allrows; set _otherrows;run;
		proc append force base=&dsforest data=_allrows;run;	
		%END;
	* remove temporary dataset;
	proc delete data=_subgroups;run;
	%END;
%mend;



%macro forest_bin(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=, estimate_statement=,  dsforest=forest);
* note the output is appended to &dsforest;

*run the interaction analysis with appropriate estimate statements to get the estimates in the subgroups;
%sw_bin_infer(ds=&ds, outcome=&outcome,cluster=&cluster, period=&period, trt=&trt,period_min=&period_min, 
covars_cat=&interaction, covars=&interaction, interaction=&interaction, period_max=&period_max, 
random_statement=%str(random intercept /subject=&cluster; ),repeated_statement=%str( ),
estimate_statement=&estimate_statement, contrast_statement=%str( ),dseff=_subgroups); 

%IF &dsforest ne %THEN %DO;
	* if subgroups are requested;
	%IF &interaction ne %THEN %DO;
		* get the variable label of the interaction variable from the master dataset;
		data _null_; set ds(obs=1 keep=&interaction); call symput('label', vlabel(&interaction));run;
		* first record to state variable determining the subgroups;
		data _firstrow; length label $50;length outcome $30;length outcome_label $100; length model $100;
		level=1; outcome="&outcome";label="&label";estimate=.;lower=.; upper=.;stderr=.;probt=.;run;
		* lines with subgroup values (in the label of the estimate statement;
		data _otherrows; set _subgroups;level=2;
		if index(label, "&trt")=0; * keep only the estimate statement records;
		run;
 		* apppend datasets;
		data _allrows; set _firstrow _otherrows;run;
		proc append force base=&dsforest data=_allrows;run;
		%END;
	%ELSE %DO;
	* if no subgroups are requested (e.g. for overall unadjusted estimate) pick from estimate statement;
	* and put at level=1 for the table;
		data _otherrows; set _subgroups;level=1;
		if index(label, "&trt")=0; * keep only the estimate statement records;
		run;
 		* apppend datasets;
		data _allrows; set _otherrows;run;
		proc append force base=&dsforest data=_allrows;run;	
		%END;
	* remove temporary dataset;
	proc delete data=_subgroups;run;
	%END;
%mend;


title2 "overall estimates (not by subgroup)";
%let interaction=; *no interaction variable;
%let estimate_statement=%str(estimate 'Overall' trt 1 /cl;);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);



title2 "admission group";
%let interaction=admission_type;
%let estimate_statement=%str(estimate 'surgical' trt 1 trt*admission_type 1 0 0 / cl;
estimate 'medical' trt 1 trt*admission_type 0 1 0 / cl;
estimate 'trauma' trt 1 trt*admission_type 0 0 1 / cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);


title2 "pre-deliric"; 
proc format; 
*upper value not included, lower is included;
value dlr_cat    0-<  0.38 = '0'
				0.38-< 0.42 = '1'
				0.42-< 0.48 = '2'
				0.48-1.00  = '3'
;
run;

proc format;
value dlr_risk 	0='E-PRE-DELIRC: 0-<38%'
				1='E-PRE-DELIRC: 38-<42%'
				2='E-PRE-DELIRC: 42-<48%'
				3='E-PRE-DELIRC: 48-100%'
;
run;

data ds; set ds; 
dlr_risk=input(put(E_PRE_DELIRIC_Calc,dlr_cat.),best12.);
label dlr_risk="Delirium risk";
format dlr_risk dlr_risk.;
run;

%let interaction=dlr_risk;
%let estimate_statement=%str(estimate '0-<38%' trt 1 trt*dlr_risk 1 0 0 0/ cl;
estimate '38-<42%' trt 1 trt*dlr_risk 0 1 0 0 / cl;
estimate '42-<48%' trt 1 trt*dlr_risk 0 0 1 0 / cl;
estimate '48-100%' trt 1 trt*dlr_risk 0 0 0 1 / cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);


title2 "apache"; 
proc format; 
*upper value not included, lower is included;
value apache_cat    0-<  60 = '0'
				60 -< 80 = '1'
				80-< 100 = '2'
				100-high  = '3'
;
run;

proc format;
value apache_risk 	0='0=APACHE: 0-< 60'
					1='1=APACHE: 60-< 80'
					2='2=APACHE: 80-< 100'
					3='3=APACHE: 100 ->'
;
run;

data ds; set ds; 
apache_risk=input(put(APACHE_IV,apache_cat.),best12.);
format apache_risk apache_risk.;
label apache_risk="APACHE IV";
run;


%let interaction=apache_risk;
%let estimate_statement=%str(estimate '0-<60' trt 1 trt*apache_risk 1 0 0 0/ cl;
estimate '60-<80' trt 1 trt*apache_risk 0 1 0 0 / cl;
estimate '80-100' trt 1 trt*apache_risk 0 0 1 0 / cl;
estimate '100-' trt 1 trt*apache_risk 0 0 0 1/ cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);


title2 "LOS_ICU risk"; 
*** note the fit here is less, residuals are systematically higher, so model tends to underestimate observed***;
proc format; 
*upper value not included, lower is included;
value los_cat    low-<  3 = '0'
				 3 -< 6 = '1'
				 6-< 12 = '2'
				 >12  = '3'
;
run;

proc format;
value los_risk 		0='0=LOS ICU: 0-< 3'
					1='1=LOS ICU: 3-< 6'
					2='2=LOS ICU: 6-< 12'
					3='3=LOS ICU: 12->'
;
run;

data ds; set ds; 
losicu_risk=input(put(los_icu,los_cat.),best12.);
format losicu_risk los_risk.;
label losicu_risk="LoS ICU";
run;


%let interaction=losicu_risk;
%let estimate_statement=%str(estimate '0-3 days' trt 1 trt*losicu_risk 1 0 0 0/ cl;
estimate '3-6 days' trt 1 trt*losicu_risk 0 1 0 0 / cl;
estimate '6-12 days' trt 1 trt*losicu_risk 0 0 1 0 / cl;
estimate '> 12 days' trt 1 trt*losicu_risk 0 0 0 1/ cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);



 
title2 "D28 survival"; 
** note: D90 survival did not converge ***;
proc format;* format name must end with character;
value d28f	 	1='dead'
				0='alive'
;
run;


data ds; set ds;
label d28='Survival status at day 28';
format d28 d28f.;
run;
%let interaction=d28;
%let estimate_statement=%str(estimate 'alive' trt 1 trt*d28 1 0 / cl;
estimate 'dead' trt 1 trt*d28 0 1/ cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);

title2 "D90 survival"; 
** note D90 as outcome did not converge because of infitite likelihood **;
proc format;* format name must end in character;
value d90f	 	1='dead'
					0='alive'
;
run;

data ds; set ds;
label d90='Survival status at day 90'; 
format d90 d90f.;
run;
%let interaction=d90;
%let estimate_statement=%str(estimate 'alive' trt 1 trt*d90 1 0 / cl;
estimate 'dead' trt 1 trt*d90 0 1/ cl;
);
%forest_cont(ds=ds, outcome=dcf,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_cont(ds=ds, outcome=dlr_d_log,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=dlr,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);
%forest_bin(ds=ds, outcome=d90,cluster=center, period=period, trt=trt,
period_min=1,period_max=14, interaction=&interaction, estimate_statement=&estimate_statement,  dsforest=forest);

* save forest data;
data dir.forest;set forest;run;


* needed for the macro forest_plot: left right labels;
proc format;
   value $interpretationfmt
     "L" = "<== intervention better"  
     "R" = "control better ==>"    
      ;                             
run;

%macro forest_plot(dsforest=dir.forest, outcome=, statistic=' ', vrefline=,x_L=-2, x_R=2,xmin=,xmax=);
data forest_outcome; set &dsforest(where=(outcome="&outcome")); 
run;
* add empty row at the bottom, and add numbering;
data empty_row;
label=" ";outcome=" ";outcome_label=" "; model=" ";estimate=.;lower=.; upper=.;stderr=.;probt=.;level=1; 
run;
data forest_outcome; set forest_outcome empty_row;
record=_n_;
run;
* add interpretation to x-as;
data interpretation; length text $ 25; 
format text $interpretationfmt.;
x1=&x_L; text="L"; output;
x1=&x_R; text="R"; output;
x1=0; text=" "; output; *empty line;
run;
data forest_outcome; set forest_outcome interpretation; record=_n_;run;
	* the interpretation records must have the same record value;
	data forest_outcome; set forest_outcome nobs=nobs;
		if _n_ ge nobs-1 then record=record-1;
	run;
* no datalines/cards statement allowed within macro code;
* see https://communities.sas.com/t5/SAS-Programming/DATALINES-statement-inside-a-macro/td-p/37960;
* therefore make the dataset record by record;
data attrmap;
* make sure variable values are not truncated to the variable length in the first record length;
length  texcolor $ 10; length textweight $ 10;
id="text"; value=1;texcolor="black"; textsize=20;textweight="bold";output;
id="text"; value=2;texcolor="blue"; textsize=5;textweight="normal";output;
run;

* get the value for the overall reference line and the label of the variable;
%local overallvalue;%local outcome_label;
data _null_; set forest_outcome(where=(label="Overall"));
call symput('overallvalue', estimate);
call symput('outcome_label', outcome_label);
run;

title1 " "; * clear titles;
proc sgplot data=forest_outcome dattrmap=attrmap noborder nowall noautolegend ; 
label label="&outcome_label";
  *--- estimates and CIs ---;
   scatter y=record x=estimate / 
      markerattrs=(symbol=squarefilled)
      ;
   highlow y=record low=lower high=upper;
  *--- adding yaxis table at left ---;
   yaxistable label / 
      location=inside
      position=left
	  textgroup=level
      textgroupid=text
      ;
   *--- primary axes ---;
   yaxis
      reverse
      display=none
      offsetmin=0.05
      ; * tune with offsetmin between 0 and 1 how much space between top and first 'subgroup';
 *--- cleaner x axis ---;
   xaxis 
      display=(nolabel) 
	  %IF &xmin ne %THEN min=&xmin;%IF &xmax ne %THEN max=&xmax;
     ;
  *---interpretation above x-axis---;
	  text x=x1 y=record text=text / 
      position=bottom 
      contributeoffsets=none 
      strip
      ;
   *--- text above x2axis to denote the summary statistic ---;
   * use a second scatter with 0 size symbols ;
   * to be able to use the x2axis statement;
   scatter y=record x=estimate / 
      markerattrs=(size=0) 
      x2axis
      ;
   x2axis 
      label=&statistic 
      display=(noline noticks novalues) 
      ;
*---- refline overall estimate ---;
	%IF &vrefline eq %THEN %DO; refline &overallvalue /axis=x;%END;
	%IF &vrefline ne %THEN %DO; refline &vrefline /axis=x; %END;
run;
%mend;

* options mprint mlogic symbolgen ;






ods rtf file="forest_plots_st200517.doc";
libname dir ".";
data ds; set dir.forest;
* keep only the dcf, dlr_d_log, dlr, and d90;
if outcome in ("dcf", "dlr_d_log", "dlr", "d90"); 
* exponentiate the log values (logodds becomes odds ratio , log becomes median ratio);
if index(outcome_label, "log") > 0 then do; estimate=exp(estimate); lower=exp(lower); upper=exp(upper);end;
* keep only the log(odds) records for dlr and d90;
If outcome="dlr" and index(outcome_label,"log(odds)") =0 and compress(model) ne " " then delete;
If outcome="d90" and index(outcome_label,"log(odds)") = 0 and compress(model) ne " " then delete;
* do not analyze d28 and d90 for the outcome d90;
if outcome = "d90" and label in ("alive","dead","Dead at day 90","Dead at day 28") then delete;
* variabele namen hernoemen;
if outcome="dlr" then outcome_label="delirium incidence";
if outcome="d90" then outcome_label="Died at day 90";
if outcome="dlr_d_log" then outcome_label="delirium duration (days)";
run;

* dcf: more delier/coma free is better;
proc format;
   value $interpretationfmt
     "L" = "<== control better"  
     "R" = "intervention better ==>"    
      ;                             
run;
%forest_plot(dsforest=ds, outcome=dcf, statistic='Difference in days', vrefline=0);
* dlr_d_log: less delier is better;
proc format;
   value $interpretationfmt
     "L" = "<==  intervention better"  
     "R" = "control better ==>"    
      ;                             
run;
%forest_plot(dsforest=ds, outcome=dlr_d_log, statistic='Median ratio', vrefline=1, x_L=0.5, x_R=1.5);
* dlr: less delier is better;
proc format;
   value $interpretationfmt
     "L" = "<==  intervention better"  
     "R" = "control better ==>"    
      ;                             
run;
%forest_plot(dsforest=ds, outcome=dlr, statistic='Odds ratio',vrefline=1, x_L=-1, x_R=3,xmin=-2);
* dead90: less dead is better;
proc format;
   value $interpretationfmt
     "L" = "<==  intervention better"  
     "R" = "control better ==>"    
      ;                             
run;
%forest_plot(dsforest=ds, outcome=d90, statistic='Odds ratio',vrefline=1,x_L=0.5, x_R=1.6);

* print the estimates etc;
proc sort data=ds; by outcome; 
proc print data=ds noobs; by outcome;var label outcome_label estimate lower upper;run;  

ods rtf close;



**************************************************************;
****** dcf and dlr_d in the patient with delier (ITT and PP)**; 
**************************************************************;
data itt_delier; set itt; if dlr=1;run;
data pp_delier; set itt; if dlr=1; 
if Per_protocol_analysis=1;;
*exclude training periods;
if center=1 and period=5 then delete;
if center=2 and period=2 then delete;
if center=3 and period=13 then delete;
if (center=4 or center=6) and period=7 then delete;
if center=5 and period=11 then delete;
if center=7 and period=9 then delete;
if center=8 and period=12 then delete;
if center=9 and period=10 then delete;
if center=10 and period=3 then delete;
run; 

ods rtf file="underpin_onlydelier.doc";
title1 "ITT patients with delier (dlr=1)";
title2 "as selection is on a post-baseline event: selection bias possible";
* delier and coma free days alive;
%sw_cont(ds=itt_delier,outcome=dcf,outcome_min=0,outcome_max=28, dseff=table3, dsdescrip=table3_cont);

*delirium days;
%sw_cont(ds=itt_delier,outcome=dlr_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
	*log;
%sw_cont(ds=itt_delier,outcome=dlr_d_log, outcome_min=0, outcome_max=4,dseff=table3, dsdescrip=table3_cont);

title1 "PP patients (no trainig periods and per_protocol_analysis=1) with delier";
title2 "as selection is on a post-baseline event: selection bias possible";
%sw_cont(ds=pp_delier,outcome=dcf,outcome_min=0,outcome_max=28, dseff=table3, dsdescrip=table3_cont);

*delirium days;
%sw_cont(ds=pp_delier,outcome=dlr_d, outcome_min=0, outcome_max=28,dseff=table3, dsdescrip=table3_cont);
	*log;
%sw_cont(ds=pp_delier,outcome=dlr_d_log, outcome_min=0, outcome_max=4,dseff=table3, dsdescrip=table3_cont);

ods rtf close;
**/

