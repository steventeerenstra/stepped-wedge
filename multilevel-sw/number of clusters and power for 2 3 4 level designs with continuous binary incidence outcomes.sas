***************************************************************************************************************;
***** 6 macros to calculate the number of clusters or the power for continuous, binary and rate  outcomes *****;
***************************************************************************************************************;
** terminology n1,n2,n3,n4, rho12, rho23, rho34, sigma2tot as in the paper **;
** various options for calculating the sigma2tot are given: 
** binary: based on p0 and p1, or based on the p closest to 0.5 in the control and in intervention condition 
           or based  on from user input **;
** rate: based on lambda0 and lambda1 or based on the maximum rate in the control and intervention condition 
           or based  on from user input **;
**note that the levels that are measured as cohort are given as 4,43,432 or 3, 32, or 2 ***;
**so a 4 level design with level 4, 3 and 2 measured as cohort is denoted as: cohort=432;
** examples to do the calculations as for the paper are found below **;



* debugging options: remove "*" in next line to get debugging code;
*options macrogen symbolgen source2 mlogic;

%macro SWpower_cont(steps=2, cohort=43, n1=1,n2=1,n3=1,n4=1,rho12=1,rho23=1,rho34=1,mu1=0, mu0=1,sigma2tot=1,alpha=0.05);
data _sw; 
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;n4=&n4;rho12=&rho12;rho23=&rho23; rho34=&rho34;
mu1=&mu1;mu0=&mu0;sigma2tot=&sigma2tot;
* calculate VIFs, total sample sizes and rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
Ntot4=n1*n2*n3*n4; Ntot3=n1*n2*n3; Ntot2=n1*n2;
* calculate rho according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;Ntot=Ntot4; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;Ntot=Ntot4;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4; vif_ml=vif4;Ntot=Ntot4;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;Ntot=Ntot2;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
var_delta=4*sigma2tot*vif_sw/Ntot;
delta=mu1-mu0;
z_alpha=probit(1-&alpha/2);
power=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
run;
proc print data=_sw noobs; 
var sigma2tot cohort s n1 n2 n3 n4 rho12 rho23 rho34 mu1 mu0 delta vif_ml rho vif_swPG1 vif_sw var_delta power; 
run;
%mend;


%macro SWpower_bin(steps=2, cohort=43, n1=1,n2=1,n3=1,n4=1,rho12=1,rho23=1,rho34=1,p1=0.2, 
p0=0.35,p1_closest05=,p0_closest05=,sigma2tot_userinput=, alpha=0.05);
* choices for sigma2tot;
data _sigma2tot;
length var_method $70.;
var_method="based on p0 and p1"; sigma2tot=(&p1*(1-&p1)+&p0*(1-&p0))/2;output;
%IF %sysevalf(&p0_closest05 ne  AND &p1_closest05 ne, boolean) %THEN 
%DO; 
	var_method="based on p0 closest to 0.5 (&p0_closest05) and p1 closest to 0.5 (&p1_closest05)"; 
	sigma2tot=(&p0_closest05*(1-&p0_closest05)+&p1_closest05*(1-&p1_closest05))/2;
	output;
%END;
%IF %sysevalf(&sigma2tot_userinput ne , boolean) %THEN
%DO;
	var_method="user input: &sigma2tot_userinput";
	sigma2tot=&sigma2tot_userinput;
	output;
%END;
run;

data _sw; set _sigma2tot;
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;n4=&n4;rho12=&rho12;rho23=&rho23; rho34=&rho34;
p1=&p1;p0=&p0;
* calculate VIFs, total sample sizes and rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
Ntot4=n1*n2*n3*n4; Ntot3=n1*n2*n3; Ntot2=n1*n2;
* calculate rho according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;Ntot=Ntot4; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;Ntot=Ntot4;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4;  vif_ml=vif4;Ntot=Ntot4;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;Ntot=Ntot2;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
var_delta=4*sigma2tot*vif_sw/Ntot;
delta=p1-p0;
z_alpha=probit(1-&alpha/2);
power=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
run;
proc print data=_sw noobs; 
var var_method sigma2tot cohort s n1 n2 n3 n4 rho12 rho23 rho34 p1 p0 delta vif_ml rho vif_swPG1 vif_sw var_delta power; 
run;
%mend;

%macro SWpower_rate(steps=4, cohort=32, n1=10,n2=4,n3=116,n4=1,rho12=0.7,rho23=0.01,rho34=1,lambda1=0.005, 
lambda0=0.011,lambda0_max=,lambda1_max=,sigma2tot_userinput=, alpha=0.05);
* choices for sigma2tot;
data _sigma2tot;
length var_method $70.;
* unlike the binomial case, the expected total variance has the correction for 1/(1-rho12) already; 
var_method="based on lambda0 and lambda1"; sigma2tot=(&lambda1+&lambda0)/(2*(1-&rho12));output;
%IF %sysevalf(&lambda0_max ne  AND &lambda1_max ne, boolean) %THEN 
%DO; 
	var_method="based on lambda0_max (&lambda0_max) and lambda1_max (&lambda1_max)"; 
	sigma2tot=(&lambda0_max+&lambda1_max)/(2*(1-&rho12));
	output;
%END;
%IF %sysevalf(&sigma2tot_userinput ne , boolean) %THEN
%DO;
	var_method="user input: &sigma2tot_userinput";
	sigma2tot=&sigma2tot_userinput;
	output;
%END;
run;

data _sw; set _sigma2tot;
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;n4=&n4;rho12=&rho12;rho23=&rho23; rho34=&rho34;
lambda1=&lambda1;lambda0=&lambda0;
* calculate VIFs, total sample sizes and rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
Ntot4=n1*n2*n3*n4; Ntot3=n1*n2*n3; Ntot2=n1*n2;
* calculate rho according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;Ntot=Ntot4; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;Ntot=Ntot4;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4;  vif_ml=vif4;Ntot=Ntot4;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;Ntot=Ntot3;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;Ntot=Ntot2;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
var_delta=4*sigma2tot*vif_sw/Ntot;
delta=lambda1-lambda0;
z_alpha=probit(1-&alpha/2);
power=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
run;
proc print data=_sw noobs; 
var var_method sigma2tot cohort s n1 n2 n3 n4 rho12 rho23 rho34 lambda1 lambda0 delta vif_ml rho vif_swPG1 vif_sw var_delta power; 
run;

%mend;

%macro SWnclusters_cont(steps=4, cohort=43, n1=1,n2=1,n3=1,rho12=1,rho23=1,rho34=1,mu1=0, mu0=1,sigma2tot=1,alpha=0.05, power=0.80);
data _sw; 
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;;rho12=&rho12;rho23=&rho23; rho34=&rho34;
mu1=&mu1;mu0=&mu0;sigma2tot=&sigma2tot; power_aim=&power;
*calculate uncorrected total sample size Ntot0 for a parallel group post-test trial;
delta=mu1-mu0;
z_alpha=probit(1-&alpha/2);
z_beta=probit(&power);
Ntot_pg1=2*2*(z_alpha+z_beta)**2*&sigma2tot/(delta**2);
* calculate VIFs using rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
* calculate rho and cluster size according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;clussize=n1*n2*n3; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;clussize=n1*n2*n3;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4;  vif_ml=vif4;clussize=n1*n2*n3;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;clussize=n1;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
Ntot=Ntot_pg1*vif_sw;
n_clusters_raw=Ntot/clussize; * raw number of clusters;
n_clusters=max(ceil(n_clusters_raw),s); *number of cluster has to be integer and minimal one cluster per step;
*calculating power based on estimated number of clusters;
var_delta=4*sigma2tot*vif_sw/(clussize*n_clusters);
var_delta_raw=4*sigma2tot*vif_sw/(clussize*n_clusters_raw);
power_check=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
*power_check_raw=probnorm(abs(delta)/sqrt(var_delta_raw)-z_alpha);
run;
proc print data=_sw noobs; 
var cohort s n1 n2 n3 rho12 rho23 rho34 mu1 mu0 sigma2tot delta Ntot_pg1 vif_ml rho vif_swPG1 vif_sw Ntot n_clusters power_aim power_check;
run;
%mend;


%macro SWnclusters_bin(steps=4, cohort=43, n1=1,n2=1,n3=1,rho12=1,rho23=1,rho34=1,p1=0, p0=1,p1_closest05=,p0_closest05=, 
sigma2tot_userinput=,alpha=0.05, power=0.80);
* choices for sigma2tot;
data _sigma2tot;
length var_method $70.;
var_method="based on p0 and p1"; sigma2tot=(&p1*(1-&p1)+&p0*(1-&p0))/2;output;
%IF %sysevalf(&p0_closest05 ne  AND &p1_closest05 ne, boolean) %THEN 
%DO; 
	var_method="based on p0 closest to 0.5 (=&p0_closest05) and p1 closest to 0.5 (=&p1_closest05)"; 
	sigma2tot=(&p0_closest05*(1-&p0_closest05)+&p1_closest05*(1-&p1_closest05))/2;
	output;
%END;
%IF %sysevalf(&sigma2tot_userinput ne , boolean) %THEN
%DO;
	var_method="user input: &sigma2tot_userinput";
	sigma2tot=&sigma2tot_userinput;
	output;
%END;
run;

data _sw; set _sigma2tot;
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;;rho12=&rho12;rho23=&rho23; rho34=&rho34;
p1=&p1;p0=&p0; power_aim=&power;
*calculate uncorrected total sample size Ntot0 for a parallel group post-test trial;
delta=p1-p0;
z_alpha=probit(1-&alpha/2);
z_beta=probit(&power);
Ntot_pg1=2*2*(z_alpha+z_beta)**2*sigma2tot/(delta**2);
* calculate VIFs using rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
* calculate rho and cluster size according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;clussize=n1*n2*n3; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;clussize=n1*n2*n3;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4;  vif_ml=vif4;clussize=n1*n2*n3;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;clussize=n1;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
Ntot=Ntot_pg1*vif_sw;
n_clusters_raw=Ntot/clussize; * raw number of clusters;
n_clusters=max(ceil(n_clusters_raw),s); *number of cluster has to be integer and minimal one cluster per step;
*calculating power based on estimated number of clusters;
var_delta=4*sigma2tot*vif_sw/(clussize*n_clusters);
var_delta_raw=4*sigma2tot*vif_sw/(clussize*n_clusters_raw);
power_check=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
*power_check_raw=probnorm(abs(delta)/sqrt(var_delta_raw)-z_alpha);
run;
proc print data=_sw noobs; 
var cohort s n1 n2 n3 rho12 rho23 rho34 p1 p0 var_method sigma2tot delta Ntot_pg1 vif_ml rho vif_swPG1 vif_sw Ntot n_clusters power_aim power_check;
run;
%mend;



%macro SWnclusters_rate(steps=4, cohort=43, n1=1,n2=1,n3=1,rho12=1,rho23=1,rho34=1,lambda1=0.005, lambda0=0.011,
lambda1_max=, lambda0_max=, sigma2tot_userinput=,alpha=0.05, power=0.80);
* choices for sigma2tot;
data _sigma2tot;
length var_method $70.;
var_method="based on lambda0 and lambda1";
* unlike the binomial case, the expected total variance has the correction for 1/(1-rho12) already; 
sigma2tot=(&lambda1+&lambda0)/(2*(1-&rho12));output;
%IF %sysevalf(&lambda1_max ne  AND &lambda0_max ne, boolean) %THEN 
%DO; 
	var_method="based on lambda0 max (=&lambda0_max) and lambda1 max (=&lambda1_max)"; 
	sigma2tot=(&lambda1_max + &lambda0_max)/(2*(1-&rho12));
	output;
%END;
%IF %sysevalf(&sigma2tot_userinput ne , boolean) %THEN
%DO;
	var_method="user input: &sigma2tot_userinput";
	sigma2tot=&sigma2tot_userinput;
	output;
%END;
run;

data _sw; set _sigma2tot;
*** set all parameters **;
cohort=&cohort; s=&steps;
n1=&n1; n2=&n2;n3=&n3;;rho12=&rho12;rho23=&rho23; rho34=&rho34;
lambda1=&lambda1;lambda0=&lambda0; power_aim=&power;
*calculate uncorrected total sample size Ntot0 for a parallel group post-test trial;
delta=lambda1-lambda0;
z_alpha=probit(1-&alpha/2);
z_beta=probit(&power);
Ntot_pg1=2*2*(z_alpha+z_beta)**2*sigma2tot/(delta**2);
* calculate VIFs using rhotilde;
vif2=1+(n1-1)*rho12; 
rhotilde23=rho23*(n1*rho12)/(1+(n1-1)*rho12);
vif3=vif2*(1+(n2-1)*rhotilde23);
rhotilde34=rho34*(n2*rhotilde23)/(1+(n2-1)*rhotilde23);
vif4=vif3*(1+(n3-1)*rhotilde34);
* calculate rho and cluster size according to the scenario ;
%IF &cohort=4 %THEN %DO; rho=rho12*rho23*rho34*n1*n2*n3/vif4; vif_ml=vif4;clussize=n1*n2*n3; %END;
%IF &cohort=43 %THEN %DO; rho=rho12*rho23*n1*n2*(1+(n3-1)*rho34)/vif4; vif_ml=vif4;clussize=n1*n2*n3;%END; 
%IF &cohort=432 %THEN %DO; rho=1 - (1-rho12)/vif4;  vif_ml=vif4;clussize=n1*n2*n3;%END;
%IF &cohort=3 %THEN %DO; rho=rho12*rho23*n1*n2/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=32 %THEN %DO; rho=rho12*n1*(1+(n2-1)*rho23)/vif3; vif_ml=vif3;clussize=n1*n2;%END;
%IF &cohort=2 %THEN %DO; rho=n1*rho12/vif2; vif_ml=vif2;clussize=n1;%END;
if rho=0 then warning="rho=0, so standard error=0";
vif_swPG1=3/2*(1-rho)/(s-1/s)*(1+s*rho)/(1+s*rho/2);
vif_sw=vif_ml*vif_swPG1;
Ntot=Ntot_pg1*vif_sw;
n_clusters_raw=Ntot/clussize; * raw number of clusters;
n_clusters=max(ceil(n_clusters_raw),s); *number of cluster has to be integer and minimal one cluster per step;
*calculating power based on estimated number of clusters;
var_delta=4*sigma2tot*vif_sw/(clussize*n_clusters);
var_delta_raw=4*sigma2tot*vif_sw/(clussize*n_clusters_raw);
power_check=probnorm(abs(delta)/sqrt(var_delta)-z_alpha);
*power_check_raw=probnorm(abs(delta)/sqrt(var_delta_raw)-z_alpha);
run;
proc print data=_sw noobs; 
var cohort s n1 n2 n3 rho12 rho23 rho34 lambda1 lambda0 var_method sigma2tot delta Ntot_pg1 vif_ml rho vif_swPG1 vif_sw Ntot n_clusters power_aim power_check;
run;
%mend;

***********************EXAMPLES ****************************;



title "Example 1 (binary endpoint): power";
title2 "total variance provided via user_input";
%SWpower_bin(steps=4, cohort=43, n1=5,n2=15,n3=5,n4=4,rho12=0.6,rho23=0.05,rho34=0.01,p1=0.2, 
p0=0.35,p0_closest05=0.40,p1_closest05=0.25,sigma2tot_userinput=(0.25*(1-0.25)+0.40*(1-0.40))/(2*(1-0.6)), alpha=0.05);

title "reversed Example 2 (rate outcome): checking power based on sample size calculated in Example 2";
title2 "total variance based on lambda1 and lambda0"; 
%SWpower_rate(steps=4, cohort=32, n1=10,n2=4,n3=116,n4=1,rho12=0.7,rho23=0.01,lambda1=0.005, 
lambda0=0.011,lambda1_max=0.007,lambda0_max=0.013,sigma2tot_userinput=(0.008+0.020)/(2*(1-0.7)), alpha=0.05);

title "reversed Example 1 (binary endpoint): calculating number of clusters size based on the power as calculated in Example 1";
title2 "total variance provided via user_input";
%SWnclusters_bin(steps=4, cohort=43, n1=5,n2=15,n3=5,rho12=0.6,rho23=0.05,rho34=0.01,p1=0.20, p0=0.35, 
sigma2tot_userinput=(0.25*(1-0.25)+0.40*(1-0.40))/(2*(1-0.6)),alpha=0.05, power=0.80);

title "Example 2 (rate outcome): calculating number of clusters ";
title2 "total variance based on lambda1 and lambda0";
%SWnclusters_rate(steps=4, cohort=32, n1=10,n2=4,rho12=0.7,rho23=0.01,rho34=1,lambda1=0.005, lambda0=0.011,
alpha=0.05, power=0.80);



