/* some non-standard stepped wedge designs */

*************** closed cohort design *************************;

/* checked against my formulas in appendix ML-SW paper;


ods rtf file="MoreelBeraad_st200324.doc";
title3 "power for closed cohort";
data a;
do ICC=0.01, 0.05, 0.06;*intra cluster correlation;
do r_c=0.3, 0.5,0.8;
do r_s=0.5,0.8;
do I=6; *number of clusters;
do T=5; * number of time points;
do n=30;*number of subjects per cluster;
do delta=0.3,0.4, 0.45,0.5;
do SD=1;
theta=delta;
n_tot=I*n;* total number of individuals;
DE_c=1+(n-1)*ICC;
var_clustermean=DE_c*SD**2/n; 
r__=r_c*(n*ICC)/DE_c + r_s*(1-ICC)/DE_C; *correlation between cluster means;
* repeated measures design;
S=15; * sum of 1's;
R=43; *squared rowsums;
C=71; * squared columsums;
DE_rm_num= I*(1-r__)*(1+(T-1)*r__);
DE_rm_denom= S*I- C + r__*(S**2 +S*I*(T-1)-C*(T-1) - I*R) ;
DE_rm=DE_rm_num/DE_rm_denom;* Appendix ML SW paper p.31.;
var_theta=DE_rm*var_clustermean; *idem;
*see Appendix p.32 in multilevel SW paper, to reconcile Hoopers formula 11;
* that has a I**2/4 wrt the sample size of a t-test and ;
* mine that has a I wrt to the cluster mean;
* cluster stepped wedge;
	se_sw=sqrt(var_theta);
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=delta/sd; *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de_c) ); 
	output;
end;end; end;end;end;end;end;end;
run;

title2 "continuous";
title3 "T=5 with baseline, minimal detectable difference seems to be "; 
proc sort data=a; by I T N delta;
proc print data=a noobs;by I T N delta;  var ICC r_c r_s es theta se_sw   power_sw_df ;
;run;
**/


/*** cross-sectional designs ***********************/


/* design;
* T=3 (the baseline measurement has then to be obtained from retrospective data);
0 1 1 
0 1 1
0 1 1
0 0 1
0 0 1
*****
* T=2  (stepped wedge without baseline data), almost no loss of power;
1 1 
1 1
1 1
0 1
0 1
*/


title "primary endpoint: binary";
data a;
do I=5; 
do N=20,30,40,50,60; *number of subjects per cluster;
do T=2,3; * total number of measurements;
do rho = 0.1,0.15, 0.20;
do mu =0.2,0.3,0.4,0.5, 0.6;
do mu_r= 0,0.1;	
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * mu*(1-mu);
	theta=mu-mu_r; * the effect;
	U=8; * sum of 1's;
    V= 14; *squared rowsums;
	W=  34; * squared columsums;
	
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	se_theta=sqrt(var_theta);
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=2*abs(arsin(sqrt(mu)) - arsin(sqrt(mu_r))); *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

title2 "no baseline, one two periods";
proc print data=a noobs;by N;  var I T rho mu mu_r es theta var_theta se_theta   power_sw_df ;
where T=2 and theta >=0.3;*power_par_df power_sw power_par ;
 ;run;


title2 "including baseline, three periods"; 
proc print data=a noobs;by N;  var I T rho mu mu_r es theta var_theta se_theta   power_sw_df ;
where T=3 and theta >=0.3;*power_par_df power_sw power_par ;
 ;run;


*********************  QoL **************************************;
title "QoL : continuous";
data a;
do I=5; 
do N=20,30,40,50,60; *number of subjects per cluster;
do T=2,3; * total number of measurements;
do rho = 0.1,0.15, 0.20;
do delta =12;
do SD= 25.6;	
	sigma2=(SD)**2 / N;
	tau2=( rho/(1-rho) ) * SD**2;
	theta=delta; * the effect;
	U=8; * sum of 1's;
    V= 14; *squared rowsums;
	W=  34; * squared columsums;
	
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	se_theta=sqrt(var_theta);
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=delta/sd; *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

title2 "including baseline, three periods"; 
proc print data=a noobs;by N;  var I T rho es theta var_theta se_theta   power_sw_df power_par_df;
where T=3 and theta >=0.3;*power_par_df power_sw power_par ;
 ;run;


*********************************************************************************;
** K before measurements, L after measurements;
*** 0 ... 0 1 1 1 ...1 1....1
*** 0 ....0 0 1 1 ...1 1....1
*** 0 ....0 0 0 1 ...1 1....1
*** ....
*** 0 ....0 0 0 ...0 1 1....1
*** K meas  step.wedge  L meas;

* title "binary outcome ";  
title2 "tau2= (rho/(1-rho) )* sigma2_e ";
title3 " I=number of clusters, N=cluster size per time unit";
title4 " standard stepped wedge with K-1 leading control measurments and L post measurements";
title6 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=5; 
do K=6;
do L=6;
do N=10; *number of subjects per cluster;
do rho = 0.1;
do mu =0.2,0.3,0.4,0.5, 0.6;
do mu_r= 0;	
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * mu*(1-mu);
	theta=mu-mu_r; * the effect;
	U=I*(I+1)/2 + I*L; * sum of 1's;
    V= (I+L)**3/3 + (I+L)**2/2 + (I+L)/6 
	   - L**3/3  - L**2/2 - L/6; *squared rowsums;
	W=  I**3/3 + I**2/2 + I/6 + L*I**2; * squared columsums;
	T=K+I+L; * total number of measurements;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	se_theta=sqrt(var_theta);
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=2*abs(arsin(sqrt(mu)) - arsin(sqrt(mu_r))); *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;end;
run;

proc print data=a noobs;  var I K L T N rho mu mu_r es theta var_theta se_theta   power_sw_df power_par_df;*power_sw power_par ;
 ;run;

/*
 *** check above with standard stepped wedge (use K=1, L=0 in  the above);
 title "binary outcome (CORRECTED calculation for somatic units)";  
title2 "tau2= (rho/(1-rho) )* sigma2_e ";
title3 " I=number of clusters, T=number of timepoints, N=cluster size";
title4 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=6; 
do T=7;
do N=10; *number of subjects per cluster;
do	rho = 0.1;
do mu =0.50;
do mu_r= 0.30;	
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * mu*(1-mu);
	theta=mu-mu_r; * the effect;
	U=(I*T)/2;
	W= (I**2 * T *(2*T-1) ) / (6*(T-1)); 
	V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=2*abs(arsin(sqrt(mu)) - arsin(sqrt(mu_r))); *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

proc print data=a noobs; var I T N rho mu mu_r  U V W theta es power_sw power_sw_df power_par power_par_df;
 ;run;
*/

 ******************************************************************************;
 /* design
 with baseline measurements T=5;
 0	1 	1	1	1
 0	0	1	1	1
 0	0	1	1	1
 0	0	0	1	1
 0	0	0	1	1
 0	0	0	0	1
 or without the baseline measurements T=4

 */

data a;
do I=6; 
do N=30; *number of subjects per cluster;
do T=5; * total number of measurements;
do rho = 0.1;
do delta =0.2,0.3,0.4,0.42,0.43,0.45,0.5;
do SD= 1;	
	sigma2=(SD)**2 / N;
	tau2=( rho/(1-rho) ) * SD**2;
	theta=delta; * the effect;
	U=15; * sum of 1's;
    V=43; *squared rowsums;
	W=71; * squared columsums;
	
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	se_theta=sqrt(var_theta);
	power_sw=probnorm(-1.96 + abs(theta) / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + abs(theta) / sqrt(var_theta), I );
	* parallel cluster;
	es=delta/sd; *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

title2 "T=5 with baseline, minimal detectable difference seems to be 0.42*SD"; 
proc print data=a noobs;by N;  var I T rho es theta var_theta se_theta   power_sw_df power_par_df;
;run;

