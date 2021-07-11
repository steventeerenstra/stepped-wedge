title "binary outcome (CORRECTED calculation for somatic units)";  
title2 "tau2= (rho/(1-rho) )* sigma2_e ";
title3 " I=number of clusters, T=number of timepoints, N=cluster size";
title4 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=14,16,18,20; 
do T=6;
do N=20; *number of subjects per cluster;
do	rho = 0.1;
do mu =0.22;
do mu_r= 0.132;	
	* the below comes down to H&H choice, so sigma2_s=mu*(1-mu) (or mu=(mu+mu_r)/2), and sigma2_c=( rho/(1-rho) ) *sigma2_s;
	* alternative would be to take sigma2_tot=sigma2_c+sigma_2_s=mu*(1-mu) and sigma2_c= rho* sigma2_tot, sigma2_s=(1-rho)*sigma2_tot;  
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * mu*(1-mu);
	theta=mu-mu_r; * the effect;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
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

proc print data=a noobs; var I T N rho mu mu_r  theta es power_sw power_sw_df power_par power_par_df;
 ;run;

******** PG afdeling ************;
title "binary outcome (CORRECTED calculation for psychiatric geriatric units)";  
title2 "tau2= (rho/(1-rho) )* sigma2_e ";
title3 " I=number of clusters, T=number of timepoints, N=cluster size";
title4 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=14,16,18,20; 
do T=6;
do N=20; *number of subjects per cluster;
do	rho = 0.1;
do mu =0.30;
do mu_r= 0.195;	
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * mu*(1-mu);
	theta=mu-mu_r; * the effect;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + theta / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + theta / sqrt(var_theta), I );
	* parallel cluster;
	es=2*abs(arsin(sqrt(mu)) - arsin(sqrt(mu_r))); *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

proc print data=a noobs; var I T N rho mu mu_r  theta es power_sw power_sw_df power_par;
 ;run;





 ******** continuous outcome ************;
title "CORRECTED calculation for continuous outcome";  
title2 "tau2= (rho)*sigma_tot**2 ";
title3 " I=number of clusters, T=number of timepoints, N=cluster size";
title4 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=14,16,18,20; 
do T=6;
do N=20; *number of subjects per cluster;
do	rho = 0.1;
do theta=0.2; *effect= difference in means;
do sigma_tot=1.0; * total sd; 
	sigma2=(1-rho)*sigma_tot**2/N; * variance of a cluster mean due to sampling N subjects;
	tau2= rho* sigma_tot**2; *between cluster variance;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + theta / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + theta / sqrt(var_theta), I );
	* parallel cluster;
	es=theta/sigma_tot; *Cohen's effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + es/sqrt(2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + es/sqrt(2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

proc print data=a noobs; var I T N rho sigma_tot theta es power_sw power_sw_df power_par;
 ;run;
 ****;
/***
title "incorrect calculation to see the difference (for the somatic case)";
title2 "the definition tau2= (rho/(1-rho) )* sigma2_e is wrong ..";
title3 " I=number of clusters, T=number of timepoints, N=cluster size";
title4 "power_sw_df=power stepped wedge, power_par=power parallel cluster trial";
data a;
do I=14; 
do T=6;
do N=20; *number of subjects per cluster;
do	rho = 0.1;
do mu =0.22;
do mu_r= 0.132;	
	sigma2=mu*(1-mu) / N;
	tau2=( rho/(1-rho) ) * sigma2;
	theta=mu-mu_r; * the effect;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + theta / sqrt(var_theta) );
	*based on I = number of clusters degrees of freedom;
	power_sw_df=probt( tinv(0.025,I) + theta / sqrt(var_theta), I );
	* parallel cluster;
	es=2*abs(arsin(sqrt(mu)) - arsin(sqrt(mu_r))); *Cohen's h i.e. the effect size;
	de= (1+ (N-1)*rho);; * the design factor;
	power_par=probnorm(-1.96 + (es/2) * sqrt( (I*N)/de) ); 
	*based on I=number of clusters as degree of freedom;
	power_par_df=probnorm( tinv(0.025,I-2) + (es/2) * sqrt( (I*N)/de) ); 
	output;
end;end;end;end;end;end;
run;

proc print data=a noobs; var I T N rho mu mu_r  theta es power_sw power_sw_df power_par;
 ;run;
***/

 /*** example  
title "checking the example of Fig 2, Hussey, Contemp Clin Trial 2007 28: 182-191";
title2 "note that reading the figure suggests that 0.375 is a more likely estimate at a power of 0.80 ";
options nocenter formdlim="";

data b;
do I=24; 
do T=5;
do N=100; *number of subjects per cluster;
do	cv= 0.3, 0.5;
do rr=0.375;
	mu=0.05;
	tau2=(cv*mu)**2;
	sigma2=mu*(1-mu) / N;
	theta=rr*mu;
	U=(I*T)/2;W= (I**2 * T *(2*T-1) ) / (6*(T-1)); V= ( I*T*(2*T-1) ) / 6;
	var_theta= ( I*sigma2*(sigma2 + T*tau2) ) 
                    /
                    ( (I*U-W)*sigma2 + (U**2 +I*T*U - T*W - I*V)*tau2  );

	* cluster stepped wedge;
	power_sw=probnorm(-1.96 + theta / sqrt(var_theta) );
	output;
end;end;end;end;end;
run;

proc print data=b;run;
***/
