data alzheimer25;
	set '~/Project 1 LDA/alzheimer25.sas7bdat';
run;


proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;


data alzheimer_long;
    if _n_ = 1 then set stats;   
    set alzheimer25;
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;

    do TIME = 0 to 6;
        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];

        if not missing(CDRSB) 
         or not missing(BPRS)
         or not missing(ABPET)
         or not missing(TAUPET) then do;

            /* dichotomize cdrsb */
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;

            /* standardize age and bmi */
            AGE_STD = (AGE - mean_age) / sd_age;
            BMI_STD = (BMI - mean_bmi) / sd_bmi;

            output;
        end;
    end;

    keep PATID SEX AGE AGE_STD BMI BMI_STD JOB ADL ABPET TAUPET
         TIME CDRSB_CAT ABPET0 TAUPET0;
run;

/* Overall mean trajectory of BPRS over time */
proc sgplot data=alzheimer_long;
    title "Mean CDRSB Trajectory Over Time";
    vline TIME / response=CDRSB_CAT stat=mean markers 
                 limitstat=stderr limits=both;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean CDRSB";
run;



/* RANDOM EFFECTS */

/* PQL */
proc glimmix data=alzheimer_long method=rspl;
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* MQL */
proc glimmix data=alzheimer_long method=rmpl;
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Risultati molto simili */

/* Laplace */
proc glimmix data=alzheimer_long method=laplace;
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Adaptive Gaussian */
proc glimmix data=alzheimer_long method=quad(qpoints = 30);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Comparison with different q
I think 30 is enough but show it
*/

proc glimmix data=alzheimer_long method=quad(qpoints = 3);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

proc glimmix data=alzheimer_long method=quad(qpoints = 50);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;


/* REDUCTION
We can have a look at the table with the Wald tests:
I would get rid of ADL and TIME*ADL and TIME*SEX
The others are more or less significant
(AB not much but it is interesting)
*/

/* Adaptive Gaussian */
proc glimmix data=alzheimer_long method=quad(qpoints = 30);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ABPET TIME * TAUPET
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Compare BIC/QIC --> take second model */

/* Just out of curiosity:
should we get rid of ABPET and TIME*ABPET?
*/

/* Adaptive Gaussian */
proc glimmix data=alzheimer_long method=quad(qpoints = 30);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX AGE_STD BMI_STD TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * TAUPET
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Yes but I would keep them just to have
a comparison also with them */



/* Model con nlmixed */

proc nlmixed data=alzheimer_long qpoints=30;
parms beta0=0.4935 beta1=0.07764 beta2=0.1167 beta3=0.1172 beta4 = -0.1710 beta5 = -0.6089
	  beta6 = -0.4547 beta7 = -0.07627 beta8 = 0.08884 beta9 = 0.2529
	  sigma=0.044;
	teta = beta0 + b + beta1*SEX + beta2*AGE_STD +
			beta3*BMI_STD + beta4*ABPET + beta5*TAUPET
			+ beta6*TIME + beta7*TIME * BMI_STD + beta8 * TIME * ABPET + beta9 * TIME * TAUPET;
	expteta = exp(teta);
	p = expteta/(1+expteta);
model CDRSB_CAT ~ binary(p);
random b ~ normal(0,sigma**2) subject=PATID;
predict b out=EB_estimates;   /* <-- qui ottieni le EB */
run;

/* Graphical analysis? */
/* --> see R */

/* Empirical Bayes estimates? --> nlmixed? */

proc sgplot data=EB_estimates;
   histogram pred;   /* "pred" è la variabile che contiene le EB stimate */
   density pred;     /* opzionale: aggiunge la curva di densità normale stimata */
   xaxis label="Stime EB (random effect b)";
   yaxis label="Frequenza";
run;



/* PROVA CON PROBIT */
/* I think we can ignore this one */
/*
proc nlmixed data=alzheimer_long qpoints=30;
   parms beta0=1.2 beta1=-0.2 beta2=-0.12 beta3=-0.3 beta4=0.05 beta5=0.05 sigma=0.02;
   eta = beta0 + b 
         + beta1*SEX 
         + beta2*BMI_STD 
         + beta3*TIME 
         + beta4*TIME*SEX 
         + beta5*TIME*BMI_STD;
   p = probnorm(eta); 
   model CDRSB_CAT ~ binary(p);
   random b ~ normal(0, sigma**2) subject=PATID;
   predict b out=EB_estimates;
run;
*/


/* Random Slope */
/* I think we can ignore this one */
/*
proc glimmix data=alzheimer_long method=quad(qpoints = 30);
	class PATID SEX;
	model CDRSB_CAT(event = '1') = TIME SEX BMI_STD TIME*BMI_STD TIME*SEX
		/ dist=binary solution;
	random intercept TIME / subject=PATID;
run;
*/
