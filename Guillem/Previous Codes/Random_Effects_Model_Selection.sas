data alzheimer25;
    set "/home/u64400254/sasuser.v94/alzheimer25.sas7bdat"; 
run;

/* Calculates mean, standard deviation for patients' Age and BMI */

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

        	TIMECLSS = put(TIME, 1.);


            output;
        end;
    end;

	keep PATID SEX AGE AGE_STD BMI BMI_STD ADL 
     	TIME TIMECLSS CDRSB CDRSB_CAT ABPET TAUPET ABPET0 TAUPET0;
run;

/* Random-Effects Model */

/* GLMM (Generalized Linear Mixed Model)


/* GLMM (MQL - Marginal Quasi-Likelihood) */
proc glimmix data=alzheimer_long method=RMPL;
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
          TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET0 TIME*TAUPET0 TIME*SEX
          / dist=binary link=logit solution;
          
    random intercept / subject=PATID;
run;
/* GLMM (PQL - Penalized Quasi-Likelihood) */
proc glimmix data=alzheimer_long method=RSPL;
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
          TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET0 TIME*TAUPET0 TIME*SEX
          / dist=binary link=logit solution;
          
    random intercept / subject=PATID;
run;

/* GLMM (QUAD - Adaptive Gaussian Quadrature)*/
proc glimmix data=alzheimer_long method=QUAD(QPOINTS=10);
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
          TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET0 TIME*TAUPET0 TIME*SEX
          / dist=binary link=logit solution;
          
    random intercept / subject=PATID;
run;

%put The coefficients almost look the same so we will choose the GLMM (QUAD) for simplicity.;
%put And also because it is theoretically the best model.;

/* Beta-Binomial Model */


data alz_betabin;
    if _n_ = 1 then set stats; 
    set alzheimer25;

    AGE_STD = (AGE - mean_age) / sd_age;
    BMI_STD = (BMI - mean_bmi) / sd_bmi;

    array raw[0:6] cdrsb0-cdrsb6;
    Y_SUM = 0;
    N_TRIALS = 0;
    
    do i = 0 to 6;
        if raw[i] ne . then do;
            N_TRIALS = N_TRIALS + 1;
            if raw[i] >= 10 then Y_SUM = Y_SUM + 1;
        end;
    end;

    if N_TRIALS > 0;
    
    keep PATID Y_SUM N_TRIALS SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0;
run;

proc fmm data=alz_betabin;
    class SEX;

    model Y_SUM / N_TRIALS = SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0
          / dist=betabinomial;
          
run;

%put Since Scale Parameter = 0, the data behaves like a standard Binomial;


