
data alzheimer25;
	set '/home/u63619091/LDA HWs/HW1/alzheimer25.sas7bdat';
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

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET; /* Calculate mean/sd for biomarkers */
    output out=bio_stats mean=m_ab m_tau std=s_ab s_tau;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    /* Create Standardized Versions */
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
run;

data alzheimer_long_centered;
    set alzheimer_long_centered;
    /* add log-transformed time; handle zeros if needed */
    LOGTIME = log(TIME);
	LOGAGE = log(AGE_STD);
run;


/* The following models do currently implement log(age), they use standardized age. */
/* fit may improve if we switch to log(age) */

proc genmod data=alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET_STD TAUPET_STD
     TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET_STD TIME*TAUPET_STD
     TIME*SEX/
        dist=binomial link=logit type3;

    repeated subject=PATID / withinsubject=timeclss type=exch logor=FULLCLUST covb corrw modelse;
run;

proc genmod data=alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET_STD TAUPET_STD
     TIME*AGE_STD TIME*BMI_STD TIME*TIME*TIME TIME*TIME*TIME*TAUPET_STD TIME*TIME
     TIME*SEX/
        dist=binomial link=logit type3;

    repeated subject=PATID / withinsubject=timeclss type=exch logor=FULLCLUST covb corrw modelse;
run;

data alzheimer_long_centered;
    set alzheimer_long_centered;
    /* add log-transformed time; handle zeros if needed */
    LOGTIME = log(TIME);
run;

proc genmod data=alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = LOGTIME SEX AGE_STD BMI_STD ADL ABPET_STD TAUPET_STD
                       LOGTIME*AGE_STD LOGTIME*BMI_STD LOGTIME*TAUPET_STD
/*                        LOGTIME*LOGTIME*LOGTIME */
/*                        LOGTIME*LOGTIME*LOGTIME*TAUPET_STD */
                       LOGTIME*LOGTIME
/*                        LOGTIME*SEX */
        / dist=binomial link=logit type3;

    repeated subject=PATID /
        withinsubject=TIMECLSS
        type=exch
        logor=FULLCLUST
        covb corrw modelse;
run;



/* Return to marginal model */

proc genmod data=alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

    /* Model: population-averaged logistic regression */
    model CDRSB_CAT = LOGTIME
                       SEX
                       AGE_STD
                       BMI_STD
                       ADL
                       ABPET_STD
                       TAUPET_STD
                       LOGTIME*AGE_STD
                       LOGTIME*BMI_STD
                       LOGTIME*SEX
                       LOGTIME*ADL
        / dist=binomial link=logit type3;

    /* Specify repeated measures (GEE) */
    repeated subject=PATID / type=exch corrw covb modelse;
run;

/* Random Effects Model */

proc glimmix data=alzheimer_long_centered method=laplace;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT(event='1') =
            LOGTIME
            SEX
            AGE_STD
            BMI_STD
            ADL
            ABPET_STD
            TAUPET_STD
            LOGTIME*AGE_STD
            LOGTIME*BMI_STD
            LOGTIME*SEX
            LOGTIME*ADL
            LOGTIME*ABPET_STD
            LOGTIME*TAUPET_STD
          
        / dist=binomial link=logit solution;

    /* Random effects: intercept + slope for logtime */
    random intercept / subject=PATID type=un;
    covtest;
run;

proc glimmix data=alzheimer_long_centered method=laplace;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT(event='1') =
            LOGTIME
            SEX
            AGE_STD
            BMI_STD
            ADL
            ABPET_STD
            TAUPET_STD
            LOGTIME*AGE_STD
            LOGTIME*BMI_STD
            LOGTIME*SEX
            LOGTIME*ADL
        / dist=binomial link=logit solution;

    random LOGTIME / subject=PATID type=un;

    covtest;
run;
