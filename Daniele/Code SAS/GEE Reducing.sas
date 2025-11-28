data alzheimer25;
	set '~/Project 1 LDA/alzheimer25.sas7bdat';
run;


proc means data=alzheimer25 noprint;
    var AGE BMI ABPET0 TAUPET0;
    output out=stats mean=mean_age mean_bmi mean_abpet0 mean_taupet0
    std=sd_age sd_bmi sd_abpet0 sd_taupet0;
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
            AB_STD = (ABPET0 - mean_abpet0) / sd_abpet0;
            TAU_STD = (TAUPET0 - mean_taupet0) / sd_taupet0;
            
            TIMECLSS = put(TIME, 1.);

            output;
        end;
    end;

    keep PATID SEX AGE AGE_STD BMI BMI_STD JOB ADL 
    		AB_STD TAU_STD ABPET TAUPET
         TIME CDRSB_CAT ABPET0 TAUPET0 TIMECLSS;
run;

/* Overall mean trajectory of BPRS over time */
proc sgplot data=alzheimer_long;
    title "Mean CDRSB Trajectory Over Time";
    vline TIME / response=CDRSB_CAT stat=mean markers 
                 limitstat=stderr limits=both;
    xaxis label="Time (years)" values=(0 to 6 by 1);
    yaxis label="Mean CDRSB";
run;


/* GEE */

/* Longitudinal */
/* Unstrucutred */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET TIME * TAUPET
     TIME * SEX /
        dist=binomial link=logit type3;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;

/* From type3 table get rid of time*sex, time*adl, time*ab
and see how it goes */

proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET TAUPET 
     TIME * AGE_STD TIME * BMI_STD TIME * TAUPET /
        dist=binomial link=logit type3;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;

/* Better QIC, ok */

/* Try to reduce also timeage ab sex adl */

proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME AGE_STD BMI_STD TAUPET 
     TIME * BMI_STD TIME * TAUPET /
        dist=binomial link=logit type3;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;

/* Final model, should be ok in terms of the QIC */

