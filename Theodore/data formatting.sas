data alzheimer25;
	set '/home/u63619091/LDA HWs/HW1/alzheimer25.sas7bdat';
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

    keep TRIAL PATID SEX AGE AGE_STD EDU BMI BMI_STD INKOMEN JOB ADL WZC 
         TIME CDRSB CDRSB_CAT BPRS ABPET TAUPET;
run;

proc genmod data=alzheimer_long;
    class PATID SEX JOB WZC TRIAL;

    model CDRSB_CAT = TIME AGE_STD BMI_STD EDU BPRS ABPET TAUPET /
        dist=binomial link=logit;

    repeated subject=PATID / type=exchangeable;
run;
