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

/* GEE */


/* We chose GEE because it is the most simple one. Now we will reduce it.*/



/* Doing step-by-step procedure (eliminating variables with p-value > 0.05) I ended up with this: */

proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME BMI_STD AGE_STD SEX
     TIME * BMI_STD TIME/
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;

/* ???? all variables are gone xd */
/* So since ABPET, TAUPET and ADL are highly correlated with age, that's why they disappeared? */






