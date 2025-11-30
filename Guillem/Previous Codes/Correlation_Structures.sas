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


/* GEE, Generalized Estimating Equations */


/* Unstrucutred */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME AGE_STD BMI_STD ABPET TAUPET /
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;


/* Exchangeable */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME AGE_STD BMI_STD ABPET TAUPET /
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=EXCH covb corrw modelse;
run;


/* AR(1) */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME AGE_STD BMI_STD ABPET TAUPET /
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=AR covb corrw modelse;
run;


/* Independence */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME AGE_STD BMI_STD ABPET TAUPET /
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=IND covb corrw modelse;
run;


/* ALR, Alternating Logistic Regression */


proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX/
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=exch logor = exch covb corrw modelse;
run;


/* Linearization Based Method */


/* Unstrucutred */
proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX TIMECLSS;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=UN;
run;


/* Exchangeable */
proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX TIMECLSS;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=CS; /* Compound Symmetry */
run;

/* AR(1) */
proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX TIMECLSS;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=AR(1);
run;

/* Independence */
proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX TIMECLSS;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=VC; /* Variance Components */
run;










