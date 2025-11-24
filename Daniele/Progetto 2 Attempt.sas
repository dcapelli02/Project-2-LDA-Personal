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

    keep PATID SEX AGE AGE_STD BMI BMI_STD JOB ADL 
         TIME CDRSB_CAT ABPET0 TAUPET0;
run;

/* Here we need to create a new variable for time */

proc genmod data=alzheimer_long;
    class PATID SEX;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX/
        dist=binomial link=logit;

    repeated subject=PATID / type=unstructured corrw;
run;

/* NOTA: con il corrw faccio gee2? */

/* ALTERNATING LOGISTIC REGRESSION */
proc genmod data=alzheimer_long;
    class PATID SEX;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX/
        dist=binomial link=logit;

    repeated subject=PATID / type=exch logor = exch;
run;

/* LINEARIZATION BASED METHOD */
/* Commented because a little bit long to compute */
/*
proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=un;
run;
*/

/* TO DO: test all of them */



/* Prova per modello con random effects */
proc glimmix data=alzheimer_long method=quad(q=15) adaptive;
	class PATID SEX;
	model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random intercept / subject=PATID;
run;

/* Controllare se codice a posto */

/* Then test all other methods (PQL...) */

/* Still need to figure out procedure NLMIXED */

/* Graphical analysis? */

/* Empirical Bayes estimates? --> nlmixed? */
