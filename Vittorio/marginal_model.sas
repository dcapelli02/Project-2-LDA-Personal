data alzheimer_long;
    if _n_ = 1 then set stats;   
    set alzheimer25;

    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6]  bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;

    /* baseline da usare come covariate fisse */
    retain BPRS_BASE ABPET_BASE TAUPET_BASE;

    do TIME = 0 to 6;

        /* risposta longitudinale */
        CDRSB = cdrsb_arr[TIME];

        /* baseline solo al tempo 0 */
        if TIME = 0 then do;
            BPRS_BASE  = bprs_arr[0];
            ABPET_BASE = abpet_arr[0];
            TAUPET_BASE= taupet_arr[0];
        end;

        /* dichotomize CDRSB, includendo NA come categoria separata */
        if CDRSB < 10 then CDRSB_CAT = 0; 
        else CDRSB_CAT = 1;

        /* standardizza age e bmi */
        AGE_STD = (AGE - mean_age) / sd_age;
        BMI_STD = (BMI - mean_bmi) / sd_bmi;

        TIMECLSS = put(TIME, 1.);

        output;
    end;

    keep TRIAL PATID SEX AGE AGE_STD EDU BMI BMI_STD INKOMEN JOB ADL WZC 
         TIME TIMECLSS CDRSB_CAT 
         BPRS_BASE ABPET_BASE TAUPET_BASE;
run;






/* --- EXCHANGEABLE --- */
ods output GEEEmpPEst=gee_exch_raw;
proc genmod data=alzheimer_long descending;
  class PATID SEX TIMECLSS;
  model CDRSB_CAT = SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
                    TIME 
                    SEX*TIME AGE_STD*TIME BMI_STD*TIME
                    ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
        / dist=binomial link=logit;
  repeated subject=PATID / withinsubject=timeclss type=exch covb corrw modelse;
run;


/* --- INDEPENDENCE --- */
ods output GEEEmpPEst=gee_ind_raw;
proc genmod data=alzheimer_long descending;
  class PATID SEX TIMECLSS;
  model CDRSB_CAT = SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
                    TIME 
                    SEX*TIME AGE_STD*TIME BMI_STD*TIME
                    ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
        / dist=binomial link=logit;
  repeated subject=PATID / withinsubject=timeclss type=ind covb corrw modelse;
run;


/* --- UNSTRUCTURED --- */
ods output GEEEmpPEst=gee_un_raw;
proc genmod data=alzheimer_long descending;
  class PATID SEX TIMECLSS;
  model CDRSB_CAT = SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
                    TIME 
                    SEX*TIME AGE_STD*TIME BMI_STD*TIME
                    ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
        / dist=binomial link=logit;
  repeated subject=PATID / withinsubject=timeclss type=un covb corrw modelse;
run;

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */


/* --- Alternating Logistic Regression EXCHANGEABLE --- */
ods output GEEEmpPEst=gee_exch_raw;
proc genmod data=alzheimer_long descending;
  class PATID SEX TIMECLSS;
  model CDRSB_CAT = SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
                    TIME 
                    SEX*TIME AGE_STD*TIME BMI_STD*TIME
                    ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
        / dist=binomial link=logit;
  repeated subject=PATID / withinsubject=timeclss logor=exch covb corrw modelse;
run;

/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/* -/* --- EXCHANGEABLE --- */
ods output ParameterEstimates = glimmix_exch_raw;

proc glimmix data=alzheimer_long method=RSPL empirical;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT(event='1') =
            SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
            TIME
            SEX*TIME AGE_STD*TIME BMI_STD*TIME
            ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
            / dist=binary link=logit solution;
    random _residual_ / subject=PATID type=cs;
run;


/* --- Linearization Based Method EXCHANGEABLE --- */
ods output ParameterEstimates = glimmix_exch_raw;

proc glimmix data=alzheimer_long method=RSPL empirical;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT(event='1') =
            SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
            TIME
            SEX*TIME AGE_STD*TIME BMI_STD*TIME
            ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
            / dist=binary link=logit solution;
    random _residual_ / subject=PATID type=cs;
run;

/* --- Linearization Based Method INDEPENDENCE --- */
ods output ParameterEstimates = glimmix_ind_raw;

proc glimmix data=alzheimer_long method=RSPL empirical;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT(event='1') =
            SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
            TIME
            SEX*TIME AGE_STD*TIME BMI_STD*TIME
            ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
            / dist=binary link=logit solution;
    random _residual_ / subject=PATID type=simple;
run;

/* --- Linearization Based Method UNSTRUCTURED --- */
ods output ParameterEstimates = glimmix_un_raw;

proc glimmix data=alzheimer_long method=RSPL empirical;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT(event='1') =
            SEX AGE_STD BMI_STD ADL ABPET_BASE TAUPET_BASE
            TIME
            SEX*TIME AGE_STD*TIME BMI_STD*TIME
            ADL*TIME ABPET_BASE*TIME TAUPET_BASE*TIME
            / dist=binary link=logit solution;
    random _residual_ / subject=PATID type=un;
run;











