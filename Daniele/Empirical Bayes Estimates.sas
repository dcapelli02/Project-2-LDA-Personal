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
    
    /* Create Variables */
   TIME_QUAD = TIME*TIME;
   TIME_SEX = TIME*SEX;
   TIME_BMI = TIME*BMI_STD;
   TIME_AGE = TIME*AGE_STD;
   TIME_ADL = TIME*ADL;
   TIME_AB = TIME*ABPET_STD;
   TIME_TAU = TIME*TAUPET_STD;
   
   TIME_SEX_QUAD = TIME_QUAD * SEX;
   TIME_BMI_QUAD = TIME_QUAD * BMI_STD;
   TIME_AGE_QUAD = TIME_QUAD * AGE_STD;
   TIME_ADL_QUAD = TIME_QUAD * ADL;
   TIME_AB_QUAD = TIME_QUAD * ABPET_STD;
   TIME_TAU_QUAD = TIME_QUAD * TAUPET_STD;
run;


/* ============================================================= */
/* 0. SALVATAGGIO DEI RANDOM EFFECTS (Empirical Bayes Estimates) */
/*    --> Usando _RESIDUAL_ invece di SolutionR                  */
/* ============================================================= */

ods select all; 
ods output ResidualPanel = EB_Estimates_Data;   /* <<< MODIFICA QUI */


/* ------------------------------------------------------------- */
/* 1. ESECUZIONE DEL MODELLO GLIMMIX                             */
/* ------------------------------------------------------------- */

proc glimmix data=alzheimer_long_centered method=quad(qpoints = 30);
class PATID SEX;

model CDRSB_CAT(event = '1') = TIME BMI_STD TAUPET_STD
                               TIME_BMI TIME_TAU
                               / dist=binary solution;

/* random intercept + slope */
random intercept TIME / subject=PATID type=un;

/* salvo predizioni */
output out=predictions_data
       pred(ilink)=P_FIXED
       predicted(ilink)=P_BLUP
       residual;     /* <<< IMPORTANTE PER _RESIDUAL_ */
run;

ods select default;


/* ============================================================= */
/* 2. PREPARAZIONE RANDOM EFFECTS                                */
/*    ResidualPanel ha una struttura diversa → sistemiamo        */
/* ============================================================= */

/* ResidualPanel contiene variabili:
   Effect, PATID, Estimate, _Residual_
   → noi usiamo _Residual_ come EB estimate
*/
data EB_Estimates_Data;
    set EB_Estimates_Data;
    Estimate = _Residual_; /* <<< MODIFICA */
    keep PATID Effect Estimate;
run;

/* Trasformazione per scatterplot */
proc transpose data=EB_Estimates_Data 
               out=EB_Estimates_Wide(drop=_name_ _label_);
    by PATID;
    id Effect;
    var Estimate;
run;


/* ============================================================= */
/* 3. GRAFICI                                                    */
/* ============================================================= */

/* 3a Intercetta */
proc sgplot data=EB_Estimates_Data;
    where Effect='Intercept';
    histogram Estimate / fillattrs=(color=lightblue);
    density Estimate / type=kernel lineattrs=(color=darkblue thickness=2);
    title "Distribuzione dei Random Intercept (Residual EB)";
run;

/* 3b Slope */
proc sgplot data=EB_Estimates_Data;
    where Effect='TIME';
    histogram Estimate / fillattrs=(color=lightgreen);
    density Estimate / type=kernel lineattrs=(color=darkgreen thickness=2);
    title "Distribuzione dei Random Slopes (Residual EB)";
run;

/* 3c Correlazione */
proc sgplot data=EB_Estimates_Wide;
    scatter x=Intercept y=TIME /
            markerattrs=(symbol=circlefilled size=7 color=maroon);
    reg x=Intercept y=TIME /
        lineattrs=(color=black thickness=2);
    title "Correlazione Intercept-Slope (Residual EB)";
run;


/* ============================================================= */
/* 4. PLOT OSSERVATI VS PREDIZIONI (Fixed Effects)               */
/* ============================================================= */

proc means data=predictions_data nway noprint;
    class TIME;
    var CDRSB_CAT P_FIXED;
    output out=plot_data mean=Observed_Prop Mean_Predicted;
run;

proc sgplot data=plot_data;
    title "Observed vs Predicted Probability Over Time";

    series x=TIME y=Observed_Prop /
           lineattrs=(color=black thickness=2)
           markerattrs=(symbol=circlefilled);

    series x=TIME y=Mean_Predicted /
           lineattrs=(color=red pattern=dash thickness=2);

    yaxis label="Probability (CDRSB ≥ 10)";
    xaxis label="Years";
run;
