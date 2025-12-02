data alzheimer25;
    set "/home/u64400254/sasuser.v94/alzheimer25.sas7bdat"; 
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
         and not missing(BPRS)
         and not missing(ABPET)
         and not missing(TAUPET) then do;

            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;

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
    var ABPET TAUPET; 
    output out=bio_stats mean=m_ab m_tau std=s_ab s_tau;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
run;

/* GEE (Quadratic) */

proc genmod data = alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

	model CDRSB_CAT = TIME TIME * TIME SEX BMI_STD ADL 
     TIME * BMI_STD TIME * ADL TIME * TIME * ADL /

        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
    output out = gee_preds p = Predicted_Prob;
run;

proc means data=gee_preds nway noprint;
    class TIME;
    var CDRSB_CAT Predicted_Prob;
    output out=plot_data mean=Observed_Prop Mean_Predicted;
run;

/* Plot */
proc sgplot data=plot_data;
    title "Observed vs. Predicted Risk over Time (GEE)";
    

    series x=TIME y=Observed_Prop / lineattrs=(color=black thickness=2) 
           markerattrs=(symbol=circlefilled color=black) 
           legendlabel="Observed Proportion";
           

    series x=TIME y=Mean_Predicted / lineattrs=(color=red pattern=dash thickness=2) 
           markerattrs=(symbol=none) 
           legendlabel="Predicted Probability (GEE)";
           
    yaxis label="Probability of Severe Dementia (CDRSB >= 10)" grid values=(0.2 to 0.7 by 0.1);
    xaxis label="Years (Time)";
run;



