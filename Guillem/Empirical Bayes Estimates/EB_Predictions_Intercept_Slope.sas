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

/* Empirical Bayes Predictions */

/* Empirical Bayes Predictions are the random effects predictions b_i */

/* GLMM (QUAD - Adaptive Gaussian Quadrature)*/
proc glimmix data=alzheimer_long method=QUAD(QPOINTS=10);
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
          TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET0 TIME*TAUPET0 TIME*SEX
          / dist=binary link=logit solution;
          
    random intercept TIME/ subject=PATID TYPE=UN solution;
    
	ods output SolutionR=eb_predictions;
run;

/* Create dataset for Random Intercepts Only */
data eb_intercepts;
    set eb_predictions;
    where Effect = 'Intercept'; 
run;

/* Create dataset for Random Slopes Only */
data eb_slopes;
    set eb_predictions;
    where Effect = 'TIME';      
run;

/* ----------------------------------------------------------- */
/* RANDOM EFFECTS ON THE INTERCEPT */
/* ----------------------------------------------------------- */

proc sort data=eb_intercepts;
    by Estimate;
run;

data plot_intercepts;
    set eb_intercepts;
    Patient_Rank = _N_; 
run;


/* INTERCEPT S-CURVE */

proc sgplot data=plot_intercepts;
    title "1. Random Intercepts: Subject Heterogeneity (Baseline Risk)";
    

    scatter x=Patient_Rank y=Estimate / markerattrs=(symbol=circlefilled size=2 color=blue);
    

    refline 0 / axis=y lineattrs=(color=black pattern=dash);
    

    yaxis label="Deviation from Average Baseline Risk (Log-Odds)";
    xaxis label="Patients (Sorted from Low Risk to High Risk)";
    

    inset "Above 0 = Starts WORSE than average" / position=topleft textattrs=(size=9 color=blue);
    inset "Below 0 = Starts BETTER than average" / position=bottomright textattrs=(size=9 color=blue);
run;

/* ----------------------------------------------------------- */
/* HISTOGRAM */
/* ----------------------------------------------------------- */
proc sgplot data=eb_intercepts;
    title "2. Distribution of Random Intercepts (Normality Check)";
    
    histogram Estimate / nbins=30 fillattrs=(color=lightblue);
    density Estimate / type=normal lineattrs=(color=darkblue thickness=2);
    
    yaxis grid;
    xaxis label="Random Intercept Estimate";
run;


/* ----------------------------------------------------------- */
/* RANDOM EFFECTS ON THE SLOPE */
/* ----------------------------------------------------------- */
proc sort data=eb_slopes;
    by Estimate;
run;

data plot_slopes;
    set eb_slopes;
    Patient_Rank = _N_; 
run;

/* ----------------------------------------------------------- */
/* SLOPE S-CURVE */
/* ----------------------------------------------------------- */
proc sgplot data=plot_slopes;
    title "3. Random Slopes: Subject Heterogeneity (Speed of Decline)";
    
    scatter x=Patient_Rank y=Estimate / markerattrs=(symbol=circlefilled size=2 color=darkred);
    
    refline 0 / axis=y lineattrs=(color=black pattern=dash);
    
    yaxis label="Deviation from Average Slope (Speed)";
    xaxis label="Patients (Sorted from Slowest to Fastest)";
    
    inset "Above 0 = Declining FASTER than average" / position=topleft textattrs=(size=9 color=darkred);
    inset "Below 0 = Declining SLOWER (More Stable)" / position=bottomright textattrs=(size=9 color=darkred);
run;

/* ----------------------------------------------------------- */
/* HISTOGRAM */
/* ----------------------------------------------------------- */
proc sgplot data=eb_slopes;
    title "4. Distribution of Random Slopes (Normality Check)";
    
    histogram Estimate / nbins=30 fillattrs=(color=pink);
    density Estimate / type=normal lineattrs=(color=darkred thickness=2);
    
    yaxis grid;
    xaxis label="Random Slope Estimate";
run;

