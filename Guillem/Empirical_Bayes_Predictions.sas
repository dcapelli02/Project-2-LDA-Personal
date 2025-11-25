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

/* Empirical Bayes Predictions */

/* Empirical Bayes Predictions are the random effects predictions b_i */

/* GLMM (QUAD - Adaptive Gaussian Quadrature)*/
proc glimmix data=alzheimer_long method=QUAD(QPOINTS=10);
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
          TIME*AGE_STD TIME*BMI_STD TIME*ADL TIME*ABPET0 TIME*TAUPET0 TIME*SEX
          / dist=binary link=logit solution;
          
    random intercept / subject=PATID solution;
    
	ods output SolutionR=eb_predictions;
run;

proc print data=eb_predictions (obs=10);
    var Subject Estimate StdErrPred;
    title "First 10 Empirical Bayes Estimates (Random Intercepts)";
run;

/* Sort them from worse to best */
proc sort data=eb_predictions;
    by Estimate;
run;

/* 2. Create a Rank variable for plotting 
data eb_sorted;
    set eb_predictions;
    rank = _N_;
    /* Calculate Approximate 95% CI 
    Lower = Estimate - 1.96*StdErrPred;
    Upper = Estimate + 1.96*StdErrPred;
    
    /* Flag significant ones (Does CI cross 0?) 
    if Lower > 0 or Upper < 0 then Significant="Yes";
    else Significant="No";
run; 


proc freq data=eb_sorted;
    tables Significant;
    title "Count of Significant vs Non-Significant Patients";
run;
*/

/* Plot 
proc sgplot data=eb_sorted;
    title "Caterpillar Plot of Random Intercepts";
    highlow x=rank low=Lower high=Upper / lineattrs=(color=gray);
    refline 0 / axis=y lineattrs=(color=black pattern=dash);
    scatter x=rank y=Estimate / group=Significant markerattrs=(symbol=circlefilled size=3);
	styleattrs datacontrastcolors=(Blue Red);
    yaxis label="Random Intercept (Log-Odds Deviation)";
    xaxis label="Patients (Sorted)";
run;

data simple_plot;
    set eb_estimates;
    Patient_Rank = _N_; 
run;
*/

/* Plot*/
proc sgplot data=simple_plot;
    title "Subject Heterogeneity (Simple View)";
    scatter x=Patient_Rank y=Estimate / markerattrs=(symbol=circlefilled size=2 color=blue);
    refline 0 / axis=y; /* The Average Line */
    yaxis label="Patient Risk (Relative to Average)";
    xaxis label="Patients (Sorted)";
run;


%put We can say, since there is an S shape, that our patients are not identical.;
%put We have a subgroup of patients that get worse than the average and another one that get better.;
%put But most people are near the average.;


/* Histogram of Empirical Bayes Predictions */

proc sgplot data=eb_estimates;
    histogram Estimate / nbins=30;
    density Estimate / type=normal;
    title "Distribution of Random Effects (Should be Normal)";
run;

%put The Random-Effects follow a normal distribution centered at 0 as seen in the histogram.;
 


