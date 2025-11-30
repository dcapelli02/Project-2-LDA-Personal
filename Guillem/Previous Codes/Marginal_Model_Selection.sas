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

/* Unstrucutred */
proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX /
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
run;

/* ALR */

proc genmod data=alzheimer_long descending;
    class PATID SEX TIMECLSS;

    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX/
        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss logor=FULLCLUST covb corrw modelse;
run;


/* Linearization Based Method */


/* Unstrucutred */

proc glimmix data=alzheimer_long method=RSPL empirical;
	class PATID SEX TIMECLSS;
	model CDRSB_CAT(descending) = TIME SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0 
     TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * ABPET0 TIME * TAUPET0
     TIME * SEX
		/ dist=binary solution;
	random _residual_ / subject=PATID type=UN;
run;


%put Between GEE, ALR and LBM, the coefficients are almost the same;
%put So we can consider the most simple one for us to reduce it;


/* Bahadur Model */

data alz_bahadur_clean;
    if _n_ = 1 then set stats; 
    set alzheimer25;

    AGE_STD = (AGE - mean_age) / sd_age;
    BMI_STD = (BMI - mean_bmi) / sd_bmi;

    array raw[0:6] cdrsb0-cdrsb6;
    array Y[0:6] Y0-Y6; 

    do i = 0 to 6;
        if raw[i] ne . then Y[i] = (raw[i] >= 10); 
        else Y[i] = .;
    end;
    
    if nmiss(SEX, AGE_STD, BMI_STD, ADL, ABPET0, TAUPET0) > 0 then delete;

    keep PATID Y0-Y6 SEX AGE_STD BMI_STD ADL ABPET0 TAUPET0;
run;

proc nlmixed data=alz_bahadur_clean qpoints=1 MAXITER=100;

    parms b0=-1 b_t=0.2 b_sex=0 b_age=0 b_bmi=0 b_adl=0 b_ab=0 b_tau=0 
          b_t_age=0 b_t_bmi=0 b_t_adl=0 b_t_ab=0 b_t_tau=0 b_t_sex=0 
          rho=0.1;

    array Y_arr[7] Y0-Y6;
    ll_indep = 0; 
    corr_sum = 0;

    do t = 0 to 6;
        idx = t + 1; 
        
        if Y_arr[idx] ne . then do;
            /* A. MARGINAL MEAN */
            eta = b0 + b_t*t + b_sex*SEX + b_age*AGE_STD + b_bmi*BMI_STD + 
                  b_adl*ADL + b_ab*ABPET0 + b_tau*TAUPET0 +
                  b_t_age*(t*AGE_STD) + b_t_bmi*(t*BMI_STD) + b_t_adl*(t*ADL) +
                  b_t_ab*(t*ABPET0)   + b_t_tau*(t*TAUPET0) + b_t_sex*(t*SEX);

            p = 1 / (1 + exp(-eta));

            /* SAFETY RAIL 1: Prevent P from being exactly 0 or 1 */
            if p < 1e-6 then p = 1e-6;
            if p > (1 - 1e-6) then p = (1 - 1e-6);

            /* Now safe to divide */
            z = (Y_arr[idx] - p) / sqrt(p*(1-p));
            
            ll_indep = ll_indep + (Y_arr[idx]*log(p) + (1-Y_arr[idx])*log(1-p));

            /* B. CORRELATION LOOP */
            do k = 0 to t-1;
                k_idx = k + 1;
                if Y_arr[k_idx] ne . then do;
                    eta_k = b0 + b_t*k + b_sex*SEX + b_age*AGE_STD + b_bmi*BMI_STD + 
                            b_adl*ADL + b_ab*ABPET0 + b_tau*TAUPET0 +
                            b_t_age*(k*AGE_STD) + b_t_bmi*(k*BMI_STD) + b_t_adl*(k*ADL) +
                            b_t_ab*(k*ABPET0)   + b_t_tau*(k*TAUPET0) + b_t_sex*(k*SEX);
                    
                    p_k = 1 / (1 + exp(-eta_k));
                    
                    /* SAFETY RAIL 2: Apply same bounds to previous point */
                    if p_k < 1e-6 then p_k = 1e-6;
                    if p_k > (1 - 1e-6) then p_k = (1 - 1e-6);

                    z_k = (Y_arr[k_idx] - p_k) / sqrt(p_k*(1-p_k));
                    
                    corr_sum = corr_sum + (rho * z * z_k);
                end;
            end;
        end;
    end;

    /* C. FINAL LIKELIHOOD */
    cf = 1 + corr_sum;
    
    /* SAFETY RAIL 3: Ensure correction factor is positive */
    if cf <= 1e-6 then cf = 1e-6; 
    
    ll = ll_indep + log(cf);
    
    model Y0 ~ general(ll);
run;


%put Correlation component of the Bahadur likelihood is overpowering the Mean component. 
%put The model is so desperate to make the correlations valid that it is breaking the regression line.
%put That`s why we should not use the Bahadur model for our data.





