/* =============================================================================
   09_KM_logrank_tableS8.do
   PRODUCES:
     - SI Table S8: Log-rank test results for equality of survival functions
       (tab:disease_outcomes_chisquare). The chi-squared statistic and p-value
       for each outcome are written to the log; transcribe into the LaTeX
       table by hand.

   INPUT:  ${UKB}/251109_sugar_data_for_KM_curves.dta
   ============================================================================= */

capture log close
log using "${LOG}/09_KM_logrank_tableS8.log", replace

use "${UKB}/251109_sugar_data_for_KM_curves.dta", clear

foreach disease in diag_hyp_both_w diab_dm_w3min ///
                   diag_LungCancer_time diag_C50_Breast_time ///
                   diag_C22_Liver_time diag_C20_Rectum_time ///
                   diag_C61_Prostate_time {
    capture drop `disease'_time
    gen `disease'_time = `disease'
    replace `disease'_time = age_max if `disease' == .
    capture drop `disease'_event
    gen `disease'_event = (`disease' != .)
    stset `disease'_time, failure(`disease'_event)

    di _n "=== Log-rank test: `disease' ==="
    sts test sugar_rationed2
    sts test sugar_rationed2, wilcoxon
}

log close
