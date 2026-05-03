/* ============================================================================ */
/*  Multiple Hypothesis Testing Correction -- Bonferroni                       */
/*                                                                             */
/*  This script does not produce a stand-alone appendix table. It computes the */
/*  Bonferroni-adjusted p-values reported INLINE within SI Table S2            */
/*  (tab:cox_bmj_polished) -- specifically the "Bonf. p" rows beneath each     */
/*  cancer column and the joint chi-squared / Bonferroni p row. The per-       */
/*  outcome and joint p-values must be transcribed from the matrix output      */
/*  below into the LaTeX source of Table S2.                                   */
/*                                                                             */
/*  Five primary cancer outcomes: lung, liver, rectum, prostate, breast        */
/*  Two testing strategies:                                                    */
/*    A. Most-exposed cohort (study==9, in-utero+24m) coefficient              */
/*    B. Joint test of all rationed cohorts (study 5-9)                        */
/* ============================================================================ */

clear all
set more off
log using "${LOG}/15_multiple_testing_correction.log", replace

/* ============================================================================ */
/*  PART 1: BONFERRONI CORRECTION                                              */
/* ============================================================================ */

use "${UKB}/251025_sugar_data_for_cox.dta", clear

local n_outcomes = 5

* --- Storage matrices ---
matrix pvals_coef  = J(5, 4, .)   // outcome × {unadj_p, bonf_p, z, HR}
matrix pvals_joint = J(5, 3, .)   // outcome × {unadj_p, bonf_p, chi2}
matrix rownames pvals_coef  = Lung Liver Rectum Prostate Breast
matrix colnames pvals_coef  = "Unadj_p" "Bonf_p" "z_stat" "HR"
matrix rownames pvals_joint = Lung Liver Rectum Prostate Breast
matrix colnames pvals_joint = "Unadj_p" "Bonf_p" "chi2"

* ------------------------------------------------------------------
* (1) Lung cancer
* ------------------------------------------------------------------
stset diag_LungCancer_time, failure(diag_LungCancer_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)

* Most-exposed cohort (study==9) coefficient
local hr = exp(_b[9.study])
local z  = _b[9.study] / _se[9.study]
local p  = 2 * normal(-abs(`z'))
matrix pvals_coef[1,1] = `p'
matrix pvals_coef[1,2] = min(1, `p' * `n_outcomes')
matrix pvals_coef[1,3] = `z'
matrix pvals_coef[1,4] = `hr'

* Joint test: all rationed cohorts (study 5-9)
testparm 5.study 6.study 7.study 8.study 9.study
matrix pvals_joint[1,1] = r(p)
matrix pvals_joint[1,2] = min(1, r(p) * `n_outcomes')
matrix pvals_joint[1,3] = r(chi2)

* ------------------------------------------------------------------
* (2) Liver cancer (C22)
* ------------------------------------------------------------------
stset diag_C22_time, id(eid) failure(diag_C22_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)

local hr = exp(_b[9.study])
local z  = _b[9.study] / _se[9.study]
local p  = 2 * normal(-abs(`z'))
matrix pvals_coef[2,1] = `p'
matrix pvals_coef[2,2] = min(1, `p' * `n_outcomes')
matrix pvals_coef[2,3] = `z'
matrix pvals_coef[2,4] = `hr'

testparm 5.study 6.study 7.study 8.study 9.study
matrix pvals_joint[2,1] = r(p)
matrix pvals_joint[2,2] = min(1, r(p) * `n_outcomes')
matrix pvals_joint[2,3] = r(chi2)

* ------------------------------------------------------------------
* (3) Rectum cancer (C20)
* ------------------------------------------------------------------
stset diag_C20_Rectum_time, failure(diag_C20_Rectum_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca5 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)

local hr = exp(_b[9.study])
local z  = _b[9.study] / _se[9.study]
local p  = 2 * normal(-abs(`z'))
matrix pvals_coef[3,1] = `p'
matrix pvals_coef[3,2] = min(1, `p' * `n_outcomes')
matrix pvals_coef[3,3] = `z'
matrix pvals_coef[3,4] = `hr'

testparm 5.study 6.study 7.study 8.study 9.study
matrix pvals_joint[3,1] = r(p)
matrix pvals_joint[3,2] = min(1, r(p) * `n_outcomes')
matrix pvals_joint[3,3] = r(chi2)

* ------------------------------------------------------------------
* (4) Prostate cancer (C61, men only)
* ------------------------------------------------------------------
stset diag_C61_Prostate_time, failure(diag_C61_Prostate_event)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==1, ///
      cluster(yearmobirth)

local hr = exp(_b[9.study])
local z  = _b[9.study] / _se[9.study]
local p  = 2 * normal(-abs(`z'))
matrix pvals_coef[4,1] = `p'
matrix pvals_coef[4,2] = min(1, `p' * `n_outcomes')
matrix pvals_coef[4,3] = `z'
matrix pvals_coef[4,4] = `hr'

testparm 5.study 6.study 7.study 8.study 9.study
matrix pvals_joint[4,1] = r(p)
matrix pvals_joint[4,2] = min(1, r(p) * `n_outcomes')
matrix pvals_joint[4,3] = r(chi2)

* ------------------------------------------------------------------
* (5) Breast cancer (C50, women only)
* ------------------------------------------------------------------
stset diag_C50_Breast_time, failure(diag_C50_Breast_event)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==0, ///
      cluster(yearmobirth)

local hr = exp(_b[9.study])
local z  = _b[9.study] / _se[9.study]
local p  = 2 * normal(-abs(`z'))
matrix pvals_coef[5,1] = `p'
matrix pvals_coef[5,2] = min(1, `p' * `n_outcomes')
matrix pvals_coef[5,3] = `z'
matrix pvals_coef[5,4] = `hr'

testparm 5.study 6.study 7.study 8.study 9.study
matrix pvals_joint[5,1] = r(p)
matrix pvals_joint[5,2] = min(1, r(p) * `n_outcomes')
matrix pvals_joint[5,3] = r(chi2)

* --- Display Bonferroni results ---
di _n "=== BONFERRONI CORRECTION: Most-exposed cohort (study==9) ==="
matrix list pvals_coef, format(%9.4f)

di _n "=== BONFERRONI CORRECTION: Joint test (study 5-9) ==="
matrix list pvals_joint, format(%9.4f)


di _n "=== Done. Results saved to: ==="
di    "  260424_multiple_testing_results.dta"
di    "  260424_multiple_testing_table.tex"

log close
