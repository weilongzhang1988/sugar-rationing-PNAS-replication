/* ============================================================================ */
/*  Confounder Sensitivity Analysis — LTL and Granzyme B (Reviewer 2, Q5)      */
/*                                                                              */
/*  Purpose: Augment the baseline OLS event-study models for leukocyte telomere */
/*  length (LTL) and Granzyme B with potential confounders: BMI, C-reactive     */
/*  protein (CRP), lymphocyte percentage, and smoking status.                   */
/*                                                                              */
/*  Note: BMI, CRP, and smoking plausibly lie ON the causal pathway from        */
/*  early-life sugar exposure to biological ageing; conditioning on them        */
/*  therefore risks attenuating genuine effects (post-treatment bias).          */
/*  We present both baseline and augmented specifications so the reader         */
/*  can assess the degree of attenuation.                                       */
/*                                                                              */
/*  Output:                                                                     */
/*    - 260423_confounder_sensitivity_LTL.tex   (Panel A)                       */
/*    - 260423_confounder_sensitivity_GZMB.tex  (Panel B)                       */
/*    - 260423_confounder_sensitivity_LTL_coefplot.pdf                          */
/*    - 260423_confounder_sensitivity_GZMB_coefplot.pdf                         */
/* ============================================================================ */

clear all
set more off
log using "${LOG}/12_tableS17_S18_confounder_LTL_GZMB.log", replace

use "${UKB}/251029_sugar_data_for_cox.dta", clear

/* -------------------------------------------------------------------------- */
/*  VARIABLE DEFINITIONS                                                      */
/*                                                                            */
/*  CRP:            UK Biobank field 30710 (mg/L)                             */
/*  Lymphocyte %:   UK Biobank field 30180 (%)                                */
/* -------------------------------------------------------------------------- */

rename (creactiveproteininstance0 lymphocytepercentageinstance0)(crp lymphocyte_pct)

* --- Log-transform CRP (right-skewed) ---
capture confirm variable crp
if !_rc {
    gen ln_crp = ln(crp) if crp > 0 & !missing(crp)
    label variable ln_crp "Log C-reactive protein (mg/L)"
}

* --- Summarise confounder availability ---
di _n "=== Confounder availability ==="
foreach var in bmi crp ln_crp smoking lymphocyte_pct {
    capture confirm variable `var'
    if !_rc {
        qui count if !missing(`var')
        di "`var':  N = " r(N)
    }
    else {
        di "`var':  *** NOT FOUND — check variable name ***"
    }
}

* --- Sample flag: non-missing on all confounders ---
gen confounder_sample = !missing(bmi) & !missing(ln_crp) & ///
    !missing(lymphocyte_pct) & !missing(smoking)
tab confounder_sample

di _n "=== Observations with complete confounder data ==="
count if confounder_sample == 1


/* ============================================================================ */
/*  PANEL A: LEUKOCYTE TELOMERE LENGTH (LTL)                                   */
/* ============================================================================ */

eststo clear

* --- Model 1: Baseline (replicates main specification) ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index, ///
    cluster(yearmobirth)
eststo ltl_base

* --- Model 2: Baseline + BMI ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    bmi, cluster(yearmobirth)
eststo ltl_bmi

* --- Model 3: Baseline + CRP ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    ln_crp, cluster(yearmobirth)
eststo ltl_crp

* --- Model 4: Baseline + Lymphocyte % ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    lymphocyte_pct, cluster(yearmobirth)
eststo ltl_lymph

* --- Model 5: Baseline + Smoking ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    smoking, cluster(yearmobirth)
eststo ltl_smoke

* --- Model 6: Full model (all confounders) ---
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    bmi ln_crp lymphocyte_pct smoking, ///
    cluster(yearmobirth)
eststo ltl_full

* --- Export LTL table ---
esttab ltl_base ltl_bmi ltl_crp ltl_lymph ltl_smoke ltl_full ///
    using "${TAB}/tableS17_LTL_confounder_sensitivity.tex", ///
    keep(*study*) ///
    cells("b(fmt(3) star) ci(par fmt(3))") ///
    label nonumbers ///
    mtitles("Baseline" "+ BMI" "+ CRP" "+ Lymphocyte" "+ Smoking" "Full") ///
    booktabs ///
    title("Panel A: Confounder sensitivity — Leukocyte telomere length (Z-adjusted)") ///
    addnotes("OLS regressions with clustered standard errors by year-month of birth." ///
             "Baseline controls: sex, month of birth, age, genetic PCs 1--10," ///
             "home coordinates, Townsend deprivation index." ///
             "BMI, CRP, and smoking are potential mediators (post-treatment);" ///
             "conditioning on them may attenuate genuine causal effects." ///
             "* p<0.10, ** p<0.05, *** p<0.01.") ///
    replace

* (LTL baseline-vs-full coefplot is not used in the published appendix.)


/* ============================================================================ */
/*  PANEL B: GRANZYME B                                                        */
/* ============================================================================ */

eststo clear

* --- Model 1: Baseline (replicates main specification, adding cluster SEs) ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index, ///
    cluster(yearmobirth)
eststo gzmb_base

* --- Model 2: Baseline + BMI ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    bmi, cluster(yearmobirth)
eststo gzmb_bmi

* --- Model 3: Baseline + CRP ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    ln_crp, cluster(yearmobirth)
eststo gzmb_crp

* --- Model 4: Baseline + Lymphocyte % ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    lymphocyte_pct, cluster(yearmobirth)
eststo gzmb_lymph

* --- Model 5: Baseline + Smoking ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    smoking, cluster(yearmobirth)
eststo gzmb_smoke

* --- Model 6: Full model (all confounders) ---
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 ///
    age home_longi home_lati Townsend_deprivation_index ///
    bmi ln_crp lymphocyte_pct smoking, ///
    cluster(yearmobirth)
eststo gzmb_full

* --- Export Granzyme B table ---
esttab gzmb_base gzmb_bmi gzmb_crp gzmb_lymph gzmb_smoke gzmb_full ///
    using "${TAB}/tableS18_GZMB_confounder_sensitivity.tex", ///
    keep(*study*) ///
    cells("b(fmt(3) star) ci(par fmt(3))") ///
    label nonumbers ///
    mtitles("Baseline" "+ BMI" "+ CRP" "+ Lymphocyte" "+ Smoking" "Full") ///
    booktabs ///
    title("Panel B: Confounder sensitivity — Granzyme B (log$_2$ NPX)") ///
    addnotes("OLS regressions with clustered standard errors by year-month of birth." ///
             "Baseline controls: sex, month of birth, age, genetic PCs 1--10," ///
             "home coordinates, Townsend deprivation index." ///
             "BMI, CRP, and smoking are potential mediators (post-treatment);" ///
             "conditioning on them may attenuate genuine causal effects." ///
             "* p<0.10, ** p<0.05, *** p<0.01.") ///
    replace

* (Granzyme B baseline-vs-full coefplot is not used in the published appendix.)
* (Selection-into-confounder-sample logit table is not used in the published appendix.)


log close
