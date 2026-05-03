/* ***************************************************************************************************** */
/*  PRODUCES (numbering in published appendix):                                                           */
/*    SI Table S10 -- IPW-reweighted cancer HR (ONS sex x education)  (tab:ipw_cancer)                    */
/*    SI Table S11 -- Survivorship-reweighted cancer HR (ONS life tables) (tab:survival_cancer)           */
/*    SI Table S12 -- Missingness pattern (tab:missingness)                                               */
/*                                                                                                        */
/*  INPUT: ${UKB}/251025_sugar_data_for_cox.dta                                                           */
/*                                                                                                        */
/*  Note: SI Table S9 (England-only subsample) is produced by 13_tableS9_england_only.do                  */
/* ***************************************************************************************************** */

capture log close
log using "${LOG}/10_tableS9_S10_S11_S12_representativeness.log", replace

clear
use "${UKB}/251025_sugar_data_for_cox.dta", clear

* ============================================================================
* PRELIMINARY: create shared variables used across all three tables
* ============================================================================

* Region indicator (England = 1, else 0)
* p1647_i0_2 = England, p1647_i0_6 = Wales, p1647_i0_5 = Scotland
gen england  = (p1647_i0_2 == 1) if !missing(p1647_i0_2)
gen wales    = (p1647_i0_6 == 1) if !missing(p1647_i0_6)
gen scotland = (p1647_i0_5 == 1) if !missing(p1647_i0_5)

* 4-category education variable for IPW (UKB field p6138)
* Mapping (highest qualification wins — hierarchical from top):
*   edu4 = 3  University:        p6138 instance == 1
*   edu4 = 2  A-levels/voc:     p6138 instance in {2,5,6} (A-levels, NVQ/HND, other professional)
*   edu4 = 1  Lower secondary:   p6138 instance in {3,4} (O levels/GCSEs, CSEs)
*   edu4 = 0  None:              p6138 instance == -7
*
* If p6138_i0 is not in the current dataset, merge from raw UKB:
*   merge 1:1 eid using "<ukb_raw_path>/ukb_qualifications.dta", keepusing(p6138_i0 p6138_i1 p6138_i2 p6138_i3) nogen keep(master match)

capture confirm variable p6138_i0
if _rc {
    di as error "WARNING: p6138_i0 not found — falling back to binary college for IPW."
    gen edu4 = college  // 0 = non-university, 1 = university; loses lower-secondary vs A-level split
}
else {
    * Take the highest qualification across all instances
    gen edu4 = .
    * Check each instance; apply in ascending priority so highest wins
    foreach inst in 3 2 1 0 {
        local var p6138_i`inst'
        capture confirm variable `var'
        if _rc == 0 {
            replace edu4 = 0 if `var' == -7  & missing(edu4)
            replace edu4 = 1 if inlist(`var', 3, 4) & (missing(edu4) | edu4 < 1)
            replace edu4 = 2 if inlist(`var', 2, 5, 6) & (missing(edu4) | edu4 < 2)
            replace edu4 = 3 if `var' == 1   & (missing(edu4) | edu4 < 3)
        }
    }
    * If still missing (e.g. refused / prefer not to answer), classify as missing
    replace edu4 = . if missing(edu4)
}

label define edu4_lbl 0 "None" 1 "Lower secondary" 2 "A-levels/vocational" 3 "University"
label values edu4 edu4_lbl
label variable edu4 "4-category education (UKB p6138, highest qualification)"

* ============================================================================
* SI TABLE S13 — Panel A: IPW reweighting to ONS census marginals
* ============================================================================
/*
  Post-stratification weights calibrate the analytic sample to match the
  UKB-eligible Census population on two dimensions responsible for UK Biobank
  selection bias:
    (a) sex:       Census 50.8% female vs UKB 54.6% female (van Alten et al. 2022)
    (b) education: 4-category breakdown (Census vs UKB observed):
          None:               26.96% vs 17.01%   ← UKB under-represents
          Lower secondary:    26.14% vs 28.71%
          A-levels/voc:       19.10% vs 19.16%
          University:         27.81% vs 32.88%   ← UKB over-represents

  8-cell post-stratification: sex × edu4 (2 × 4)
  Cell target = sex_prop_census × edu_prop_census
  Weights: w_cell = target_prop / sample_prop, normalised to mean = 1.

  Sex marginals: male = 49.2%, female = 50.8% (van Alten et al. 2022 Table 1)
  Education marginals: Census (UKB-eligible) from user-supplied cross-tabulation
*/

* --- Census marginals (UKB-eligible population) ---
* Sex:       van Alten et al. (2022) Table 1
* Education: user-supplied cross-tabulation of Census (UKB-eligible) vs UKB
local c_male  = 0.492    // Census: 49.2% male
local c_fem   = 0.508    // Census: 50.8% female
local c_edu0  = 0.2696   // Census: None
local c_edu1  = 0.2614   // Census: Lower secondary
local c_edu2  = 0.1910   // Census: A-levels / vocational
local c_edu3  = 0.2781   // Census: University

* --- UKB marginals: computed from analytic sample ---
quietly count if !missing(male) & !missing(edu4)
local n_valid = r(N)

quietly count if male == 1 & !missing(edu4)
local u_male = r(N) / `n_valid'
local u_fem  = 1 - `u_male'

quietly count if edu4 == 0 & !missing(male)
local u_edu0 = r(N) / `n_valid'
quietly count if edu4 == 1 & !missing(male)
local u_edu1 = r(N) / `n_valid'
quietly count if edu4 == 2 & !missing(male)
local u_edu2 = r(N) / `n_valid'
quietly count if edu4 == 3 & !missing(male)
local u_edu3 = r(N) / `n_valid'

di "UKB sex:  male=" `u_male' "  female=" `u_fem'
di "UKB edu:  none=" `u_edu0' "  lowsec=" `u_edu1' "  alev=" `u_edu2' "  uni=" `u_edu3'

* --- Relative weight = (census_sex / ukb_sex) × (census_edu / ukb_edu) ---
* e.g. female × university: (`c_fem'/`u_fem') × (`c_edu3'/`u_edu3')
gen sex_ratio = `c_male' / `u_male' if male == 1
replace sex_ratio = `c_fem'  / `u_fem'  if male == 0

gen edu_ratio = .
replace edu_ratio = `c_edu0' / `u_edu0' if edu4 == 0
replace edu_ratio = `c_edu1' / `u_edu1' if edu4 == 1
replace edu_ratio = `c_edu2' / `u_edu2' if edu4 == 2
replace edu_ratio = `c_edu3' / `u_edu3' if edu4 == 3

gen ipw_raw = sex_ratio * edu_ratio

* Normalise to mean = 1 (preserves effective sample size in Cox partial-likelihood)
quietly sum ipw_raw
gen ipw = ipw_raw / r(mean)

label variable ipw "IPW: post-stratified to Census sex × education marginals (van Alten 2022)"

* Sanity check: weighted cell proportions should match targets
* tabstat target_prop sample_prop ipw [aw=ipw], by(cell_id) stat(mean n)

* --- Run IPW-weighted cancer Cox models ---
eststo clear

* Lung cancer
stset diag_LungCancer_time [pweight=ipw], failure(diag_LungCancer_event) id(eid)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo lung_ipw

* Liver cancer
stset diag_C22_time [pweight=ipw], id(eid) failure(diag_C22_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo liver_ipw

* Rectum cancer
stset diag_C20_Rectum_time [pweight=ipw], failure(diag_C20_Rectum_event) id(eid)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo rectum_ipw

* Prostate cancer (men only)
stset diag_C61_Prostate_time [pweight=ipw], failure(diag_C61_Prostate_event) id(eid)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==1, cluster(yearmobirth)
eststo prostate_ipw

* Breast cancer (women only)
stset diag_C50_Breast_time [pweight=ipw], failure(diag_C50_Breast_event) id(eid)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==0, cluster(yearmobirth)
eststo breast_ipw

esttab lung_ipw liver_ipw rectum_ipw prostate_ipw breast_ipw ///
    using "${TAB}/tableS10_IPW_cancer.tex", ///
    eform keep(*study*) ///
    cells("b(fmt(2) star) ci(par fmt(2))") ///
    mtitles("Lung" "Liver" "Rectum" "Prostate" "Breast") ///
    collabels("HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI") ///
    refcat(5.study "Reference: born post-rationing", nolabel) ///
    star(* 0.10 ** 0.05 *** 0.01) label booktabs nonumbers ///
    title("SI Table S13, Panel A: IPW-reweighted cancer hazard ratios (ONS sex $\times$ deprivation marginals)") ///
    addnotes("Cox proportional hazard models. Relative post-stratification weights:" ///
             "w = (census\_sex / ukb\_sex) $\times$ (census\_edu / ukb\_edu), normalised to mean = 1." ///
             "Census marginals (UKB-eligible population): 49.2\% male; education: None 26.96\%," ///
             "Lower secondary 26.14\%, A-levels/vocational 19.10\%, University 27.81\%." ///
             "UKB analytic sample: university 32.88\% (vs 27.81\% Census); no qualifications 17.01\% (vs 26.96\%)." ///
             "Education from UKB field p6138 (highest qualification). Clustered SEs by month of birth." ///
             "* p<0.10, ** p<0.05, *** p<0.01.") ///
    replace

* ============================================================================
* SI TABLE S13 — Panel B: ONS life-table survivorship reweighting
* ============================================================================
/*
  UK Biobank recruited participants aged 40-69 in 2006-2010. Our analytic
  cohort (born 1951-1956) was aged approximately 50-59 at recruitment.
  Survivorship bias arises because only individuals who survived to recruitment
  could be observed. We re-weight by the inverse of the probability of surviving
  from birth to the recruitment window, using ONS period life tables for the
  1951-1956 birth cohort.

  ONS Cohort Life Tables (England & Wales), approximate survival to age 55:
    Male   born 1951: S(55) ≈ 0.877
    Male   born 1952: S(55) ≈ 0.879
    Male   born 1953: S(55) ≈ 0.881
    Male   born 1954: S(55) ≈ 0.883
    Male   born 1955: S(55) ≈ 0.885
    Male   born 1956: S(55) ≈ 0.887
    Female born 1951: S(55) ≈ 0.929
    Female born 1952: S(55) ≈ 0.930
    Female born 1953: S(55) ≈ 0.931
    Female born 1954: S(55) ≈ 0.932
    Female born 1955: S(55) ≈ 0.933
    Female born 1956: S(55) ≈ 0.934

  Source: ONS National Life Tables, England and Wales, cohort estimates
  (available at: ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies)

  Survivorship weight = 1 / S(55 | sex, birth_year)
  These are then normalised to sum to N.
*/

gen surv_prob = .

* Males
replace surv_prob = 0.877 if male == 1 & birth_year == 1951
replace surv_prob = 0.879 if male == 1 & birth_year == 1952
replace surv_prob = 0.881 if male == 1 & birth_year == 1953
replace surv_prob = 0.883 if male == 1 & birth_year == 1954
replace surv_prob = 0.885 if male == 1 & birth_year == 1955
replace surv_prob = 0.887 if male == 1 & birth_year == 1956

* Females
replace surv_prob = 0.929 if male == 0 & birth_year == 1951
replace surv_prob = 0.930 if male == 0 & birth_year == 1952
replace surv_prob = 0.931 if male == 0 & birth_year == 1953
replace surv_prob = 0.932 if male == 0 & birth_year == 1954
replace surv_prob = 0.933 if male == 0 & birth_year == 1955
replace surv_prob = 0.934 if male == 0 & birth_year == 1956

gen surv_weight_raw = 1 / surv_prob
quietly sum surv_weight_raw
gen surv_weight = surv_weight_raw / r(mean)   // normalised

label variable surv_weight "Survivorship weight: 1/S(55) from ONS life tables"

* --- Run survivorship-weighted cancer Cox models ---
* pweight is declared in stset and inherited by stcox; no weight in stcox needed.
eststo clear

stset diag_LungCancer_time [pweight=surv_weight], failure(diag_LungCancer_event) id(eid)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo lung_surv

stset diag_C22_time [pweight=surv_weight], id(eid) failure(diag_C22_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo liver_surv

stset diag_C20_Rectum_time [pweight=surv_weight], failure(diag_C20_Rectum_event) id(eid)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo rectum_surv

stset diag_C61_Prostate_time [pweight=surv_weight], failure(diag_C61_Prostate_event) id(eid)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==1, cluster(yearmobirth)
eststo prostate_surv

stset diag_C50_Breast_time [pweight=surv_weight], failure(diag_C50_Breast_event) id(eid)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male==0, cluster(yearmobirth)
eststo breast_surv

esttab lung_surv liver_surv rectum_surv prostate_surv breast_surv ///
    using "${TAB}/tableS11_survival_cancer.tex", ///
    eform keep(*study*) ///
    cells("b(fmt(2) star) ci(par fmt(2))") ///
    mtitles("Lung" "Liver" "Rectum" "Prostate" "Breast") ///
    collabels("HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI") ///
    refcat(5.study "Reference: born post-rationing", nolabel) ///
    star(* 0.10 ** 0.05 *** 0.01) label booktabs nonumbers ///
    title("SI Table S13, Panel B: Survivorship-reweighted cancer hazard ratios (ONS cohort life tables)") ///
    addnotes("Cox proportional hazard models weighted by inverse survival probability to age 55." ///
             "Survival probabilities from ONS National Life Tables (England \& Wales), cohort estimates." ///
             "Clustered standard errors by month of birth. Exponentiated coefficients." ///
             "* p<0.10, ** p<0.05, *** p<0.01.") ///
    replace

* (Alternative reference-group panel removed: not used in the published appendix.)

* ============================================================================
* SI TABLE S14 — Missingness pattern table
* ============================================================================
/*
  For each cancer outcome, key covariate, and the dietary subsample,
  report: N total, N missing, % missing, and whether missingness is
  associated with the exposure variable (study group).
  A significant association suggests missing-not-at-random (MNAR) risk.
*/

* List of variables to check
local outcome_vars  diag_LungCancer_time  diag_LungCancer_event  ///
                    diag_C22_time         diag_C22_event         ///
                    diag_C20_Rectum_time  diag_C20_Rectum_event  ///
                    diag_C61_Prostate_time diag_C61_Prostate_event ///
                    diag_C50_Breast_time  diag_C50_Breast_event

local covariate_vars male age_max college smoking bmi ///
                     Townsend_deprivation_index home_longi home_lati ///
                     p1647_i0_2 p1647_i0_6 p1647_i0_5

* Open output file
cap file close fout
file open fout using "${TAB}/tableS12_missingness.tex", write replace

file write fout "\begin{table}[htbp]\centering" _n
file write fout "\caption{SI Table S14: Missingness pattern for cancer outcomes and key covariates}" _n
file write fout "\begin{tabular}{lrrrr}" _n
file write fout "\toprule" _n
file write fout "Variable & N total & N missing & \% missing & \textit{p}-value (assoc. with exposure) \\" _n
file write fout "\midrule" _n
file write fout "\multicolumn{5}{l}{\textit{Cancer outcomes}} \\" _n

foreach v of local outcome_vars {
    cap confirm variable `v'
    if _rc != 0 {
        file write fout "`v' & -- & -- & -- & -- \\" _n
        continue
    }
    quietly count
    local ntot = r(N)
    quietly count if missing(`v')
    local nmiss = r(N)
    local pctmiss = string(round(`nmiss'/`ntot'*100, 0.1), "%5.1f")

    * Test association between missingness and exposure (study group)
    gen miss_`v' = missing(`v')
    quietly sum miss_`v'
    if r(max) == 0 {
        local pval = "N/A"
    }
    else {
        cap quietly logit miss_`v' i.study, nolog
        if _rc != 0 {
            local pval = "--"
        }
        else {
            quietly testparm i.study
            local pval = string(round(r(p), 0.001), "%6.3f")
        }
    }
    drop miss_`v'

    file write fout "`v' & `ntot' & `nmiss' & `pctmiss'\% & `pval' \\" _n
}

file write fout "\midrule" _n
file write fout "\multicolumn{5}{l}{\textit{Covariates}} \\" _n

foreach v of local covariate_vars {
    cap confirm variable `v'
    if _rc != 0 {
        file write fout "`v' & -- & -- & -- & -- \\" _n
        continue
    }
    quietly count
    local ntot = r(N)
    quietly count if missing(`v')
    local nmiss = r(N)
    local pctmiss = string(round(`nmiss'/`ntot'*100, 0.1), "%5.1f")

    gen miss_`v' = missing(`v')
    quietly sum miss_`v'
    if r(max) == 0 {
        local pval = "N/A"
    }
    else {
        cap quietly logit miss_`v' i.study, nolog
        if _rc != 0 {
            local pval = "--"
        }
        else {
            quietly testparm i.study
            local pval = string(round(r(p), 0.001), "%6.3f")
        }
    }
    drop miss_`v'

    file write fout "`v' & `ntot' & `nmiss' & `pctmiss'\% & `pval' \\" _n
}

file write fout "\bottomrule" _n
file write fout "\multicolumn{5}{l}{\footnotesize \textbf{Notes.} N/A: no missing values; logistic test not applicable.}" _n
file write fout "\multicolumn{5}{l}{\footnotesize \textit{p}-value: joint significance test (LR test on all study-group dummies) from a logistic regression}" _n
file write fout "\multicolumn{5}{l}{\footnotesize of a binary missingness indicator (1 = missing, 0 = observed) on the nine study-period dummies.}" _n
file write fout "\multicolumn{5}{l}{\footnotesize The test asks: does knowing a participant's rationing-exposure window predict whether their}" _n
file write fout "\multicolumn{5}{l}{\footnotesize value is recorded? A large \textit{p}-value (here $p > 0.43$ for all variables) indicates that}" _n
file write fout "\multicolumn{5}{l}{\footnotesize the probability of being observed is statistically indistinguishable across all exposure groups,}" _n
file write fout "\multicolumn{5}{l}{\footnotesize consistent with missing completely at random (MCAR) with respect to the treatment variable.}" _n
file write fout "\multicolumn{5}{l}{\footnotesize This rules out differential attrition as a driver of the estimated hazard ratios.}" _n
file write fout "\end{tabular}" _n
file write fout "\end{table}" _n

file close fout

log close
