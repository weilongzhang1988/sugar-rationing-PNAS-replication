/* =============================================================================
   13_tableS9_england_only.do
   PRODUCES:
     - SI Table S9 (Cancer incidence hazard ratios -- England-only subsample)
       -> ${TAB}/tableS9_cancer_HR_england_only.tex

   INPUT:  ${UKB}/251025_sugar_data_for_cox.dta
   ============================================================================= */

capture log close
log using "${LOG}/13_tableS9_england_only.log", replace

use "${UKB}/251025_sugar_data_for_cox.dta", clear

* Clean row labels for the published Table S9 layout
label define study_clean ///
    1 "\hspace{0.3em}$-$27 months" ///
    2 "\hspace{0.3em}$-$21 months" ///
    3 "\hspace{0.3em}$-$15 months" ///
    4 "\hspace{0.3em}$-$9 months (ref.)" ///
    5 "\hspace{0.3em}In utero" ///
    6 "\hspace{0.3em}In utero $+$ 6 months" ///
    7 "\hspace{0.3em}In utero $+$ 12 months" ///
    8 "\hspace{0.3em}In utero $+$ 18 months" ///
    9 "\hspace{0.3em}In utero $+$ 24 months", replace
label values study study_clean

* England-only subsample (UK Biobank field p1647 instance == England == 1)
stset diag_LungCancer_time, failure(diag_LungCancer_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index ///
      if p1647_i0_2 == 1, cluster(yearmobirth)
eststo lung_eng

stset diag_C22_time, id(eid) failure(diag_C22_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index ///
      if p1647_i0_2 == 1, cluster(yearmobirth)
eststo liver_eng

stset diag_C20_Rectum_time, failure(diag_C20_Rectum_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index ///
      if p1647_i0_2 == 1, cluster(yearmobirth)
eststo rectum_eng

stset diag_C61_Prostate_time, failure(diag_C61_Prostate_event)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index ///
      if male == 1 & p1647_i0_2 == 1, cluster(yearmobirth)
eststo prostate_eng

stset diag_C50_Breast_time, failure(diag_C50_Breast_event)
stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index ///
      if male == 0 & p1647_i0_2 == 1, cluster(yearmobirth)
eststo breast_eng

esttab lung_eng liver_eng rectum_eng prostate_eng breast_eng ///
    using "${TAB}/tableS9_cancer_HR_england_only.tex", replace                ///
    eform keep(*study*)                                                       ///
    cells("b(fmt(2) star) ci(par fmt(2))")                                    ///
    mtitles("Lung" "Liver" "Rectum" "Prostate" "Breast")                      ///
    collabels("HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI" "HR" "95\% CI") ///
    refcat(5.study "\emph{Rationing exposure cohorts}", nolabel)              ///
    star(* 0.10 ** 0.05 *** 0.01) label booktabs nonumbers                    ///
    title("SI Table S9: Cancer incidence hazard ratios -- England-only subsample") ///
    addnotes("Cox proportional hazard models. England-only participants (n = 53,207 full sample;" ///
             "23,375 men for prostate; 29,832 women for breast)." ///
             "Controls: smoking, college, PCA1-10, home longitude/latitude, Townsend deprivation index." ///
             "Reference cohort: born July--December 1954." ///
             "Clustered standard errors by month of birth. Exponentiated coefficients." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

log close
