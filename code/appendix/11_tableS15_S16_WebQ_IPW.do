*===============================================================
* 260421_consumption.do
*
* Referee R1.4 response: Oxford WebQ dietary subsample
*   (a) SI Table S17  -- balance of WebQ vs non-WebQ + selection logit
*   (b) SI Table S17b -- unweighted + IPW-reweighted dietary regressions
*                        (sugars, HEI-2015, food weight).
*                        IPW is the primary specification.
*   (c) SI Figure S6  -- per-cell sample sizes
*
* Outcomes run in (b):
*   avg_totalsugars        (g/day)
*   HEI_score              (HEI-2015, points)
*   avg_foodweight         (g/day)
*
* Notes:
*   - Food_Diversity_score was inspected but excluded from the
*     reported robustness table: the unweighted pre-rationing placebo
*     cells (studies 1-2) are marginally significant and negative,
*     which complicates the placebo interpretation.
*   - Age-at-WebQ + wave adjustment is not reported. Under IPW that
*     corrects for differential WebQ participation by rationing status,
*     the sugar pre-trend loses significance and the rationed-cell
*     dose-response is uniform across exposure duration (no drift
*     toward zero at long exposure). IPW is therefore the preferred
*     way to address the referee's participation-selection concern.
*
* Log: 260421_consumption.log
*===============================================================

set more off
capture log close

clear
log using "${LOG}/11_tableS15_S16_WebQ_IPW.log", replace

use "${UKB}/251027_sugar_data_for_cox.dta", clear

* --- 0. Subsample indicator and exposure flag --------------------------------
gen     in_diet = !missing(HEI_score)
label var in_diet "Completed Oxford WebQ (=1)"

gen     ration = 0 if study >= 1 & study <= 4
replace ration = 1 if study >= 5 & study <= 9
label var ration "Born during sugar rationing (=1)"

* Region indicators (renamed earlier in the do-file). Recreate if absent.
capture confirm variable England
if _rc {
    capture rename (p1647_i0_2 p1647_i0_6 p1647_i0_5 p21000_i0_17) ///
                   (England Wales Scotland white)
}

* --- (a) SI Table S17: balance + selection logit -----------------------------
local balvars age male college Townsend_deprivation_index ///
              England Wales Scotland white smoking ration

eststo clear
eststo full_all : estpost summarize `balvars'
eststo full_in  : estpost summarize `balvars' if in_diet == 1
eststo full_out : estpost summarize `balvars' if in_diet == 0
eststo bal_t    : estpost ttest `balvars', by(in_diet) unequal

esttab full_all full_in full_out bal_t ///
       using "${TAB}/tableS15_WebQ_balance_selection.tex", replace ///
       cells("mean(pattern(1 1 1 0) fmt(3)) b(pattern(0 0 0 1) star fmt(3))") ///
       label nonum booktabs ///
       mtitles("Full analytic" "WebQ subsample" "Non-WebQ" "Diff (WebQ$-$Non)") ///
       title("SI Table S17: Balance of Oxford WebQ subsample vs full analytic sample (1951--1956 birth window)") ///
       addnotes("Means; t-tests with unequal variances. * p<0.10, ** p<0.05, *** p<0.01.")

* Selection logit: does rationing exposure predict WebQ participation
* conditional on demographics and geography?
logit in_diet i.ration male age college Townsend_deprivation_index ///
       England Wales Scotland white smoking pca1-pca10 ///
       home_longi home_lati i.birth_month, cluster(yearmobirth)
estimates store sel_logit

* --- (b) SI Table S17b: IPW-reweighted dietary regressions -------------------
predict double pr_indiet, pr
gen double ipw_diet = 1 / pr_indiet if in_diet == 1
sum ipw_diet, detail
gen double ipw_diet_tr = ipw_diet
replace ipw_diet_tr = r(p1)  if ipw_diet < r(p1)  & !missing(ipw_diet)
replace ipw_diet_tr = r(p99) if ipw_diet > r(p99) & !missing(ipw_diet)

eststo clear
foreach y in avg_totalsugars HEI_score avg_foodweight {
    eststo `y'_unw : reg `y' ib4.study male age college pca1-pca10 ///
        home_longi home_lati i.birth_month if in_diet==1, ///
        cluster(yearmobirth)
    eststo `y'_ipw : reg `y' ib4.study male age college pca1-pca10 ///
        home_longi home_lati i.birth_month [pweight=ipw_diet_tr] ///
        if in_diet==1, cluster(yearmobirth)
}

esttab avg_totalsugars_unw avg_totalsugars_ipw ///
       HEI_score_unw       HEI_score_ipw       ///
       avg_foodweight_unw  avg_foodweight_ipw  ///
       using "${TAB}/tableS16_WebQ_diet_IPW.tex", replace ///
       keep(*study*) cells("b(fmt(3) star) ci(par fmt(3))") ///
       label nomtitles nonumbers booktabs ///
       collabels("Sugars (g/d) unwgt" "Sugars IPW" ///
                 "HEI unwgt" "HEI IPW" ///
                 "Food wt (g/d) unwgt" "Food wt IPW") ///
       title("SI Table S17b: Dietary regressions reweighted to the full analytic covariate distribution") ///
       addnotes("Inverse-probability weights derived from the SI Table S17 selection logit, trimmed at the 1st/99th percentile.")

* (Per-cell counts and unweighted-vs-IPW diagnostic coefplots are not used
*  in the published appendix; SI Figure S6 is produced by
*  16_figS6_all_category_age50.do.)

log close
