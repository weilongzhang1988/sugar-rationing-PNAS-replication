/* =============================================================================
   02_fig1_cancer_event_study.do
   PRODUCES:
     - Figure 1 (main paper) component panels (combined externally into
       sugar_fig_cancer.png):
           Lung   -> 251024_esd_Lung_cancer_coefplot.pdf
           Breast -> 251024_esd_Breast_cancer_coefplot.pdf
           Liver  -> 251024_esd_Liver_cancer_coefplot.pdf
           Prostate -> 251024_esd_Prostate_cancer_coefplot.pdf
           Rectum -> 251024_esd_Rectum_cancer_coefplot.pdf
     - SI Table S2 (cancer hazard ratios)  -> tableS2_cancer_HR.tex

   INPUT:  ${UKB}/251025_sugar_data_for_cox.dta
   ============================================================================= */

capture log close
log using "${LOG}/02_fig1_cancer_event_study.log", replace

use "${UKB}/251025_sugar_data_for_cox.dta", clear

* -------------------------------------------------------------------------
* (1) Lung cancer
* -------------------------------------------------------------------------
stset diag_LungCancer_time, failure(diag_LungCancer_event)

stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo lung

di _n "=== Joint pre-trend test: Lung cancer (study 1-3 vs ref study 4) ==="
testparm 1.study 2.study 3.study

coefplot lung, omitted keep(*study*) vertical eform                          ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.00(0.2)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig1A_Lung_cancer_coefplot.pdf", replace

* -------------------------------------------------------------------------
* (2) Breast cancer (women only)
* -------------------------------------------------------------------------
stset diag_C50_Breast_time, failure(diag_C50_Breast_event)

stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male == 0, ///
      cluster(yearmobirth)
eststo breast

di _n "=== Joint pre-trend test: Breast cancer ==="
testparm 1.study 2.study 3.study

coefplot breast, omitted keep(*study*) vertical eform                        ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.50(0.2)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig1B_Breast_cancer_coefplot.pdf", replace

* -------------------------------------------------------------------------
* (3) Liver/IHBD cancer (ICD-10 C22.0-C22.9)
* -------------------------------------------------------------------------
stset diag_C22_time, id(eid) failure(diag_C22_event)

stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo liver

di _n "=== Joint pre-trend test: Liver/IHBD cancer ==="
testparm 1.study 2.study 3.study

coefplot liver, omitted keep(*study*) vertical eform                         ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.00(0.2)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig1C_Liver_cancer_coefplot.pdf", replace

* -------------------------------------------------------------------------
* (4) Prostate cancer (men only)
* -------------------------------------------------------------------------
stset diag_C61_Prostate_time, failure(diag_C61_Prostate_event)

stcox ib4.study i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index if male == 1, ///
      cluster(yearmobirth)
eststo prostate

di _n "=== Joint pre-trend test: Prostate cancer ==="
testparm 1.study 2.study 3.study

coefplot prostate, omitted keep(*study*) vertical eform                      ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.00(0.2)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig1D_Prostate_cancer_coefplot.pdf", replace

* -------------------------------------------------------------------------
* (5) Rectum cancer (ICD-10 C20)
* -------------------------------------------------------------------------
stset diag_C20_Rectum_time, failure(diag_C20_Rectum_event)

stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo rectum

di _n "=== Joint pre-trend test: Rectum cancer ==="
testparm 1.study 2.study 3.study

coefplot rectum, omitted keep(*study*) vertical eform                        ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.00(0.2)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall)) ///
    xtitle("Months of exposure to rationing") ///
    ytitle("Hazard ratio (95% CI)") ///
    graphregion(color(white)) scheme(s1color) ///
    base mlabposition(1)
graph export "${FIG}/Fig1E_Rectum_cancer_coefplot.pdf", replace

* -------------------------------------------------------------------------
* SI Table S2: Cancer HR table (used in appendix.tex tab:cox_bmj_polished)
* -------------------------------------------------------------------------
esttab lung liver rectum prostate breast using "${TAB}/tableS2_cancer_HR.tex", ///
    eform keep(*study*) replace                                               ///
    cells("b(fmt(2) star) ci(par fmt(2))")                                    ///
    collabels("Lung" "Liver/IHBD" "Rectum" "Prostate" "Breast")               ///
    varlabels(*study* "Exposure duration (months)")                           ///
    refcat(5.study "Reference: born post-rationing", nolabel)                 ///
    star(* 0.10 ** 0.05 *** 0.01)                                             ///
    label booktabs nomtitles nonumbers                                        ///
    title("SI Table S2: Hazard ratios (95\% CI) for cancer incidence by exposure duration to sugar rationing") ///
    addnotes("Cox proportional hazard estimates with clustered standard errors by month of birth." ///
             "Exponentiated coefficients reported as hazard ratios."          ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

log close
