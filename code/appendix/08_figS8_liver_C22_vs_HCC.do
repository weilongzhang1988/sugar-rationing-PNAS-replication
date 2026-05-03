/* =============================================================================
   08_figS8_liver_C22_vs_HCC.do
   PRODUCES:
     - SI Figure S8 panel (a): Liver/IHBD cancer (C22.0-C22.9, broad).
       Reuses the Cox specification from 02_fig1_cancer_event_study.do
       (251024_esd_Liver_cancer_coefplot.pdf).
     - SI Figure S8 panel (b): Hepatocellular carcinoma only (C22.0),
       Gompertz parametric hazard model
       -> FigS8b_HCC_C220_gompertz_coefplot.pdf
     - Pre-trend joint F-tests (liver C22 and rectum C20)

   INPUT:  ${UKB}/260422_sugar_cox_data.dta
   ============================================================================= */

capture log close
log using "${LOG}/08_figS8_liver_C22_vs_HCC.log", replace

use "${UKB}/260422_sugar_cox_data.dta", clear

* -------------------------------------------------------------------------
* PART A: C22 (broad), Cox PH — main outcome (panel a of Fig S8)
* -------------------------------------------------------------------------
stset diag_C22_time, id(eid) failure(diag_C22_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo c22_broad

coefplot c22_broad, omitted keep(*study*) vertical eform                     ///
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
    title("(a) Liver/IHBD cancer C22.0-C22.9", size(medium))                  ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS8a_C22_broad_coefplot.pdf", replace

* Pre-trend joint test (parallel-trends placebo)
testparm 1.study 2.study 3.study
local f_liver  = r(F)
local p_liver  = r(p)
local df_liver = r(df)

* -------------------------------------------------------------------------
* PART B: C22.0 only (HCC), Gompertz parametric hazard (panel b of Fig S8)
*         - replicates Zheng et al. (2025) outcome definition.
* -------------------------------------------------------------------------
stset diag_C220_time, id(eid) failure(diag_C220_event)
streg ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index,                      ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo hcc_gomp

coefplot hcc_gomp, omitted keep(*study*) vertical eform                      ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.00(0.5)5.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    title("(b) Hepatocellular carcinoma only (C22.0)", size(medium))         ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS8b_HCC_C220_gompertz_coefplot.pdf", replace

* -------------------------------------------------------------------------
* PART C: Pre-trend joint F-test summary (Liver C22, Rectum C20)
* -------------------------------------------------------------------------
stset diag_C20_Rectum_time, failure(diag_C20_Rectum_event)
stcox ib4.study i.male i.month_birth smoking college pca1-pca10 ///
      home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
testparm 1.study 2.study 3.study
local f_rectum  = r(F)
local p_rectum  = r(p)

di _n "============================================="
di    "  Pre-trend joint F-tests (studies 1-3 = 0)"
di    "============================================="
di    "  Liver cancer (C22):  F = `: di %5.2f `f_liver'',  p = `: di %5.3f `p_liver''"
di    "  Rectum cancer (C20): F = `: di %5.2f `f_rectum'', p = `: di %5.3f `p_rectum''"
di    "============================================="

log close
