/* =============================================================================
   07_figS7_cardiometabolic.do
   PRODUCES:
     - SI Figure S7 component panels (combined externally as fig:replication_exercise):
           Hypertension   -> 251024_esd_hypertension_coefplot.pdf
           Type 2 diabetes -> 251024_esd_T2D_coefplot.pdf
           Obesity        -> 251024_esd_Obesity_coefplot.pdf
           Heart failure  -> 251024_esd_Heart_failure_coefplot.pdf
           Type 1 diabetes (placebo) -> 251024_esd_T1D_coefplot.pdf
     - SI Table S3 (cardiometabolic hazard ratios) -> tableS3_cardiometabolic_HR.tex

   INPUT:  ${UKB}/251027_sugar_data_for_cox.dta  (hypertension, T2D, obesity, heart failure)
           ${UKB}/251029_sugar_data_for_cox.dta  (T1D, cataract placebo)
   ============================================================================= */

capture log close
log using "${LOG}/07_figS7_cardiometabolic.log", replace

* -------------------------------------------------------------------------
* I. Hypertension, T2D, Obesity, Heart failure (251027 sample)
* -------------------------------------------------------------------------
use "${UKB}/251027_sugar_data_for_cox.dta", clear

* Hypertension
stset diag_hyp_time, failure(diag_hyp_event)
streg ib4.study i.male i.month_birth p1647_i0_6 p1647_i0_5 ///
      pca1-pca10 home_longi home_lati Townsend_deprivation_index, ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo hypertension

coefplot hypertension, omitted keep(*study*) vertical eform                  ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.60(0.1)1.40, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS7A_hypertension_coefplot.pdf", replace

* Type 2 diabetes
gen     diag_dm2_time  = diab_dm_w3min
replace diag_dm2_time  = age_max if diab_dm_w3min == .
gen     diag_dm2_event = (diab_dm_w3min != .)

stset diag_dm2_time, failure(diag_dm2_event)
streg ib4.study i.male i.month_birth p1647_i0_6 p1647_i0_5 ///
      pca1-pca10 home_longi home_lati Townsend_deprivation_index, ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo T2D

coefplot T2D, omitted keep(*study*) vertical eform                           ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(red%40))                                   ///
    xline(4.5, lpattern(dash_dot) lcolor(gs10%60))                            ///
    ylabel(0.60(0.1)1.40, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS7B_T2D_coefplot.pdf", replace

* Obesity
gen     diag_obesity_time  = diag_obesitymin_prim
replace diag_obesity_time  = age_max if missing(diag_obesitymin_prim)
gen     diag_obesity_event = (diag_obesitymin_prim != .)

stset diag_obesity_time, failure(diag_obesity_event)
streg ib4.study i.male i.month_birth p1647_i0_6 p1647_i0_5 ///
      pca1-pca10 home_longi home_lati Townsend_deprivation_index, ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo obesity

coefplot obesity, omitted keep(*study*) vertical eform                       ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(red%40))                                   ///
    xline(4.5, lpattern(dash_dot) lcolor(gs10%60))                            ///
    ylabel(0.20(0.2)1.60, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS7C_Obesity_coefplot.pdf", replace

* Heart failure
stset timeToEvent2_I50_Heart_failure, failure(I50_Heart_failure)
streg ib4.study i.male i.month_birth p1647_i0_6 p1647_i0_5 ///
      pca1-pca10 home_longi home_lati Townsend_deprivation_index, ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo heart_failure

coefplot heart_failure, omitted keep(*study*) vertical eform                 ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(0.40(0.1)1.60, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS7D_Heart_failure_coefplot.pdf", replace

* -------------------------------------------------------------------------
* II. Type 1 diabetes placebo (251029 sample)
* -------------------------------------------------------------------------
use "${UKB}/251029_sugar_data_for_cox.dta", clear

stset diag_t1d_time, failure(diag_t1d_event)
streg ib4.study i.male i.college i.month_birth p1647_i0_6 p1647_i0_5 ///
      pca1-pca10 home_longi home_lati Townsend_deprivation_index, ///
      distribution(gompertz) nolog cluster(yearmobirth)
eststo T1D

coefplot T1D, omitted keep(*study*) vertical eform                           ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(1, lpattern(dash) lcolor(red%40))                                   ///
    xline(4.5, lpattern(dash_dot) lcolor(gs10%60))                            ///
    ylabel(0.20(0.3)2.00, angle(horizontal))                                  ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Hazard ratio (95% CI)")                                           ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/FigS7E_T1D_placebo_coefplot.pdf", replace

* -------------------------------------------------------------------------
* SI Table S3: Cardiometabolic hazard ratios (Hypertension, T2D, Obesity,
* Heart failure, Type 1 diabetes placebo)
* -------------------------------------------------------------------------
esttab hypertension T2D obesity heart_failure T1D ///
       using "${TAB}/tableS3_cardiometabolic_HR.tex", replace                 ///
       eform keep(*study*)                                                    ///
       cells("b(fmt(2) star) ci(par fmt(2))")                                 ///
       collabels("Hypertension" "Type 2 diabetes" "Obesity" "Heart failure" "Type 1 diabetes") ///
       varlabels(*study* "Exposure duration (months)")                        ///
       refcat(5.study "Reference: born post-rationing", nolabel)              ///
       star(* 0.10 ** 0.05 *** 0.01)                                          ///
       label booktabs nomtitles nonumbers                                     ///
       title("SI Table S3: Hazard ratios (95\% CI) for cardiometabolic disease incidence by exposure duration to sugar rationing") ///
       addnotes("Gompertz parametric hazard estimates with clustered standard errors by month of birth." ///
                "* p<0.10, ** p<0.05, *** p<0.01.")

log close
