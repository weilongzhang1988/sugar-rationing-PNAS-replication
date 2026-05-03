/* =============================================================================
   04_fig3_bioaging_event_study.do
   PRODUCES:
     - Figure 3 (main paper) component panels (combined externally into
       sugar_fig_ltl.png):
           Leukocyte telomere length -> Fig3A_LTL_coefplot.pdf
           Granzyme B               -> Fig3B_GZMB_coefplot.pdf
     - SI Table S4 (bioaging markers) -> tableS4_bioaging_markers.tex

   INPUT:  ${UKB}/251029_sugar_data_for_cox.dta
   ============================================================================= */

capture log close
log using "${LOG}/04_fig3_bioaging_event_study.log", replace

use "${UKB}/251029_sugar_data_for_cox.dta", clear

* -------------------------------------------------------------------------
* 1. Leukocyte telomere length (Z-adjusted T/S ratio)
* -------------------------------------------------------------------------
reg zadjustedtsloginstance0 ib4.study i.male i.month_birth pca1-pca10 age   ///
    home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo ltl

coefplot ltl, omitted keep(*study*) vertical                                 ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(0, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(-0.05(0.05)0.15, angle(horizontal))                                ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing", size(small))                    ///
    ytitle("Leukocyte telomere length (Z-adjusted)", size(small))             ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig3A_LTL_coefplot.pdf", replace

* -------------------------------------------------------------------------
* 2. Granzyme B (log2 NPX)
* -------------------------------------------------------------------------
reg gzmbgranzymeb ib4.study i.male i.month_birth pca1-pca10 age              ///
    home_longi home_lati Townsend_deprivation_index, cluster(yearmobirth)
eststo gzmb

coefplot gzmb, omitted keep(*study*) vertical                                ///
    recast(connected) msymbol(O) mcolor(eltblue) lcolor(eltblue)             ///
    ciopts(recast(rarea rline) color(eltblue%30) lcolor(eltblue%60))         ///
    yline(0, lpattern(dash) lcolor(gs10%60))                                  ///
    xline(4.5, lpattern(dash_dot) lcolor(red%40))                             ///
    ylabel(-0.50(0.10)0.30, angle(horizontal))                                ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing", size(small))                    ///
    ytitle("Granzyme B (log{sub:2} NPX)", size(small))                        ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig3B_GZMB_coefplot.pdf", replace

* -------------------------------------------------------------------------
* SI Table S4: Bioaging marker event-study estimates
* -------------------------------------------------------------------------
esttab ltl gzmb using "${TAB}/tableS4_bioaging_markers.tex", replace          ///
    keep(*study*)                                                             ///
    cells("b(fmt(3) star) ci(par fmt(3))")                                    ///
    label nomtitles nonumbers booktabs                                        ///
    collabels("Leukocyte telomere length (Z-adjusted)" "Granzyme B")          ///
    title("SI Table S4: Coefficient estimates (95\% CI) for biological ageing markers by exposure duration to sugar rationing") ///
    addnotes("OLS regressions with clustered standard errors by month of birth." ///
             "Each column reports exposure-duration cohort effects relative to post-rationing births." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

log close
