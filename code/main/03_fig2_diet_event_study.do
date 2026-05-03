/* =============================================================================
   03_fig2_diet_event_study.do
   PRODUCES:
     - Figure 2 (main paper) component panels (combined externally into
       sugar_fig_consumption.png):
           Healthy Eating Index   -> Fig2A_HEI_coefplot.pdf
           Food Diversity Score   -> Fig2B_FoodDiversity_coefplot.pdf
           Food Weight (g/day)    -> Fig2C_FoodWeight_coefplot.pdf
     - SI Table S5 (food weight, HEI, food diversity event-study coefficients)
       -> tableS5_food_weight_HEI_diversity.tex

   INPUT:  ${UKB}/251027_sugar_data_for_cox.dta
           (Oxford WebQ subsample with HEI_score, Food_Diversity_score,
            avg_foodweight)
   ============================================================================= */

capture log close
log using "${LOG}/03_fig2_diet_event_study.log", replace

use "${UKB}/251027_sugar_data_for_cox.dta", clear

* -------------------------------------------------------------------------
* 1. Healthy Eating Index (HEI)
* -------------------------------------------------------------------------
eststo HEI:  reg HEI_score ib4.study male age college pca1-pca10           ///
                home_longi home_lati i.birth_month, cluster(yearmobirth)

coefplot HEI, omitted keep(*study*) vertical                                 ///
    recast(connected) msymbol(O) color(navy) ciopts(color(navy%50))          ///
    xline(4.5 5.5, lpattern(dash_dot) lcolor(black%30))                       ///
    yline(0, lpattern(dash) lcolor(red%45))                                   ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Coefficient (HEI points)")                                        ///
    title("Healthy Eating Index", size(medium))                               ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig2A_HEI_coefplot.pdf", replace

* -------------------------------------------------------------------------
* 2. Food Diversity Score
* -------------------------------------------------------------------------
eststo FD:   reg Food_Diversity_score ib4.study male age college pca1-pca10 ///
                home_longi home_lati i.birth_month, cluster(yearmobirth)

coefplot FD, omitted keep(*study*) vertical                                  ///
    recast(connected) msymbol(D) color(maroon) ciopts(color(maroon%50))       ///
    xline(4.5 5.5, lpattern(dash_dot) lcolor(black%30))                       ///
    yline(0, lpattern(dash) lcolor(red%45))                                   ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Coefficient (food groups)")                                       ///
    title("Food Diversity Score", size(medium))                               ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig2B_FoodDiversity_coefplot.pdf", replace

* -------------------------------------------------------------------------
* 3. Average food weight
* -------------------------------------------------------------------------
eststo FW:   reg avg_foodweight ib4.study male age college pca1-pca10        ///
                home_longi home_lati i.birth_month, cluster(yearmobirth)

coefplot FW, omitted keep(*study*) vertical                                  ///
    recast(connected) msymbol(O) color(forest_green) ciopts(color(forest_green%50)) ///
    xline(4.5 5.5, lpattern(dash_dot) lcolor(black%30))                       ///
    yline(0, lpattern(dash) lcolor(red%45))                                   ///
    xlabel(1 "-27" 2 "-21" 3 "-15" 4 "-9"                                     ///
           5 "In-utero" 6 "In-utero+6" 7 "In-utero+12"                        ///
           8 "In-utero+18" 9 "In-utero+24", labsize(vsmall))                  ///
    xtitle("Months of exposure to rationing")                                 ///
    ytitle("Coefficient (g/day)")                                             ///
    title("Average Food Weight", size(medium))                                ///
    graphregion(color(white)) scheme(s1color)                                 ///
    base mlabposition(1)
graph export "${FIG}/Fig2C_FoodWeight_coefplot.pdf", replace

* -------------------------------------------------------------------------
* SI Table S5: Food weight, HEI, and food diversity by exposure duration
* -------------------------------------------------------------------------
esttab FW HEI FD using "${TAB}/tableS5_food_weight_HEI_diversity.tex", ///
    keep(*study*) replace                                                     ///
    cells("b(fmt(3) star) ci(par fmt(3))")                                    ///
    label nomtitles nonumbers booktabs                                        ///
    collabels("Food weight (g/day)" "Healthy eating index" "Food diversity")  ///
    title("SI Table S5: Coefficient estimates (95\% CI) for nutrient intake by exposure duration to sugar rationing") ///
    addnotes("OLS regressions with clustered standard errors by month of birth." ///
             "Each column reports estimates for exposure-duration cohorts relative to post-rationing births." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

log close
