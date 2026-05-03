/* =============================================================================
   05_table_S6_nutrient_intake.do
   PRODUCES:
     - SI Table S6 (Sugars, Fibre, Starch, Fat, Protein by exposure duration)
       -> tableS6_nutrient_intake.tex

   INPUT:  ${UKB}/251027_sugar_data_for_cox.dta
   ============================================================================= */

capture log close
log using "${LOG}/05_table_S6_nutrient_intake.log", replace

use "${UKB}/251027_sugar_data_for_cox.dta", clear

eststo sugars : reg avg_totalsugars        ib4.study male age college pca1-pca10 ///
                home_longi home_lati Townsend_deprivation_index i.birth_month, ///
                cluster(yearmobirth)
eststo fibre  : reg avg_englystdietaryfibre ib4.study male age college pca1-pca10 ///
                home_longi home_lati Townsend_deprivation_index i.birth_month, ///
                cluster(yearmobirth)
eststo starch : reg avg_starch              ib4.study male age college pca1-pca10 ///
                home_longi home_lati Townsend_deprivation_index i.birth_month, ///
                cluster(yearmobirth)
eststo fat    : reg avg_fat                 ib4.study male age college pca1-pca10 ///
                home_longi home_lati Townsend_deprivation_index i.birth_month, ///
                cluster(yearmobirth)
eststo protein: reg avg_protein             ib4.study male age college pca1-pca10 ///
                home_longi home_lati Townsend_deprivation_index i.birth_month, ///
                cluster(yearmobirth)

esttab sugars fibre starch fat protein using "${TAB}/tableS6_nutrient_intake.tex", ///
    keep(*study*) replace                                                     ///
    cells("b(fmt(3) star) ci(par fmt(3))")                                    ///
    label nomtitles nonumbers booktabs                                        ///
    collabels("Sugars" "Fibre" "Starch" "Fat" "Protein")                      ///
    title("SI Table S6: Coefficient estimates (95\% CI) for nutrient intake by exposure duration to sugar rationing") ///
    addnotes("OLS regressions with clustered standard errors by month of birth." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

log close
