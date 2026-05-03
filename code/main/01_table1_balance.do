/* =============================================================================
   01_table1_balance.do
   PRODUCES:
     - Table 1 (main paper): Characteristics of survey participants by exposure
                             to rationing.
   INPUT:  ${UKB}/251029_sugar_data_for_cox.dta
   OUTPUT: ${TAB}/Table1_balance.doc        (raw t-test output via asdoc)
           ${LOG}/01_table1_balance.log
   ============================================================================= */

capture log close
log using "${LOG}/01_table1_balance.log", replace

use "${UKB}/251029_sugar_data_for_cox.dta", clear

* Define the rationing exposure indicator used in the paper
*   before == 0 : never-rationed (study 1-4, born July 1954 onwards)
*   before == 1 : rationed       (study 5-9, born before June 1954)
gen     before = 0 if study > 0 & study < 5
replace before = 1 if study > 4 & study < 10

rename (p1647_i0_2 p1647_i0_6 p1647_i0_5 p21000_i0_17) ///
       (England    Wales      Scotland   white)

* Two-sample t-tests with unequal variances; output to a Word file using asdoc
local savepath "${TAB}/Table1_balance.doc"

asdoc ttest age_max, by(before) unequal replace save(`savepath')

local variable_list2 male England Wales Scotland white college smoking bmi   ///
    fa_Hypertension  fa_Diabetes  fa_Heart_disease fa_Stroke                  ///
    fa_Prostate_cancer fa_Breast_cancer fa_Bowel_cancer fa_Lung_cancer        ///
    mo_Hypertension  mo_Diabetes  mo_Heart_disease mo_Stroke                  ///
    mo_Prostate_cancer mo_Breast_cancer mo_Bowel_cancer mo_Lung_cancer

foreach var2 of local variable_list2 {
    asdoc ttest `var2', by(before) unequal rowappend save(`savepath')
}

log close
