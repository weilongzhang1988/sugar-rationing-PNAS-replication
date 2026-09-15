******************************************************************************
*  01_trend_altcutoff.do  --  Reply to Cirillo (PNAS Letter, 2026): Table 1
*
*  Stata twin of 01_trend_altcutoff.py. Same specification as SI Table S2
*  (Cox, cluster by month of birth). For each cancer and cutoff:
*    (a) categorical event-study cells (replication check of SI Table S2)
*    (b) HR per six months of the first 1,000 days under rationing
*        (continuous exposure, pre-period cells kept as dummies)
*    (c) HR per ordered exposure category (score trend)
*    (d) Wald test of departure from a linear gradient across the rationed
*        cells (equal successive increments in log HR), plus successive
*        differences
*  Cutoffs: k=0 published (26 Sep 1953); k=1 one quarter earlier (end June
*  1953, after the 17 May ration increase); k=2 two quarters earlier (end
*  March 1953, after sweets derationing 5 Feb and the 19 April increase).
*
*  Data: ${UKB}/251025_sugar_data_for_cox.dta (UK Biobank application 89068;
*        not redistributed). Requires the path globals set by master.do.
*  Output: ${TAB}/reply_cirillo_*.tex, ${INTER}/reply_cirillo_summary.csv
******************************************************************************
clear all
set more off
capture log close
log using "${LOG}/reply_cirillo_trend_altcutoff.log", replace text

use "${UKB}/251025_sugar_data_for_cox.dta", clear

* ---- quarter index rv (rv = quarter of birth + 41; 8 = 1951q4)
capture drop rv
gen rv = yearqbirth + 41
tab rv study

* ---- months of the first 1,000 days (conception to 24 months) spent under
*      rationing, monthly precision, for cutoff month C: E = C - (birth - 9),
*      clipped to [0, 33]
local C0 = tm(1953m9)
gen expm0 = `C0' - yearmobirth + 9
replace expm0 = 0  if expm0 < 0
replace expm0 = 33 if expm0 > 33
gen exp6_0 = expm0/6
label var exp6_0 "Exposure (6-month units), cutoff Sep 1953"

* pre-period cells kept as separate dummies; all rationed cells pooled onto
* the reference so that the slope captures them
gen preg0 = study
replace preg0 = 4 if study >= 5
gen score0 = max(study - 4, 0)          // 0 ref, 1 in-utero, ..., 5 +24m

* ---- alternative cutoffs: shift bins k quarters earlier
foreach k in 1 2 {
    gen rv`k' = rv + `k'
    gen study`k' = 1 if rv`k' > 24
    replace study`k' = 2 if inlist(rv`k', 23, 24)
    replace study`k' = 3 if inlist(rv`k', 21, 22)
    replace study`k' = 4 if inlist(rv`k', 19, 20)
    replace study`k' = 5 if inlist(rv`k', 16, 17, 18)
    replace study`k' = 6 if inlist(rv`k', 14, 15)
    replace study`k' = 7 if inlist(rv`k', 12, 13)
    replace study`k' = 8 if inlist(rv`k', 10, 11)
    replace study`k' = 9 if inlist(rv`k', 8, 9)
    replace study`k' = . if rv`k' < 8
    local Ck = `C0' - 3*`k'
    gen expm`k' = `Ck' - yearmobirth + 9
    replace expm`k' = 0  if expm`k' < 0
    replace expm`k' = 33 if expm`k' > 33
    gen exp6_`k' = expm`k'/6
    gen preg`k' = study`k'
    replace preg`k' = 4 if study`k' >= 5 & study`k' < .
    gen score`k' = max(study`k' - 4, 0) if study`k' < .
}
label define studylbl 1 "-27m" 2 "-21m" 3 "-15m" 4 "ref -9m" 5 "in utero" ///
    6 "+6m" 7 "+12m" 8 "+18m" 9 "+24m"
foreach v in study study1 study2 { label values `v' studylbl }
tab study study1
tab study study2

* ---- outcomes
local outcomes lung liver rectum prostate breast
local t_lung     diag_LungCancer_time
local e_lung     diag_LungCancer_event
local s_lung     "1"
local t_liver    diag_C22_time
local e_liver    diag_C22_event
local s_liver    "1"
local t_rectum   diag_C20_Rectum_time
local e_rectum   diag_C20_Rectum_event
local s_rectum   "1"
local t_prostate diag_C61_Prostate_time
local e_prostate diag_C61_Prostate_event
local s_prostate "male==1"
local t_breast   diag_C50_Breast_time
local e_breast   diag_C50_Breast_event
local s_breast   "male==0"

local covs i.month_birth smoking college pca1-pca10 home_longi home_lati Townsend_deprivation_index

* results matrix: rows = outcomes; per cutoff k the columns are
*   HR per 6 months, p ; HR per category step, p ; p nonlinearity ;
*   p pre-trend joint ; HR +24 cell ; p joint rationed
tempname R
matrix `R' = J(5, 24, .)
matrix rownames `R' = lung liver rectum prostate breast
local cn ""
foreach k in 0 1 2 {
    local cn `cn' HR6_k`k' p6_k`k' HRcat_k`k' pcat_k`k' pnonlin_k`k' ppre_k`k' HR24_k`k' pjoint_k`k'
}
matrix colnames `R' = `cn'

local i = 0
foreach o of local outcomes {
    local ++i
    local sexvar i.male
    if "`s_`o''" != "1" local sexvar ""
    stset `t_`o'', failure(`e_`o'')

    foreach k in 0 1 2 {
        local sv study
        if `k' > 0 local sv study`k'
        di _n(2) "==================== `o' | cutoff shift k=`k' ===================="
        * (a) categorical event study
        eststo cat_`o'_`k': stcox ib4.`sv' `sexvar' `covs' if `s_`o'', cluster(yearmobirth) nolog
        local j = 8*`k'
        capture matrix `R'[`i', `j'+7] = exp(_b[9.`sv'])
        testparm 5.`sv' 6.`sv' 7.`sv' 8.`sv' 9.`sv'
        matrix `R'[`i', `j'+8] = r(p)
        testparm 1.`sv' 2.`sv' 3.`sv'
        matrix `R'[`i', `j'+6] = r(p)
        * (d) departure from linearity in the ordered cells (equal increments)
        capture test (6.`sv' - 5.`sv' = 7.`sv' - 6.`sv') (7.`sv' - 6.`sv' = 8.`sv' - 7.`sv') (8.`sv' - 7.`sv' = 9.`sv' - 8.`sv')
        if _rc == 0 matrix `R'[`i', `j'+5] = r(p)
        lincom 5.`sv'
        lincom 6.`sv' - 5.`sv'
        lincom 7.`sv' - 6.`sv'
        lincom 8.`sv' - 7.`sv'
        capture lincom 9.`sv' - 8.`sv'

        * (b) linear trend, continuous months of exposure (6-month units)
        eststo lin_`o'_`k': stcox c.exp6_`k' ib4.preg`k' `sexvar' `covs' if `s_`o'', cluster(yearmobirth) nolog
        matrix `R'[`i', `j'+1] = exp(_b[exp6_`k'])
        matrix `R'[`i', `j'+2] = 2*normal(-abs(_b[exp6_`k']/_se[exp6_`k']))

        * (c) linear trend, ordered category score (0..5)
        eststo sc_`o'_`k': stcox c.score`k' ib4.preg`k' `sexvar' `covs' if `s_`o'', cluster(yearmobirth) nolog
        matrix `R'[`i', `j'+3] = exp(_b[score`k'])
        matrix `R'[`i', `j'+4] = 2*normal(-abs(_b[score`k']/_se[score`k']))
    }
}

di _n(2) "==================== SUMMARY MATRIX ===================="
matrix list `R', format(%9.4f)
matrix RES = `R'
preserve
clear
svmat RES, names(col)
gen outcome = ""
replace outcome = "lung"     in 1
replace outcome = "liver"    in 2
replace outcome = "rectum"   in 3
replace outcome = "prostate" in 4
replace outcome = "breast"   in 5
order outcome
export delimited using "${INTER}/reply_cirillo_summary.csv", replace
restore

* ---- tables: cell HRs under the three cutoffs; continuous trend at k=0
foreach k in 0 1 2 {
    local sv study
    if `k' > 0 local sv study`k'
    esttab cat_lung_`k' cat_liver_`k' cat_rectum_`k' cat_prostate_`k' cat_breast_`k' ///
        using "${TAB}/reply_cirillo_altcutoff_k`k'_cancer_HR.tex", replace eform keep(*`sv'*) ///
        b(%5.2f) ci(%5.2f) nostar mtitles("Lung" "Liver" "Rectum" "Prostate" "Breast") ///
        title("Cutoff shifted `k' quarter(s) earlier")
}
esttab lin_lung_0 lin_liver_0 lin_rectum_0 lin_prostate_0 lin_breast_0 ///
    using "${TAB}/reply_cirillo_trend_continuous_k0.tex", replace eform keep(exp6_0) b(%5.3f) ci(%5.3f) ///
    mtitles("Lung" "Liver" "Rectum" "Prostate" "Breast")

log close
exit
