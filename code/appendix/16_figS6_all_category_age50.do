/* =============================================================================
   16_figS6_all_category_age50.do
   PRODUCES:
     - SI Figure S6 (Impact of rationing exposure on food intake by category)
       -> ${FIG}/FigS6_all_category_age50.eps and .pdf

   This script reproduces the manually entered coefficients and standard
   errors from the IPW-reweighted dietary regression in
   `11_tableS15_S16_WebQ_IPW.do` (column avg_totalsugars_ipw and the
   pre-revision unweighted estimates for protein, fat, starch, fibre).
   No raw data are required: the do-file builds the figure data inline.
   ============================================================================= */

capture log close
log using "${LOG}/16_figS6_all_category_age50.log", replace

clear
* Step 1: Input coefficients and SEs
* ---------------------------------------------------------------------
* totalsugars column updated to IPW-reweighted estimates from
* 260421_consumption.log (primary specification after referee R1.4
* response). Other nutrient columns retain their pre-revision values
* and should be re-run under IPW for full consistency before publication.
* ---------------------------------------------------------------------
input period protein se_protein fat se_fat totalsugars se_totalsugars starch se_starch fibre se_fibre
1  -0.5234 0.9081  -0.0498 1.0674  -1.9028 1.2583  -0.4113 1.6084  -0.2369 0.2377
2   0.1395 0.6407  -0.8535 0.7532  -1.3317 1.3697  -0.0086 1.1349  -0.1019 0.1678
3  -0.6641 0.7194  -0.6404 0.8456  -2.7169 1.5169  -0.7223 1.2741   0.0570 0.1883
4   0       0        0      0        0       0        0       0        0       0
5  -1.2683 0.6360  -1.5908 0.7476  -3.6464 1.0359  -2.1739 1.1265  -0.1442 0.1665
6  -0.1841 0.6874  -0.6443 0.8080  -2.6015 1.4151  -0.1909 1.2174   0.2264 0.1910
7  -0.7417 0.7339  -0.3742 0.8627  -3.8648 1.1303  -1.0820 1.2998   0.0450 0.1921
8  -0.2884 0.7551   0.2469 0.8876  -3.3684 1.4279   0.0426 1.3374   0.1921 0.1977
9  -0.0464 0.8116   0.6031 0.9541  -3.4388 1.3607  -1.1389 1.3935   0.3831 0.2125
end

* Step 2: Create upper and lower bounds
foreach var in protein fat totalsugars starch fibre {
    gen ub_`var' = `var' + 1.96*se_`var'
    gen lb_`var' = `var' - 1.96*se_`var'
}

* Step 3: p-values and significance stars
foreach var in protein fat totalsugars starch fibre {
    gen pval_`var' = 2*normal(-abs(`var'/se_`var'))
    gen star_`var' = cond(pval_`var'<0.01,"***", cond(pval_`var'<0.05,"**", cond(pval_`var'<0.1,"*","")))
}

* Step 4: Shift periods to avoid overlap
gen period_protein     = period - 0.2
gen period_fat         = period - 0.1
gen period_totalsugars = period
gen period_starch      = period + 0.1
gen period_fibre       = period + 0.2

* Step 5: Labels
label define per 1 "-21" 2 "-15" 3 "-9" 4 "0" 5 "In-utero" 6 "In-utero+6" 7 "In-utero+12" 8 "In-utero+18" 9 "In-utero+24"
label values period per
label values period_protein per
label values period_fat per
label values period_totalsugars per
label values period_starch per
label values period_fibre per

* Step 6: Plot with stars at bottom of CI
twoway /// 
    (rcap ub_protein lb_protein period_protein, lcolor(blue%60)) ///
    (line protein period_protein, lcolor(blue%60) lpattern(solid) lwidth(medthin)) ///
    (scatter protein period_protein, msymbol(O) msize(small) mcolor(blue%60)) ///
    (scatter lb_protein period_protein, msymbol(i) mlabel(star_protein) mlabpos(6) mlabcolor(blue%60) mlabsize(small)) ///
    (rcap ub_fat lb_fat period_fat, lcolor(red%60)) ///
    (line fat period_fat, lcolor(red%60) lpattern(dash) lwidth(medthin)) ///
    (scatter fat period_fat, msymbol(O) msize(small) mcolor(red%60)) ///
    (scatter lb_fat period_fat, msymbol(i) mlabel(star_fat) mlabpos(6) mlabcolor(red%60) mlabsize(small)) ///
    (rcap ub_totalsugars lb_totalsugars period_totalsugars, lcolor(green%60)) ///
    (line totalsugars period_totalsugars, lcolor(green%60) lpattern(solid) lwidth(thick)) ///  <-- bold main line
    (scatter totalsugars period_totalsugars, msymbol(O) msize(medium) mcolor(green%60)) ///
    (scatter lb_totalsugars period_totalsugars, msymbol(i) mlabel(star_totalsugars) mlabpos(6) mlabcolor(green%60) mlabsize(small)) ///
    (rcap ub_starch lb_starch period_starch, lcolor(orange%60)) ///
    (line starch period_starch, lcolor(orange%60) lpattern(shortdash) lwidth(medthin)) ///
    (scatter starch period_starch, msymbol(O) msize(small) mcolor(orange%60)) ///
    (scatter lb_starch period_starch, msymbol(i) mlabel(star_starch) mlabpos(6) mlabcolor(orange%60) mlabsize(small)) ///
    (rcap ub_fibre lb_fibre period_fibre, lcolor(purple%60)) ///
    (line fibre period_fibre, lcolor(purple%60) lpattern(dot) lwidth(thin)) ///
    (scatter fibre period_fibre, msymbol(O) msize(small) mcolor(purple%60)) ///
    (scatter lb_fibre period_fibre, msymbol(i) mlabel(star_fibre) mlabpos(6) mlabcolor(purple%60) mlabsize(small)), ///
    yline(0, lcolor(gs8) lpattern(shortdash)) ///
	xline(4.5, lcolor(gs2) lpattern(shortdash) lwidth(thin)) ///
    xtitle("Exposure Window Relative to Rationing", size(small)) ///
    ytitle("Impact on Food Intake (grams per day)", size(small)) ///
	ylabel(-8.0(2.0)2.0, labsize(small)) ///
    xlabel(1 2 3 4 5 6 7 8 9, valuelabel labsize(vsmall)) ///
    legend(order(2 "Protein" 6 "Fat" 10 "Sugars" 14 "Starch" 18 "Dietary fibre") ///
           cols(3) size(small) region(lcolor(none))) ///
    graphregion(color(white)) bgcolor(white) scheme(s2mono)

	
	
graph export "${FIG}/FigS6_all_category_age50.eps", as(eps) replace
graph export "${FIG}/FigS6_all_category_age50.pdf", as(pdf) replace

log close
