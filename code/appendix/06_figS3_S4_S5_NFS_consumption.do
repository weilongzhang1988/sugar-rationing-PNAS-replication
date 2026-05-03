/* =============================================================================
   06_figS3_S4_S5_NFS_consumption.do
   PRODUCES (using digitised National Food Survey data from Gracner et al. 2024):
     - SI Figure S3 panels: Fats, Dairy, Fruits & Vegetables, Cereals, Meats
     - SI Figure S4 (calorie consumption: total vs sugar)
     - SI Figure S5 (food prices vs CPI / vs wage rate)

   INPUT:  ${DATA}/NFS/Dataset_1000days_Fig_S1.dta       (calories)
           ${DATA}/NFS/Dataset_1000days_NUTRITION_Fig_S3.dta (consumption sub-categories)
           ${DATA}/NFS/Dataset_1000days_SALES_Fig_S2.dta  (sugar/sweets sales)
           ${DATA}/NFS/Dataset_1000days_foodafford.dta    (food affordability)
   ============================================================================= */

capture log close
log using "${LOG}/06_figS3_S4_S5_NFS_consumption.log", replace

* -------------------------------------------------------------------------
* SI Figure S4 (referred to as fig:sugar_pattern_a in appendix.tex):
* Calorie consumption between 1950 and 1957 (total vs sugar).
* Original NFS digitisation: Dataset_1000days_Fig_S1.dta
* -------------------------------------------------------------------------
preserve
    use "${DATA}/NFS/Dataset_1000days_Fig_S1.dta", clear

    twoway                                                                    ///
    (line Kcal_diff           qdate, lcolor(cranberry) lwidth(medthick))      ///
    (line Sugars_gdaycal_diff qdate, lcolor(blue) lpattern(longdash) lwidth(medthick)), ///
        ylabel(-100(50)200, grid)                                             ///
        xtitle("Quarter")                                                     ///
        ytitle("Changes in calories (vs 1950q1-1953q2)")                      ///
        xlabel(`=yq(1950,1)'(4)`=yq(1957,1)', format(%tq) angle(45))          ///
        legend(order(1 "Total Calories" 2 "Sugar Calories") rows(2)           ///
               position(10) ring(0) region(lstyle(none)) size(small))        ///
        xline(`=yq(1953,3)', lpattern(dash) lcolor(gs8))                     ///
        text(180 `=yq(1953,3)' "End of sugar rationing", size(small) place(e)) ///
        scheme(s1color)
    graph export "${FIG}/FigS4_suger_calories.eps", as(eps) replace
    graph export "${FIG}/FigS4_suger_calories.pdf", as(pdf) replace
restore

* (Sales figures from Dataset_1000days_SALES_Fig_S2 are not used in the
*  published appendix and have been removed; SI Figure S2 in the appendix
*  is the sample-selection flowchart, an image asset.)

* -------------------------------------------------------------------------
* SI Figure S3: Quarterly consumption of fats, dairy, fruits & vegetables,
* cereals, and meat (1952q1-1958q4).
* -------------------------------------------------------------------------
preserve
    use "${DATA}/NFS/Dataset_1000days_NUTRITION_Fig_S3.dta", clear

    * (A) Fats
    twoway line butter_gday    qdate, lp(solid)    lwidth(medthin)            ///
        xtitle("") ytitle("")                                                  ///
        xlabel(-32(4)-4, valuelabel labsize(small) angle(45))                 ///
        xaxis(1 2)                                                            ///
        xla(-23 "Butter rationing ends 5/1954",                               ///
            labcolor(black) grid noticks axis(2)                              ///
            glwidth(thin) glcolor(blue) labsize(vsmall))                      ///
        xtitle("", axis(2)) xlabel(, nogrid) ylabel(0(20)60, nogrid)          ///
        legend(region(lcolor(none) fcolor(none))                              ///
               order(1 "Butter" 2 "Margarine" 3 "Lard"))                      ///
     || line margarine_gday qdate, lwidth(medthin) lpattern(solid)            ///
     || line lard_gday      qdate, lwidth(medthin) lpattern(longdash)         ///
        legend(ring(0) position(10) size(vsmall) row(1)) scheme(s1color)      ///
        title("(A) Fats consumption", size(small))                            ///
        saving("${INTER}/figS3_a", replace)
    graph export "${FIG}/FigS3a_Fats.pdf", replace

    * (B) Dairy
    twoway line MilkCream_gday10 qdate, lp(solid) lwidth(medthin)             ///
        xtitle("") ytitle("")                                                  ///
        xlabel(-32(4)-4, valuelabel labsize(small) angle(45))                 ///
        xaxis(1 2)                                                            ///
        xla(-23 "Cheese rationing ends 5/1954",                               ///
            labcolor(black) grid noticks axis(2)                              ///
            glwidth(thin) glcolor(blue) labsize(vsmall))                      ///
        xtitle("", axis(2)) xlabel(, nogrid) ylabel(0(20)80, nogrid)          ///
        legend(region(lcolor(none) fcolor(none))                              ///
               order(1 "Milk" 2 "Welfare milk" 3 "Cheese"))                   ///
     || line lmilkwelfareandschoolpt_gday qdate, lwidth(medthin) lpattern(solid)  ///
     || line Cheese_gday qdate, lwidth(medthin) lpattern(longdash)            ///
        legend(ring(0) position(10) size(vsmall) row(1)) scheme(s1color)      ///
        title("(B) Dairy consumption", size(small))                           ///
        saving("${INTER}/figS3_b", replace)
    graph export "${FIG}/FigS3b_Dairy.pdf", replace

    * (C) Fruits and vegetables
    twoway line freshfruit_gday10 qdate, lp(solid) lwidth(medthin)            ///
        xtitle("") ytitle("")                                                  ///
        xlabel(-32(4)-4, valuelabel labsize(small) angle(45))                 ///
        xaxis(1 2)                                                            ///
        xla(-22 "Fruits rationing ends 7/1954",                               ///
            labcolor(black) grid noticks axis(2)                              ///
            glwidth(thin) glcolor(blue) labsize(vsmall))                      ///
        xtitle("", axis(2)) xlabel(, nogrid) ylabel(0(20)60, nogrid)          ///
        legend(region(lcolor(none) fcolor(none))                              ///
               order(1 "Fresh fruits" 2 "Canned fruit" 3 "Vegetables"))       ///
     || line otherfruit_gday qdate, lwidth(medthin) lpattern(solid)           ///
     || line Vegetables_gday10 qdate, lwidth(medthin) lpattern(longdash)      ///
        legend(ring(0) position(10) size(vsmall) row(1)) scheme(s1color)      ///
        title("(C) Fruits and vegetables consumption", size(small))           ///
        saving("${INTER}/figS3_c", replace)
    graph export "${FIG}/FigS3c_Fruits_vegetables.pdf", replace

    * (D) Cereal
    twoway line bread_gday10 qdate, lp(solid) lwidth(medthin)                 ///
        xtitle("") ytitle("")                                                  ///
        xlabel(-32(4)-4, valuelabel labsize(small) angle(45))                 ///
        xaxis(1 2)                                                            ///
        xla(-22 "Cereal rationing ends 7/1954",                               ///
            labcolor(black) grid noticks axis(2)                              ///
            glwidth(thin) glcolor(blue) labsize(vsmall))                      ///
        xtitle("", axis(2)) xlabel(, nogrid) ylabel(0(20)60, nogrid)          ///
        legend(region(lcolor(none) fcolor(none))                              ///
               order(1 "Bread" 2 "Flour" 3 "Cereals"))                        ///
     || line flour_gday qdate, lwidth(medthin) lpattern(solid)                ///
     || line totalcereals_gday10 qdate, lwidth(medthin) lpattern(longdash)    ///
        legend(ring(0) position(10) size(vsmall) row(1)) scheme(s1color)      ///
        title("(D) Cereal consumption", size(small))                          ///
        saving("${INTER}/figS3_d", replace)
    graph export "${FIG}/FigS3d_Cereals.pdf", replace

    * (E) Meats
    twoway line Meats_gday10 qdate, lp(solid) lwidth(medthin)                 ///
        xtitle("") ytitle("")                                                  ///
        xlabel(-32(4)-4, valuelabel labsize(small) angle(45))                 ///
        xaxis(1 2)                                                            ///
        xla(-22 "Meat rationing ends 7/1954",                                 ///
            labcolor(black) grid noticks axis(2)                              ///
            glwidth(thin) glcolor(blue) labsize(vsmall))                      ///
        xtitle("", axis(2)) xlabel(, nogrid) ylabel(0(20)60, nogrid)          ///
        legend(region(lcolor(none) fcolor(none))                              ///
               order(1 "Other meats" 2 "Fish" 3 "Bacon and ham"))             ///
     || line Fish_gday qdate                                                  ///
     || line baconandhamuncooked_gday qdate, lwidth(medthin) lpattern(solid)  ///
        legend(ring(0) position(10) size(vsmall) row(1)) scheme(s1color)      ///
        title("(E) Meat consumption", size(small))                            ///
        saving("${INTER}/figS3_e", replace)
    graph export "${FIG}/FigS3e_Meats.pdf", replace

    * Combined panel S3 (other_conusmption.eps in original appendix)
    graph combine "${INTER}/figS3_a.gph" "${INTER}/figS3_b.gph"               ///
                  "${INTER}/figS3_c.gph" "${INTER}/figS3_d.gph"               ///
                  "${INTER}/figS3_e.gph",                                     ///
        col(2) imargin(2 2 2 2)                                               ///
        graphregion(color(white)) scheme(s1color) xsize(6) ysize(8)
    graph export "${FIG}/FigS3_other_consumption.eps", as(eps) replace
    graph export "${FIG}/FigS3_other_consumption.pdf", as(pdf) replace
restore

* -------------------------------------------------------------------------
* SI Figure S5: Food prices vs CPI / vs wage rate (1952q1-1956q4).
* -------------------------------------------------------------------------
preserve
    use "${DATA}/NFS/Dataset_1000days_foodafford.dta", clear
    tsset qdate

    twoway connected rfood_priceq3 qdate, lp(solid) lwidth(medthin)           ///
        xtitle("") ytitle("Price Index (1953q3 = 100)")                       ///
        xlabel(, valuelabel labsize(small))                                   ///
        legend(order(1 "Food prices/CPI" 2 "Food prices/wage rate"))          ///
     || connected rfoodprice_wageq3 qdate, lwidth(medthin) lpattern(longdash) ///
        ylabel(0(10)130) legend(ring(0) position(6) size(normal) row(1))      ///
        scheme(s1color)
    graph export "${FIG}/FigS5_food_price.eps", as(eps) replace
    graph export "${FIG}/FigS5_food_price.pdf", as(pdf) replace
restore

log close
