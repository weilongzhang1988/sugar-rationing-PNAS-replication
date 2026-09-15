******************************************************************************************
*   master.do
*   Replication package for:
*   "Early-life sugar restriction causally reduces adult cancer incidence and
*    slows biological ageing"
*    Zhu & Zhang (2026), submitted to PNAS.

*   This master file sets all paths and runs every numbered do-file in
*   PNAS_replication/code/. Each do-file is self-contained and may also be
*   run individually; the path globals defined below are required.

*   Usage:
*       1. Open Stata 17+ (commands use post-Stata 16 syntax).
*       2. Edit the two globals "ROOT" and "UKB" below.
*       3. Run this file.

*   Outputs:
*       - Tables: ${ROOT}/output/tables/*.tex
*       - Figures: ${ROOT}/output/figures/*.pdf, *.eps
*       - Logs:    ${ROOT}/output/logs/*.log

*   Required Stata user-written packages:
*       asdoc, coefplot, estout, reghdfe, ftools, rwolf2 (optional for SI S6)
*       ssc install asdoc
*       ssc install coefplot
*       ssc install estout
*       ssc install reghdfe
*       ssc install ftools
*       ssc install rwolf2
******************************************************************************************

clear all
set more off
capture log close _all
version 17

* -------------------------------------------------------------------------
* USER-CONFIGURABLE PATHS
* -------------------------------------------------------------------------
* ROOT: location of this PNAS_replication folder
global ROOT  "E:/CamUni SHSS-SAH Team Dropbox/weilong zhang/UK biobank/PNAS_replication"

* UKB: location of the large UK Biobank derivative .dta files supplied
*      under UK Biobank application 89068. These files are NOT redistributable
*      and remain in their original folder. They include:
*           251025_sugar_data_for_cox.dta   (cancer Cox sample)
*           251027_sugar_data_for_cox.dta   (cardiometabolic / dietary sample)
*           251029_sugar_data_for_cox.dta   (bioaging / placebo sample)
*           251109_sugar_data_for_KM_curves.dta (Kaplan-Meier sample)
*           260422_sugar_cox_data.dta       (C22 vs C22.0 sensitivity)
global UKB   "E:/CamUni SHSS-SAH Team Dropbox/weilong zhang/UK biobank"

* Derived paths used by the child do-files
global CODE   "${ROOT}/code"
global DATA   "${ROOT}/data"
global FIG    "${ROOT}/output/figures"
global TAB    "${ROOT}/output/tables"
global LOG    "${ROOT}/output/logs"
global INTER  "${ROOT}/output/intermediate"

* -------------------------------------------------------------------------
* RUN ORDER
* -------------------------------------------------------------------------

* --- Main paper -----------------------------------------------------------
do "${CODE}/main/01_table1_balance.do"          // Table 1
do "${CODE}/main/02_fig1_cancer_event_study.do" // Figure 1 + SI Table S2
do "${CODE}/main/03_fig2_diet_event_study.do"   // Figure 2 + SI Table S5
do "${CODE}/main/04_fig3_bioaging_event_study.do" // Figure 3 + SI Table S4
do "${CODE}/main/05_table_S6_nutrient_intake.do"  // SI Table S6 (sugars, fibre, ...)

* --- Appendix figures -----------------------------------------------------
do "${CODE}/appendix/06_figS3_S4_S5_NFS_consumption.do" // SI Figs S3, S4, S5
do "${CODE}/appendix/07_figS7_cardiometabolic.do"       // SI Fig S7 + SI Table S3
do "${CODE}/appendix/08_figS8_liver_C22_vs_HCC.do"      // SI Fig S8 (panels a,b)
do "${CODE}/appendix/09_KM_logrank_tableS8.do"          // SI Table S8 (KM + log-rank)
do "${CODE}/appendix/16_figS6_all_category_age50.do"    // SI Fig S6 (food intake by category)

* --- Reply to Cirillo (PNAS Letter, 2026) --------------------------------
* Not part of the paper's exhibits. The two Python scripts in
* code/reply_to_cirillo/ produce the reply's Tables 1 and 2 (see README, section 6).

* --- Appendix tables ------------------------------------------------------
do "${CODE}/appendix/13_tableS9_england_only.do"        // SI Table S9
do "${CODE}/appendix/10_tableS9_S10_S11_S12_representativeness.do" // SI Tables S10, S11, S12
do "${CODE}/appendix/11_tableS15_S16_WebQ_IPW.do"       // SI Tables S15, S16
do "${CODE}/appendix/12_tableS17_S18_confounder_LTL_GZMB.do" // SI Tables S17, S18

* --- Inline Bonferroni p-values for SI Table S2 --------------------------
do "${CODE}/appendix/15_multiple_testing_correction.do" // values transcribed into SI Table S2

di as result _n "==========================================================="
di as result    "Replication completed. Inspect ${ROOT}/output/"
di as result    "==========================================================="
