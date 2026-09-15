# Replication package

**Paper:** Early-life sugar restriction causally reduces adult cancer incidence and slows biological ageing
**Authors:** Chen Zhu (China Agricultural University), Weilong Zhang (University of Cambridge)
**Submitted to:** *PNAS*, 2026
**Manuscript file:** `Research_report_sugar_rationing_260429_black.tex` (`Apps/Overleaf/Sugar_Paper_PANS/`)

---

## 1. Folder layout

```
PNAS_replication/
├── master.do                       # runs every script; user must set 2 globals
├── README.md                       # this file
├── data/
│   └── NFS/                        # National Food Survey (digitised by Gracner et al. 2024)
│       ├── Dataset_1000days_Fig_S1.dta        # calorie series (Fig S4)
│       ├── Dataset_1000days_NUTRITION_Fig_S3.dta  # food category panel (Fig S3)
│       ├── Dataset_1000days_SALES_Fig_S2.dta  # sugar/sweets sales
│       └── Dataset_1000days_foodafford.dta    # food prices (Fig S5)
├── results/
│   └── reply_to_cirillo/           # aggregate estimates behind the reply to Cirillo (2026)
├── code/
│   ├── reply_to_cirillo/           # reply to Cirillo (2026): Tables 1 and 2
│   │   ├── 01_trend_altcutoff.py   # Table 1 (Python; needs UKB data)
│   │   ├── 01_trend_altcutoff.do   # Table 1 (Stata twin)
│   │   └── 02_table2_NFS_harmonised.py  # Table 2 (NFS data only)
│   ├── main/                       # produces every figure and table in the main text
│   │   ├── 01_table1_balance.do
│   │   ├── 02_fig1_cancer_event_study.do
│   │   ├── 03_fig2_diet_event_study.do
│   │   ├── 04_fig3_bioaging_event_study.do
│   │   └── 05_table_S6_nutrient_intake.do
│   └── appendix/                   # produces every figure and table in the SI Appendix
│       ├── 06_figS3_S4_S5_NFS_consumption.do        # Figs S3, S4, S5
│       ├── 07_figS7_cardiometabolic.do              # Fig S7 + Table S3
│       ├── 08_figS8_liver_C22_vs_HCC.do             # Fig S8
│       ├── 09_KM_logrank_tableS8.do                 # Table S8 + KM curves
│       ├── 10_tableS9_S10_S11_S12_representativeness.do  # Tables S10, S11, S12
│       ├── 11_tableS15_S16_WebQ_IPW.do              # Tables S15, S16
│       ├── 12_tableS17_S18_confounder_LTL_GZMB.do   # Tables S17, S18
│       ├── 13_tableS9_england_only.do               # Table S9
│       ├── 15_multiple_testing_correction.do        # Bonferroni values inline in Table S2
│       └── 16_figS6_all_category_age50.do           # Fig S6
└── output/
    ├── figures/                    # PDF and EPS figures, named by paper figure number
    ├── tables/                     # LaTeX fragments, named by paper table number
    ├── intermediate/               # .gph panels, .csv counts, etc.
    └── logs/                       # one .log per do-file
```

---

## 2. Software requirements

* **Stata 17** or higher (uses post-Stata 16 syntax for `lpoly`, `coefplot`, `frame`).
* User-written packages (install once via `ssc install`):

  ```stata
  ssc install asdoc
  ssc install coefplot
  ssc install estout
  ssc install reghdfe
  ssc install ftools
  ```

---

## 3. Data access

The analytic UK Biobank derivative datasets (~300 MB each) are not redistributable
under the UK Biobank Material Transfer Agreement. They live in the same parent
folder as `PNAS_replication/`:

| File | Purpose |
| --- | --- |
| `251025_sugar_data_for_cox.dta`  | Cancer Cox sample (lung, breast, liver, prostate, rectum) |
| `251027_sugar_data_for_cox.dta`  | Cardiometabolic and Oxford WebQ dietary sample |
| `251029_sugar_data_for_cox.dta`  | Bioaging markers (LTL, Granzyme B), T1D placebo |
| `251109_sugar_data_for_KM_curves.dta` | Kaplan-Meier survival curves |
| `260422_sugar_cox_data.dta`      | C22 vs C22.0 sensitivity analyses |

Researchers wishing to reproduce these results must apply to UK Biobank in
their own right (see the data availability statement in the paper). The
National Food Survey data redistributed in `data/NFS/` are public and were
digitised by Gracner et al. (2024).

**Codebook note (NFS food-category panel).** In
`Dataset_1000days_NUTRITION_Fig_S3.dta` the six variables whose names end in
`_gday10` (`MilkCream_gday10`, `freshfruit_gday10`, `Vegetables_gday10`,
`bread_gday10`, `totalcereals_gday10`, `Meats_gday10`) are stored in **tens of
grams per day**; multiply by 10 to obtain g/day. All other `_gday` variables
are in grams per day. The variable label of `bread_gday10` omits the
"(in 10s)" flag carried by the other five; the name suffix is authoritative.
`totalcereals_gday10` already includes bread and flour.

---

## 4. How to run

1. Open `master.do` and set the two path globals near the top:

   ```stata
   global ROOT  "<path to>/PNAS_replication"
   global UKB   "<path to UK Biobank derivative dta files>"
   ```

2. Run `master.do`. Each script can also be executed in isolation provided
   the globals `ROOT`, `UKB`, `CODE`, `DATA`, `FIG`, `TAB`, `LOG`, `INTER`
   are defined first (the master file does this).

3. Outputs land in `output/figures/`, `output/tables/`, `output/intermediate/`,
   and `output/logs/`.

---

## 5. Mapping: paper exhibits to do-files

### Main paper

| Paper exhibit | Source script | Output file(s) |
|---|---|---|
| Table 1 — balance by exposure | `main/01_table1_balance.do` | `output/tables/Table1_balance.doc` |
| Figure 1 — cancer event study | `main/02_fig1_cancer_event_study.do` | `Fig1A..E_*_coefplot.pdf` (combined into `sugar_fig_cancer.png`) |
| Figure 2 — diet event study   | `main/03_fig2_diet_event_study.do` | `Fig2A..C_*_coefplot.pdf` (combined into `sugar_fig_consumption.png`) |
| Figure 3 — bioaging           | `main/04_fig3_bioaging_event_study.do` | `Fig3A_LTL_coefplot.pdf`, `Fig3B_GZMB_coefplot.pdf` (combined into `sugar_fig_ltl.png`) |

### SI Appendix

| SI exhibit | Source script | Output file(s) |
|---|---|---|
| Fig S3 — quarterly food consumption | `appendix/06_figS3_S4_S5_NFS_consumption.do` | `FigS3_other_consumption.eps/.pdf` |
| Fig S4 — calorie series         | `appendix/06_figS3_S4_S5_NFS_consumption.do` | `FigS4_suger_calories.eps/.pdf` |
| Fig S5 — food prices            | `appendix/06_figS3_S4_S5_NFS_consumption.do` | `FigS5_food_price.eps/.pdf` |
| Fig S6 — food intake by category | `appendix/16_figS6_all_category_age50.do` | `FigS6_all_category_age50.eps/.pdf` |
| Fig S7 — cardiometabolic event study | `appendix/07_figS7_cardiometabolic.do` | `FigS7A..E_*_coefplot.pdf` |
| Fig S8 — liver C22 vs C22.0     | `appendix/08_figS8_liver_C22_vs_HCC.do` | `FigS8a_C22_broad_coefplot.pdf`, `FigS8b_HCC_C220_gompertz_coefplot.pdf` |
| Table S2 — cancer HR table      | `main/02_fig1_cancer_event_study.do` | `tableS2_cancer_HR.tex` |
| Table S3 — cardiometabolic HR   | `appendix/07_figS7_cardiometabolic.do` | `tableS3_cardiometabolic_HR.tex` |
| Table S4 — bioaging coefficients | `main/04_fig3_bioaging_event_study.do` | `tableS4_bioaging_markers.tex` |
| Table S5 — food weight, HEI, FD  | `main/03_fig2_diet_event_study.do` | `tableS5_food_weight_HEI_diversity.tex` |
| Table S6 — nutrient intake       | `main/05_table_S6_nutrient_intake.do` | `tableS6_nutrient_intake.tex` |
| Table S8 — log-rank tests       | `appendix/09_KM_logrank_tableS8.do` | written to log; manually transcribed |
| Table S9 — England-only HR      | `appendix/13_tableS9_england_only.do` | `tableS9_cancer_HR_england_only.tex` |
| Table S10 — IPW HR              | `appendix/10_tableS9_S10_S11_S12_representativeness.do` | `tableS10_IPW_cancer.tex` |
| Table S11 — Survivorship HR     | `appendix/10_tableS9_S10_S11_S12_representativeness.do` | `tableS11_survival_cancer.tex` |
| Table S12 — Missingness         | `appendix/10_tableS9_S10_S11_S12_representativeness.do` | `tableS12_missingness.tex` |
| Table S15 — WebQ balance + selection | `appendix/11_tableS15_S16_WebQ_IPW.do` | `tableS15_WebQ_balance_selection.tex` |
| Table S16 — WebQ unweighted vs IPW   | `appendix/11_tableS15_S16_WebQ_IPW.do` | `tableS16_WebQ_diet_IPW.tex` |
| Table S17 — LTL confounder sensitivity | `appendix/12_tableS17_S18_confounder_LTL_GZMB.do` | `tableS17_LTL_confounder_sensitivity.tex` |
| Table S18 — Granzyme B confounder sensitivity | `appendix/12_tableS17_S18_confounder_LTL_GZMB.do` | `tableS18_GZMB_confounder_sensitivity.tex` |
| Table S2 inline Bonferroni p-values | `appendix/15_multiple_testing_correction.do` | (matrix output transcribed into Table S2) |

### Manually formatted tables (no replication code required)

* Table S1 — Disease outcome definitions (ICD-10 codes, from medical-records linkage)
* Table S7 — Sample distribution by birth year-month (descriptive summary)
* Table S13 — Disease outcome definitions, sample construction, estimation approach
* Table S14 — Effect-size benchmarking (drawn from external meta-analyses)
* Figure S1 — Food rationing timeline (image asset)
* Figure S2 — Flowchart of sample selection (image asset)

### Composite figures

`sugar_fig_cancer.png`, `sugar_fig_consumption.png`, and `sugar_fig_ltl.png`
shown in the main text are composite layouts assembled outside Stata
(Adobe Illustrator / Inkscape) from the individual coefficient plots emitted
by scripts 02, 03, and 04. The component PDFs are written to
`output/figures/` and can be combined exactly as in the published version.

---

## 6. Reply to Cirillo (PNAS Letter, 2026)

`code/reply_to_cirillo/` holds the code behind the two tables in our reply to
N. Cirillo, "Reassessing the causal interpretation of the UK sugar-rationing
natural experiment" (PNAS, 2026). `results/reply_to_cirillo/` holds the
aggregate output (hazard ratios, confidence intervals, test statistics, sample
and event counts); no individual-level UK Biobank data are included.

| Reply exhibit | Script | Output |
|---|---|---|
| Table 1 — HR per six months under rationing, three cutoffs, departure-from-linearity test | `reply_to_cirillo/01_trend_altcutoff.py` (or the Stata twin `01_trend_altcutoff.do`) | `results/reply_to_cirillo/trend_altcutoff_summary.csv`, `altcutoff_cells.csv` |
| Table 2 — NFS purchase changes, as reported vs harmonised units | `reply_to_cirillo/02_table2_NFS_harmonised.py` | `results/reply_to_cirillo/table2_NFS_harmonised.csv` |

Table 1 needs `251025_sugar_data_for_cox.dta` (UK Biobank, application 89068,
not redistributed; see section 3). Run the Python version with the folder
holding that file in the environment variable `UKB` (or as the first argument);
it uses Breslow ties and Lin-Wei cluster-robust standard errors by month of
birth, and reproduces SI Table S2 to three decimals. Table 2 runs from the
redistributed NFS data alone. The `k` column in the summary and cell files
indexes the cutoff: 0 = 26 September 1953 (published), 1 = end June 1953,
2 = end March 1953.

---

## 7. Citation

If you use this replication package, please cite both the paper and the
underlying data sources:

* Zhu, C. & Zhang, W. (2026). *Early-life sugar restriction causally reduces
  adult cancer incidence and slows biological ageing.* Submitted to *PNAS*.
* UK Biobank (application 89068). [https://www.ukbiobank.ac.uk](https://www.ukbiobank.ac.uk)
* Gracner, T. et al. (2024). *Exposure to sugar rationing in the first 1000
  days of life protected against chronic disease.* Science 386:1043-1048.
