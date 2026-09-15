"""
02_table2_NFS_harmonised.py  --  Reply to Cirillo (PNAS Letter, 2026): Table 2

Change in household purchases, mean of the four quarters after (1953Q4-1954Q3)
versus the four quarters before (1952Q4-1953Q3) the 26 September 1953 end of
sugar rationing, from the National Food Survey series redistributed in
data/NFS/ (digitised by Gracner et al. 2024).

Codebook note. In Dataset_1000days_NUTRITION_Fig_S3.dta the variables whose
names end in "_gday10" are stored in TENS of grams per day (multiply by 10 to
obtain g/day); the others are in grams per day. The variable label of
bread_gday10 omits the "(in 10s)" flag that the other five carry.

"As reported" reproduces the comment's Table 1 by adding the _gday10 series
at face value to the gram series; "harmonised" rescales every series to
g/day first. The component series of each category are inferred from
reproducing the comment's figures (the comment does not list them).

Input:  data/NFS/Dataset_1000days_NUTRITION_Fig_S3.dta,
        data/NFS/Dataset_1000days_Fig_S1.dta (sugar, kcal/day)
Output: output/reply_to_cirillo/table2_NFS_harmonised.csv
Requires: pandas.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
NFS = os.path.join(ROOT, "data", "NFS")
OUTDIR = os.path.join(ROOT, "output", "reply_to_cirillo")
os.makedirs(OUTDIR, exist_ok=True)

PRE = (pd.Period("1952Q4"), pd.Period("1953Q3"))
POST = (pd.Period("1953Q4"), pd.Period("1954Q3"))


def window_means(df):
    q = pd.PeriodIndex(df["qdate"], freq="Q")
    pre = df[(q >= PRE[0]) & (q <= PRE[1])].mean(numeric_only=True)
    post = df[(q >= POST[0]) & (q <= POST[1])].mean(numeric_only=True)
    return pre, post


def pct(pre, post, expr):
    a = eval(expr, {}, pre.to_dict())
    b = eval(expr, {}, post.to_dict())
    return 100 * (b / a - 1)


nut = pd.read_stata(os.path.join(NFS, "Dataset_1000days_NUTRITION_Fig_S3.dta"))
pre, post = window_means(nut)

# category: (as-reported expression, harmonised expression)
CATS = {
    "Fats (butter, margarine, lard)":
        ("butter_gday + margarine_gday + lard_gday",
         "butter_gday + margarine_gday + lard_gday"),
    "Fruit and vegetables (canned fruit, fresh fruit, vegetables)":
        ("otherfruit_gday + freshfruit_gday10 + Vegetables_gday10",
         "otherfruit_gday + 10*freshfruit_gday10 + 10*Vegetables_gday10"),
    "Dairy (milk, welfare/school milk, cheese)":
        ("MilkCream_gday10 + lmilkwelfareandschoolpt_gday + Cheese_gday",
         "10*MilkCream_gday10 + lmilkwelfareandschoolpt_gday + Cheese_gday"),
    "Cereals (as reported: total cereals + bread + flour; harmonised: total cereals)":
        ("totalcereals_gday10 + bread_gday10 + flour_gday",
         "10*totalcereals_gday10"),
    "Meat (meat, fish, bacon and ham)":
        ("Meats_gday10 + Fish_gday + baconandhamuncooked_gday",
         "10*Meats_gday10 + Fish_gday + baconandhamuncooked_gday"),
}

rows = []
for cat, (rep, harm) in CATS.items():
    rows.append(dict(category=cat,
                     as_reported_pct=round(pct(pre, post, rep), 2),
                     harmonised_pct=round(pct(pre, post, harm), 2)))

# sugar (kcal/day) from the calorie series, as in the comment's first row
sug = pd.read_stata(os.path.join(NFS, "Dataset_1000days_Fig_S1.dta"))
spre, spost = window_means(sug)
rows.insert(0, dict(category="Sugar (kcal/day)",
                    as_reported_pct=round(pct(spre, spost, "Sugars_gdaycal"), 2),
                    harmonised_pct=round(pct(spre, spost, "Sugars_gdaycal"), 2)))

# single-series changes (all in g/day) for reference
singles = []
for c in nut.columns:
    if c == "qdate":
        continue
    singles.append(dict(category=f"series: {c}" + (" (x10 = g/day)" if c.endswith("gday10") else ""),
                        as_reported_pct=round(pct(pre, post, c), 2),
                        harmonised_pct=round(pct(pre, post, c), 2)))

out = pd.DataFrame(rows + singles)
out.to_csv(os.path.join(OUTDIR, "table2_NFS_harmonised.csv"), index=False)
print(out.to_string(index=False))
