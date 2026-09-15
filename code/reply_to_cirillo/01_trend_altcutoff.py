"""
01_trend_altcutoff.py  --  Reply to Cirillo (PNAS Letter, 2026): Table 1

Cox proportional hazards (Breslow ties) with cluster-robust (Lin-Wei)
standard errors clustered by month of birth, mirroring
    stcox ..., cluster(yearmobirth)
in Stata (the same specification as SI Table S2 of the paper).

For each of the five cancers and three cutoff dates (published 26 Sep 1953;
end June 1953; end March 1953) the script reports
  (a) the categorical event-study cells of SI Table S2 (replication check),
  (b) HR per six months of the first 1,000 days spent under rationing
      (continuous exposure, pre-period cells kept as dummies),
  (c) HR per ordered exposure category (score trend),
  (d) a 3-df Wald test of departure from a linear gradient across the five
      rationed cells (equal successive increments in log HR).

Input:  251025_sugar_data_for_cox.dta  (UK Biobank derivative, application
        89068; NOT redistributed -- see README section 3). Set the folder
        holding it in the environment variable UKB or pass it as the first
        command-line argument.
Output: output/reply_to_cirillo/trend_altcutoff_summary.csv  (one row per
        cancer x cutoff) and altcutoff_cells.csv (cell-level HRs).

Requires: numpy, pandas, scipy. Runtime ~10-20 minutes.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

# ---------- paths -----------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
UKB = os.environ.get("UKB") or (sys.argv[1] if len(sys.argv) > 1 else None)
if UKB is None:
    sys.exit("Set the UKB environment variable (folder holding "
             "251025_sugar_data_for_cox.dta) or pass it as the first argument.")
OUTDIR = os.path.join(ROOT, "output", "reply_to_cirillo")
os.makedirs(OUTDIR, exist_ok=True)

COVS = ["smoking", "college"] + [f"pca{i}" for i in range(1, 11)] + \
       ["home_longi", "home_lati", "Townsend_deprivation_index"]
OUTS = {"lung": ("diag_LungCancer_time", "diag_LungCancer_event", None),
        "liver": ("diag_C22_time", "diag_C22_event", None),
        "rectum": ("diag_C20_Rectum_time", "diag_C20_Rectum_event", None),
        "prostate": ("diag_C61_Prostate_time", "diag_C61_Prostate_event", 1),
        "breast": ("diag_C50_Breast_time", "diag_C50_Breast_event", 0)}
NEEDED = ["yearqbirth", "yearmobirth", "study", "month_birth", "male"] + COVS + \
         [v for tv, ev, _ in OUTS.values() for v in (tv, ev)]

d = pd.read_stata(os.path.join(UKB, "251025_sugar_data_for_cox.dta"), columns=NEEDED)

# ---------- cohort indices --------------------------------------------------
# rv = quarter of birth + 41 (8 = 1951q4), as in the package do-files.
yq = pd.DatetimeIndex(d["yearqbirth"])
d["rv"] = (yq.year - 1960) * 4 + (yq.quarter - 1) + 41
ym = pd.DatetimeIndex(d["yearmobirth"])
d["m"] = (ym.year - 1960) * 12 + (ym.month - 1)                  # Stata %tm
C0 = (1953 - 1960) * 12 + 8                                       # tm(1953m9)


def make_study(rv):
    """Exposure-cohort cells of SI Table S2 (4 = reference, born Jul-Dec 1954)."""
    s = pd.Series(np.nan, index=rv.index)
    s[rv > 24] = 1
    s[rv.isin([23, 24])] = 2
    s[rv.isin([21, 22])] = 3
    s[rv.isin([19, 20])] = 4
    s[rv.isin([16, 17, 18])] = 5
    s[rv.isin([14, 15])] = 6
    s[rv.isin([12, 13])] = 7
    s[rv.isin([10, 11])] = 8
    s[rv.isin([8, 9])] = 9
    return s


assert (make_study(d["rv"]) == d["study"]).all(), "study reconstruction mismatch"

# k = 0: published cutoff (Sep 1953); k = 1: one quarter earlier (end Jun 1953);
# k = 2: two quarters earlier (end Mar 1953). Cells shift with the cutoff.
for k in (0, 1, 2):
    d[f"study{k}"] = make_study(d["rv"] + k)
    Ck = C0 - 3 * k
    e = (Ck - d["m"] + 9).clip(0, 33)      # months of the first 1,000 days under rationing
    d[f"exp6_{k}"] = e / 6.0
    d[f"score{k}"] = (d[f"study{k}"] - 4).clip(lower=0)


# ---------- Cox with Breslow ties -----------------------------------------
def _risk_sums(X, w, inv, J):
    """S0, S1, S2 for risk sets {t >= ut[j]} via per-time group sums and
    reverse cumulative sums (J unique times, inv = time index per row)."""
    p = X.shape[1]
    g0 = np.zeros(J); g1 = np.zeros((J, p)); g2 = np.zeros((J, p, p))
    np.add.at(g0, inv, w)
    np.add.at(g1, inv, w[:, None] * X)
    for j in range(J):
        m = inv == j
        Xj = X[m]
        g2[j] = Xj.T @ (w[m][:, None] * Xj)
    S0 = np.cumsum(g0[::-1])[::-1]
    S1 = np.cumsum(g1[::-1], axis=0)[::-1]
    S2 = np.cumsum(g2[::-1], axis=0)[::-1]
    return S0, S1, S2


def cox_breslow(X, t, ev, cluster, tol=1e-9, maxit=50):
    n, p = X.shape
    ut, inv = np.unique(t, return_inverse=True)         # ascending unique times
    J = len(ut)
    dj = np.zeros(J); np.add.at(dj, inv, ev)
    ev_idx = ev == 1
    xsum_ev = X[ev_idx].sum(0)
    beta = np.zeros(p)
    for it in range(maxit):
        eta = X @ beta
        w = np.exp(eta)
        S0, S1, S2 = _risk_sums(X, w, inv, J)
        xbar = S1 / S0[:, None]
        g = xsum_ev - (dj[:, None] * xbar).sum(0)
        H = np.zeros((p, p))
        for j in range(J):
            if dj[j] == 0:
                continue
            H += dj[j] * (S2[j] / S0[j] - np.outer(xbar[j], xbar[j]))
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    eta = X @ beta; w = np.exp(eta)
    S0, S1, S2 = _risk_sums(X, w, inv, J)
    xbar = S1 / S0[:, None]
    H = np.zeros((p, p))
    for j in range(J):
        if dj[j] == 0:
            continue
        H += dj[j] * (S2[j] / S0[j] - np.outer(xbar[j], xbar[j]))
    Hinv = np.linalg.inv(H)
    # score residuals (Lin-Wei): U_i = ev_i (x_i - xbar(T_i)) - w_i [x_i A(T_i) - B(T_i)]
    a = dj / S0
    A = np.cumsum(a)
    B = np.cumsum(a[:, None] * xbar, axis=0)
    Ai = A[inv]; Bi = B[inv]; xbari = xbar[inv]
    U = ev[:, None] * (X - xbari) - w[:, None] * (X * Ai[:, None] - Bi)
    cl_ids, cl_inv = np.unique(cluster, return_inverse=True)
    G = len(cl_ids)
    Ug = np.zeros((G, p))
    np.add.at(Ug, cl_inv, U)
    meat = Ug.T @ Ug
    V = Hinv @ meat @ Hinv * (G / (G - 1))              # Stata small-sample factor
    return beta, V, int(ev.sum()), n, G


def wald(beta, V, R):
    R = np.atleast_2d(R)
    diff = R @ beta
    chi2 = float(diff @ np.linalg.solve(R @ V @ R.T, diff))
    return chi2, R.shape[0], 1 - stats.chi2.cdf(chi2, R.shape[0])


def build(df, cohort_terms, sexvar):
    """Design matrix: cohort terms, male dummy (if both sexes), month-of-birth
    dummies, covariates. Returns X and column names."""
    cols, names = [], []
    for nm, s in cohort_terms:
        cols.append(s.values.astype(float)); names.append(nm)
    if sexvar:
        cols.append(df["male"].values.astype(float)); names.append("male")
    for mth in range(2, 13):
        cols.append((df["month_birth"].values == mth).astype(float)); names.append(f"mob{mth}")
    for c in COVS:
        cols.append(df[c].values.astype(float)); names.append(c)
    return np.column_stack(cols), names


rows, cells = [], []
for o, (tv, evv, sex) in OUTS.items():
    base = d.dropna(subset=COVS + ["home_longi"])
    if sex is not None:
        base = base[base["male"] == sex]
    t = base[tv].values.astype(float); ev = base[evv].values.astype(float)
    cl = base["m"].values
    for k in (0, 1, 2):
        sv = base[f"study{k}"]
        ok = sv.notna().values
        df = base[ok]; tt = t[ok]; ee = ev[ok]; cc = cl[ok]
        sv = sv[ok]
        # (a) categorical event study
        present = [j for j in (1, 2, 3, 5, 6, 7, 8, 9) if (sv == j).sum() > 0]
        rat = [j for j in (5, 6, 7, 8, 9) if j in present]
        terms = [(f"c{j}", (sv == j).astype(float)) for j in present]
        X, names = build(df, terms, sex is None)
        b, V, nev, n, G = cox_breslow(X, tt, ee, cc)
        idx = {nm: i for i, nm in enumerate(names)}
        se = np.sqrt(np.diag(V))
        hr = {j: (np.exp(b[idx[f"c{j}"]]),
                  np.exp(b[idx[f"c{j}"]] - 1.96 * se[idx[f"c{j}"]]),
                  np.exp(b[idx[f"c{j}"]] + 1.96 * se[idx[f"c{j}"]])) for j in present}
        Rj = np.zeros((len(rat), len(b)))
        for i, j in enumerate(rat):
            Rj[i, idx[f"c{j}"]] = 1
        chi_j, _, p_j = wald(b, V, Rj)
        Rp = np.zeros((3, len(b)))
        for i, j in enumerate((1, 2, 3)):
            Rp[i, idx[f"c{j}"]] = 1
        chi_p, _, p_pre = wald(b, V, Rp)
        # (d) departure from linearity: (c6-c5)-(c7-c6)=0, etc.
        trip = [(rat[i], rat[i + 1], rat[i + 2]) for i in range(len(rat) - 2)]
        Rl = np.zeros((len(trip), len(b)))
        for i, (a1, a2, a3) in enumerate(trip):
            Rl[i, idx[f"c{a1}"]] = 1; Rl[i, idx[f"c{a2}"]] = -2; Rl[i, idx[f"c{a3}"]] = 1
        chi_l, _, p_lin = wald(b, V, Rl)
        # (b) continuous trend (6-month units) with pre-cells as dummies
        terms = [("exp6", df[f"exp6_{k}"])] + [(f"c{j}", (sv == j).astype(float)) for j in (1, 2, 3)]
        X, names = build(df, terms, sex is None)
        b2, V2, *_ = cox_breslow(X, tt, ee, cc)
        hr6 = np.exp(b2[0]); se6 = np.sqrt(V2[0, 0]); p6 = 2 * stats.norm.sf(abs(b2[0] / se6))
        # (c) ordered score trend
        terms = [("score", df[f"score{k}"])] + [(f"c{j}", (sv == j).astype(float)) for j in (1, 2, 3)]
        X, names = build(df, terms, sex is None)
        b3, V3, *_ = cox_breslow(X, tt, ee, cc)
        hrs = np.exp(b3[0]); ses = np.sqrt(V3[0, 0]); ps = 2 * stats.norm.sf(abs(b3[0] / ses))
        h9 = hr.get(9, (np.nan,) * 3)
        rows.append(dict(outcome=o, k=k, n=n, events=nev, clusters=G,
                         HR24=h9[0], HR24_lo=h9[1], HR24_hi=h9[2],
                         HRmax=hr[rat[-1]][0], maxcell=rat[-1],
                         chi2_joint=chi_j, p_joint=p_j, p_pretrend=p_pre,
                         chi2_nonlin=chi_l, p_nonlin=p_lin,
                         HR_per6m=hr6, HR_per6m_lo=np.exp(b2[0] - 1.96 * se6),
                         HR_per6m_hi=np.exp(b2[0] + 1.96 * se6), p_trend6m=p6,
                         HR_per_cell=hrs, HR_per_cell_lo=np.exp(b3[0] - 1.96 * ses),
                         HR_per_cell_hi=np.exp(b3[0] + 1.96 * ses), p_trend_cell=ps))
        for j in present:
            cells.append(dict(outcome=o, k=k, cell=j, HR=hr[j][0], lo=hr[j][1], hi=hr[j][2]))
        print(f"{o:9s} k={k} n={n} ev={nev} G={G} | HR+24={h9[0]:.3f} [{h9[1]:.2f},{h9[2]:.2f}] "
              f"joint chi2={chi_j:.2f} p={p_j:.2e} | pre p={p_pre:.3f} | "
              f"trend/6m HR={hr6:.3f} p={p6:.2e} | per-cell HR={hrs:.3f} p={ps:.2e} | "
              f"nonlin chi2({len(trip)})={chi_l:.2f} p={p_lin:.3f}")
        sys.stdout.flush()

pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "trend_altcutoff_summary.csv"), index=False)
pd.DataFrame(cells).to_csv(os.path.join(OUTDIR, "altcutoff_cells.csv"), index=False)
print("\nsaved to", OUTDIR)
