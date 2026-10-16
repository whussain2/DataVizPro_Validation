"""Independent Python references for Table 4.

1. fastmixture 1.3 (ADMIXTURE likelihood model) on the 10,000-marker matrix, K = 5
2. statsmodels 0.15 OLS and logistic regression on the rows DataVizPro used

Run from the repository root after scripts 02 and 10 have written their inputs:
    pip install fastmixture==1.3.0 bed-reader statsmodels==0.15.0
    python python/references.py
"""
import os, subprocess, sys
import numpy as np
import pandas as pd

TMP = os.path.join("results", "intermediate")


def run_fastmixture():
    src = os.path.join(TMP, "prep10.tsv")
    if not os.path.exists(src):
        sys.exit("Run scripts/02_rice_qc.R first.")
    from bed_reader import to_bed
    d = pd.read_csv(src, sep="\t", index_col=0)
    X = d.values.astype("float32")
    X[X == -9] = np.nan
    bed = os.path.join(TMP, "rice10.bed")
    to_bed(bed, X, properties={"iid": list(d.index), "fid": list(d.index), "sid": list(d.columns)})
    for seed in range(1, 6):   # five seeds; all reach the same likelihood
        subprocess.run(["fastmixture", "--bfile", bed[:-4], "--K", "5", "--seed", str(seed),
                        "--out", os.path.join(TMP, "fm_k5")], check=True)


def run_statsmodels():
    import statsmodels.api as sm
    out = {}
    for name, fn in (("ols", sm.OLS), ("logit", sm.Logit)):
        path = os.path.join(TMP, f"{name}.csv")
        if not os.path.exists(path):
            print(f"skip {name}: run scripts/10_phenotypes.R first"); continue
        df = pd.read_csv(path)
        model = fn(df.iloc[:, 0], sm.add_constant(df.iloc[:, 1:]))
        res = model.fit() if name == "ols" else model.fit(disp=0, tol=1e-12, maxiter=200)
        out[name] = res.params.values
    pd.DataFrame(dict((k, pd.Series(v)) for k, v in out.items())).to_csv(
        os.path.join(TMP, "statsmodels_coefs.csv"), index=False)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "fastmixture"):
        run_fastmixture()
    if what in ("all", "statsmodels"):
        run_statsmodels()
