"""Figures of the ALF benchmark from results/*.csv (run after the C++ benchmarks). Output: img/*.png."""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = pathlib.Path(__file__).parent
res, img = here / "results", here / "img"
img.mkdir(exist_ok=True)


def load(p):
    return np.genfromtxt(p, delimiter=",", names=True)


# 1. Sum-rule identity: worst relative error over all orders d, per latitude.
fig, ax = plt.subplots(figsize=(8, 4))
for alg, n in (("sfc", 1800), ("msfc", 2700), ("era", 2700)):
    a = load(res / f"alf_identities_{alg}_n{n}.csv")
    lat = np.unique(a["lat_deg"])
    worst = np.array([a["max_rel_error"][a["lat_deg"] == l].max() for l in lat])
    ax.semilogy(lat, np.maximum(worst, 1e-18), ".-", ms=3, label=f"{alg.upper()} (n={n})")
ax.set_xlabel("latitude [deg]")
ax.set_ylabel("max relative error of the sum rule")
ax.set_title("ALF sum-rule identity: worst case over orders d")
ax.grid(alpha=.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(img / "alf_sum_rule_vs_latitude.png", dpi=150)

# 2. Production kernel against the sequential oracle: |V_kernel - V_oracle| along the meridian.
k = load(res / "kernel_vs_oracle_V_egm2008_n2190.csv")
fig, ax = plt.subplots(figsize=(8, 4))
ax.semilogy(k["lat_deg"], np.maximum(k["abs_err"], 1e-12), ".", ms=3)
ax.axhline(np.spacing(6.3e7), color="r", ls="--", lw=1, label="1 ULP of V (6.3e7 m²/s²)")
ax.set_xlabel("latitude [deg]")
ax.set_ylabel("|V_kernel − V_oracle| [m²/s²]")
ax.set_title("fast kernel vs oracle, EGM2008 n=2190")
ax.grid(alpha=.3, which="both")
ax.legend()
fig.tight_layout()
fig.savefig(img / "kernel_vs_oracle_V.png", dpi=150)
print("max |dV| =", k["abs_err"].max(), "ULP =", np.spacing(6.3e7))
