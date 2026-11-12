"""Figure of the meridian benchmark: |gsh - GrafLab| / max|GrafLab| per functional vs latitude."""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = pathlib.Path(__file__).parent
r = here / "results"
g = np.genfromtxt(r / "gsh_meridian.csv", delimiter=",", names=True)
f = np.genfromtxt(r / "graflab_era_meridian.csv", delimiter=",", names=True)
mgal, E = 1e-5, 1e-9
pairs = [("T", g["T_m2s2"], f["T_m2s2"]), ("zeta_ell", g["zeta_ell_m"], f["zeta_ell_m"]),
         ("Dg_sa", g["Dg_sa_ms2"], f["Dg_sa_mgal"] * mgal), ("dg", g["dg_ms2"], f["dg_mgal"] * mgal),
         ("V_xx", g["V_xx_s2"], f["Vxx_E"] * E), ("V_yy", g["V_yy_s2"], f["Vyy_as_exported_E"] * E),
         ("V_zz", g["V_zz_s2"], f["Vzz_E"] * E), ("V_xz", g["V_xz_s2"], f["Vxz_E"] * E)]
fig, ax = plt.subplots(figsize=(9, 5))
for name, a, b in pairs:
    d = np.abs(a - b) / np.nanmax(np.abs(b))
    ax.semilogy(g["lat_deg"], np.maximum(d, 1e-17), ".-", ms=3, lw=.7, label=name)
ax.axhline(1e-12, color="k", ls="--", lw=1, label="criterion 1e-12")
ax.set_xlabel("geodetic latitude [deg]")
ax.set_ylabel("|gsh − GrafLab| / max|GrafLab|")
ax.set_title("Meridian, EGM2008 n=2190, 181 points")
ax.grid(alpha=.3, which="both")
ax.legend(ncol=3, fontsize=8)
fig.tight_layout()
fig.savefig(here / "img" / "meridian_gsh_vs_graflab.png", dpi=150)
