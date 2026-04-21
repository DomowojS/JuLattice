import pandas as pd
import matplotlib.pyplot as plt
import math

# ---- Parameter (have to be same as JuLattice!!!) ----
Mach_Number = 0.05
# U_lat = Ma * cs
lattice_U_inf = Mach_Number / math.sqrt(3)
Radius      = 0.0115
D           = 2 * Radius
length_Y    = 0.6
nu          = 1e-6
Re          = 2760
U_inf       = Re * nu / D

# ---- read .CSV ----
df = pd.read_csv("wake_profil.csv", skipinitialspace=True)
# delete "space"
df.columns = df.columns.str.strip()

# ---- take last timestep (= cumulativ mean) ----
t_last = df["t_phys"].max()
last = df[df["t_phys"] == t_last].copy()

# ---- compute non-dimensional coordinates ----
last["y_D"] = (last["y_phys"] - length_Y / 2) / D
last["U_norm"] = last["mean_u"] / lattice_U_inf

# ---- Plot ----
fig, ax = plt.subplots(figsize=(5, 8))

ax.plot(last["U_norm"], last["y_D"], "k-o", markersize=3, linewidth=1.5)
ax.set_xlabel(r"$U_m / U_\infty$ [-]")
ax.set_ylabel(r"$y / D$ [-]")
ax.set_title(f"Wake velocity profile\n(x = cylinder_x + 3D,  t = {t_last:.2f} s)")
ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")   # Zylindermitte
ax.axhline( 0.5, color="lightgray", linewidth=0.5, linestyle=":")  # Zylinderrand
ax.axhline(-0.5, color="lightgray", linewidth=0.5, linestyle=":")
ax.grid(True, linestyle="--", alpha=0.4)
#ax.set_xlim(0.5, 1.05)

plt.tight_layout()
#plt.savefig("wake_profil.png", dpi=150)
plt.show()