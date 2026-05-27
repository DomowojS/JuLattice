import pandas as pd
import matplotlib.pyplot as plt
import math
import tkinter as tk
import os
from tkinter import filedialog, simpledialog

# ---- Parameter (have to be same as JuLattice!!!) ----
Radius      = 0.0115
D           = 2 * Radius
# length_Y    = 0.6 <-- large Domain
# length_Y    = 10 * D
length_Y    = 15 * D


# ---- read .CSV ----
tk.Tk().withdraw()
Re = simpledialog.askfloat(
    "Reynolds Number",
    "Enter Re: ",
    initialvalue=2760.0
)
if Re is None:
    raise SystemExit("No Re number entered! :(")

Ma = simpledialog.askfloat(
    "Mach Number",
    "Enter Ma: ",
    initialvalue=0.1
)
if Ma is None:
    raise SystemExit("No Ma number entered! :(")

lattice_U_inf = Ma / math.sqrt(3)


file_path = filedialog.askopenfilename(
    title="Select wake profile CSV",
    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
)
if not file_path:
    raise SystemExit("No file selected.")

fname = os.path.basename(file_path)
distance = "6D" if "_6D_" in fname else "3D"

df = pd.read_csv(file_path, skipinitialspace=True)
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
ax.set_title(f"Re = {Re:.0f} | Wake velocity profile\n(x = cylinder_x + {distance},  t = {t_last:.2f} s)")
ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")   # Zylindermitte
ax.axhline( 0.5, color="lightgray", linewidth=0.5, linestyle=":")  # Zylinderrand
ax.axhline(-0.5, color="lightgray", linewidth=0.5, linestyle=":")
ax.grid(True, linestyle="--", alpha=0.4)
ax.set_xlim(left=0.2)

plt.tight_layout()
#plt.savefig("wake_profil.png", dpi=150)
plt.show()