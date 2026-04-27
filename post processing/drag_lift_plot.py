import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Parameter
t_transient = 0.0

# Read
df = pd.read_csv("simulation_data//forces.csv", skipinitialspace=True)
df.columns = df.columns.str.strip()

t   = df["t_phys"].values
Cd  = df["Cd"].values
Cl  = df["Cl"].values


# Mean
mask = t >= t_transient
Cd_mean = np.mean(Cd[mask])
Cl_rms = np.sqrt(np.mean(Cl[mask]**2))

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

ax1.plot(t, Cd, "b-", linewidth = 0.8)
ax1.axhline(Cd_mean, color="b", linewidth=1.5, linestyle="--", label=f"Cd_mean = {Cd_mean:.4f}")
ax1.set_ylabel("Cd [-]")
ax1.legend(fontsize=9)
ax1.grid(True, linestyle="--", alpha=0.4)

ax2.plot(t, Cl, "r-", linewidth=0.8)
ax2.axhline(0, color="gray", linewidth=0.5, linestyle="--")
ax2.set_ylabel("Cl [-]")
ax2.set_xlabel("t [s]")
ax2.legend([f"Cl (rm = {Cl_rms:.4f})"], fontsize=9)
ax2.grid(True, linestyle="--", alpha=0.4)

fig.suptitle(f"Drag & Lift | Cd_mean = {Cd_mean:.4f} | Cl_rms = {Cl_rms:.4f}", fontsize=11)
plt.tight_layout()
plt.show()