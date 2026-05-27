import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog, simpledialog

# Parameter
tk.Tk().withdraw()

t_transient = 40.0
Re = simpledialog.askfloat("Reynolds Number", "Enter Re:", initialvalue=2760.0)
if Re is None:
    raise SystemExit("no Re entered!!!")

# Read
file_path = filedialog.askopenfilename(
    title="Select forces csv",
    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
)
if not file_path:
    raise SystemExit("No file selected.")

df = pd.read_csv(file_path, skipinitialspace=True)
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

# FFT of Cl (post-transient only)
Cl_fft = Cl[mask]
t_fft  = t[mask]
dt     = np.mean(np.diff(t_fft))        # calculate timestep size
N      = len(Cl_fft)                    # number of samples

freqs  = np.fft.rfftfreq(N, d=dt)       # frequency axis
power  = np.abs(np.fft.rfft(Cl_fft))    # amplitude at each frequency

# Calculate Strouhal number
D = 0.023                               # cylinder diameter
U = Re * 1e-6 / D                       # inflow vel
St_axis = freqs * D / U                 

peak_idx = np.argmax(power[1:]) + 1     # index of largest amplitude
f_peak   = freqs[peak_idx]              # f at largest amplitude
St_peak  = St_axis[peak_idx]            # Strouhal nubmer at largest amplitude

fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(15, 5))

# Strouhal axis
ax3.plot(St_axis, power, "k-", linewidth=0.8,
         label=f"U = {U:.4f} m/s | D = {D:.3f} m")
ax3.axvline(St_peak, color="r", linestyle="--", 
            label=f"St = {St_peak:.4f}  (f = {f_peak:.4f} Hz)")
ax3.set_xlabel("St = f·D/U [-]")
ax3.set_ylabel("spectral amplitude of Cl[-]")
ax3.set_title("FFT of Cl (Strouhal)")
ax3.set_xlim(left=0)
#ax3.set_ylim(0, 100)
ax3.legend(fontsize=9)
ax3.grid(True, linestyle="--", alpha=0.4)

# frequency axis
ax4.plot(freqs, power, "k-", linewidth=0.8,
         label=f"U = {U:.4f} m/s | D = {D:.3f} m")
ax4.axvline(f_peak, color="r", linestyle="--",
            label=f"f = {f_peak:.4f} Hz (St = {St_peak:.4f})")
ax4.set_xlabel("f [Hz]")
ax4.set_ylabel("spectral amplitude of Cl [-]")
ax4.set_title("FFT of Cl (Frequency)")
ax4.set_xlim(left=0)
ax4.legend(fontsize=9)
ax4.grid(True, linestyle="--", alpha=0.4)

fig2.suptitle(f"FFT of Cl  |  Re = {Re:.0f}  |  Peak: St = {St_peak:.4f}, f = {f_peak:.4f} Hz", fontsize=11)
plt.tight_layout()
plt.show()