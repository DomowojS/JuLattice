# JuLattice — MRT D2Q9 with Grid Refinement

![Demo](media/Vorticity_Cylinder.gif)

**JuLattice** is a 2D Lattice Boltzmann Method (LBM) solver written in Julia. This branch implements a **two-level grid refinement** scheme on top of the D2Q9 MRT collision operator, enabling higher spatial resolution in a sub-domain of interest without refining the entire mesh.

## Features

- D2Q9 MRT collision with configurable BGK and acoustic relaxation frequencies
- Two-level grid refinement (coarse + fine) with acoustic scaling (Δx/2, Δt/2)
- Bouzidi curved boundary conditions for immersed objects
- Rectangular obstacle with arbitrary angle
- Lift and drag force output (cL, cD) at every time step
- Contour plots of vorticity, u-velocity, v-velocity, and velocity magnitude via GLMakie
- Post-processing script (`output/Analysis.jl`) with Welch FFT, Strouhal number, convergence diagnostics, and TikZ data export

## Requirements

![Julia version](https://img.shields.io/badge/julia-1.9%2B-blue)

**Simulation dependencies** (defined in `Project.toml`):
- `GLMakie`
- `MeshGrid`
- `Revise`

**Analysis script dependencies** (not in `Project.toml` — install separately if needed):
- `FFTW`
- `Plots`

## Getting Started

### 1. Clone the repository and check out this branch

```bash
git clone https://github.com/DomowojS/JuLattice.git
cd JuLattice
git checkout MRT_D2Q9_gridRefinement
```

### 2. Activate the project environment

Start Julia in the project directory and activate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### 3. Configure the simulation

All simulation parameters are set directly at the top of `JuLattice_main.jl` inside the `run()` function. Open the file and edit the **User Settings** block:

```julia
# Domain
lengthX = 8.0          # m  (channel length)
lengthY = 3.0          # m  (channel height)

# Fine grid sub-domain
lengthXFine       = 4.0    # m  (width of fine region)
lengthYFine       = 2.4    # m  (height of fine region)
positionFineGridX = 1.5    # m  (lower-left x anchor)
positionFineGridY = 0.3    # m  (lower-left y anchor)

# Obstacle (rectangle)
d        = 0.5             # m  (characteristic length, used for Re)
angleDeg = 30.0            # degrees
positionX = 3.0            # m
positionY = lengthY / 2    # m

# Fluid properties
reynoldsNumber = 300
machNumber     = 0.1       # Ma = U / c_s
viscosity      = 1e-4      # m²/s

# Discretisation
deltaX         = 0.05      # m per coarse lattice unit
simulationTime = 3600.0    # s
```

The fine grid spacing and time step are derived automatically:
```
deltaXFine = deltaX / 2
deltaTFine = deltaT / 2     (acoustic scaling)
```

The fine box is snapped to the nearest integer number of coarse cells so that the coarse-fine coupling nodes align correctly.

### 4. Run the simulation

```julia
include("JuLattice_main.jl")
JuLattice.run()
```

This will:
1. Print a log header with discretisation settings
2. Initialise both grids to equilibrium at the inflow velocity
3. Advance the simulation (2 fine sub-steps per coarse step)
4. Write force data to `output/forces.txt` (columns: `t  cL  cD`)
5. Save contour snapshots to `output/` at regular intervals

### 5. Blockage ratio study

The cases studied in this branch are defined in `output/CaseTable.txt`. Each case changes `lengthY` (and `positionFineGridY`) to vary the blockage ratio while keeping everything else fixed:

| Case | B (%) | lengthY | positionFineGridY |
|------|--------|---------|-------------------|
| 1    | 4.9    | 24.0    | 10.8              |
| 2    | 9.9    | 12.0    | 4.8               |
| 3    | 14.8   | 8.0     | 2.8               |
| 4    | 19.7   | 6.0     | 1.8               |
| 5    | 29.6   | 4.0     | 0.8               |
| 6    | 45.6   | 2.6     | 0.1               |

Fixed parameters across all cases: `Re = 300`, `Ma = 0.1`, `ν = 1e-4 m²/s`, `d = 0.5 m`, `angleDeg = 30°`, `deltaX = 0.02 m`.

Run each case by adjusting `lengthY` and `positionFineGridY` in `run()` and saving the output force file as `output/B_XX_X_forces.txt` (e.g. `B_19_7_forces.txt`).

### 6. Post-process results

With one or more `B_*_forces.txt` files in `output/`, run the analysis script from the Julia REPL:

```julia
include("output/Analysis.jl")
```

This produces:

| Output directory | Contents |
|---|---|
| `output/cL/` | Raw cL time-series plots per case |
| `output/cD/` | Raw cD time-series plots per case |
| `output/mean/` | Cumulative mean (from t = 0) per case |
| `output/std/` | Cumulative std (from t = 0) per case |
| `output/FFT/` | Welch amplitude spectrum per case |
| `output/Blockage/` | Mean, std, frequency, Strouhal vs. blockage ratio |
| `output/tikz/` | Two-column `.txt` files for `\addplot table` in pgfplots |

The analysis constants at the top of `Analysis.jl` can be adjusted:

```julia
const T_START     = 4000.0   # discard transient before this time [s]
const T_WINDOW    = 4000.0   # window for mean/std blockage plots [s]
const CHAR_LENGTH = 0.5      # d [m]  — for Strouhal St = f·d/U
const INFLOW_VEL  = 0.06     # U [m/s]
```

The Strouhal number is computed as:

```
St = f_dominant × d / U
```

where `f_dominant` is the peak frequency from the Welch spectrum over the full `[T_START, end]` window.

## Project Structure

```
JuLattice/
├── JuLattice_main.jl     # Entry point — user settings and main loop
├── src/
│   ├── Kernel.jl         # MRT collision + streaming (shared by coarse and fine)
│   ├── GridSetup.jl      # Domain and node classification (coarse)
│   ├── GridRefinement.jl # Fine grid classification and C↔F synchronisation
│   ├── Plotter.jl        # GLMakie live/snapshot plots
│   └── IO.jl             # Force output and logging
├── output/
│   ├── Analysis.jl       # Post-processing: Welch FFT, Strouhal, TikZ export
│   ├── CaseTable.txt     # Blockage ratio case definitions
│   └── B_*_forces.txt    # Force output files (one per case)
├── Project.toml
└── README.md
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Other Branches

### [`main`](https://github.com/DomowojS/JuLattice/tree/main)
Single-grid D2Q9 BGK solver with JSON configuration files. Use this for simpler setups.

### [`JuLattice_for_teaching`](https://github.com/DomowojS/JuLattice/tree/JuLattice_for_teaching)
Hard-coded, minimal version intended for educational use.
