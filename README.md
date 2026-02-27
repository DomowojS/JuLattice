# JuLattice

![Demo](media/Vorticity_Cylinder.gif)

**JuLattice** is a 2D Lattice Boltzmann Method (LBM) solver written in Julia. It implements the D2Q9 lattice scheme with the BGK collision operator (single relaxation rate).

## Features

- D2Q9 BGK collision with user-defined relaxation rate
- Rectangular obstacle with arbitrary angle, Bouzidi curved boundary conditions
- Lift and drag force output (cL, cD) at every time step
- Real-time contour plots of vorticity and velocity via GLMakie

## Requirements

![Julia version](https://img.shields.io/badge/julia-1.9%2B-blue)

**Dependencies** (defined in `Project.toml`):
- `GLMakie`
- `MeshGrid`
- `Revise`

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/DomowojS/JuLattice.git
cd JuLattice
```

### 2. Activate the project environment

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### 3. Configure the simulation

All parameters are set at the top of `JuLattice_main.jl` inside the `run()` function:

```julia
# Domain
lengthX = 8.0          # m
lengthY = 3.0          # m

# Obstacle (rectangle)
d        = 0.5         # m  (characteristic length, used for Re)
angleDeg = 30.0        # degrees
positionX = 3.0        # m
positionY = lengthY/2  # m

# Fluid properties
reynoldsNumber = 300
machNumber     = 0.1   # Ma = U / c_s  (keep < 0.1 for incompressible)
viscosity      = 1e-4  # m²/s

# Discretisation
deltaX         = 0.05  # m per lattice unit
simulationTime = 3600.0  # s
```

### 4. Run the simulation

```julia
include("JuLattice_main.jl")
JuLattice.run()
```

This will:
1. Print a log header with discretisation settings
2. Initialise the grid to equilibrium at the inflow velocity
3. Advance the simulation, updating the real-time GLMakie plots
4. Write force data to `output/forces.txt` (columns: `t  cL  cD`)

## Project Structure

```
JuLattice/
├── JuLattice_main.jl   # Entry point — user settings and main loop
├── src/
│   ├── Kernel.jl       # BGK collision + streaming
│   ├── GridSetup.jl    # Domain and node classification
│   ├── Plotter.jl      # GLMakie real-time plots
│   └── IO.jl           # Force output and logging
├── Project.toml
└── README.md
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Other Branches

### [`MRT_D2Q9_gridRefinement`](https://github.com/DomowojS/JuLattice/tree/MRT_D2Q9_gridRefinement)
Two-level grid refinement with MRT collision. Use this for higher-resolution studies.

### [`JuLattice_for_teaching`](https://github.com/DomowojS/JuLattice/tree/JuLattice_for_teaching)
Hard-coded, minimal version intended for educational use.
