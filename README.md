# JuLattice — MRT D2Q9 with Grid Refinement

![Demo](media/VelMag_Rectangle.gif)

**JuLattice** is a 2D Lattice Boltzmann Method (LBM) solver written in Julia. This branch implements a single D2Q9 grid with MRT collision operator.

## Features

- D2Q9 MRT collision with configurable BGK and acoustic relaxation rates
- Bouzidi curved boundary conditions for immersed objects
- Rectangular obstacle with arbitrary angle
- Lift and drag force output (cL, cD) at every time step
- Contour plots of vorticity, u-velocity, v-velocity via GLMakie

## Requirements

![Julia version](https://img.shields.io/badge/julia-1.9%2B-blue)

**Dependencies** (all defined in `Project.toml`):
- `GLMakie`
- `MeshGrid`
- `Revise`

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

The fine box is snapped to the nearest integer number of coarse cells so that the coarse-fine coupling nodes align correctly.

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

### [`main`](https://github.com/DomowojS/JuLattice/tree/main)
Single-grid D2Q9 BGK solver with JSON configuration files. Use this for simpler setups.

### [`MRT_D2Q9`](https://github.com/DomowojS/JuLattice/tree/MRT_D2Q9_gridRefinement)
Single- grid D2Q9 solver with MRT collision.

### [`JuLattice_for_teaching`](https://github.com/DomowojS/JuLattice/tree/JuLattice_for_teaching)
Hard-coded, minimal version intended for educational use.
