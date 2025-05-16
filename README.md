# JuLattice
![Demo](media/Vorticity_Cylinder.gif)
**JuLattice** is a 2D Lattice Boltzmann Method (LBM) solver written in Julia. It implements the D2Q9 lattice scheme and uses the Bhatnagar-Gross-Krook (BGK) approximation for collision, where a single relaxation time (defined by the user) governs the simulation dynamics.

## Features

- D2Q9 LBM implementation
- BGK collision model with user-defined relaxation time
- Configurable simulation parameters via JSON configuration files
- Real-time plotting of simulation results
- Written in pure Julia

## Requirements

- ![Julia version](https://img.shields.io/badge/julia-1.9%2B-blue)

## Getting Started

### 1. Clone the repository

Open a terminal and run:

```bash
git clone https://github.com/yourusername/JuLattice.git](https://github.com/DomowojS/JuLattice.git
cd JuLattice
```
## 2. Activate the project environment
Start Julia inside the project directory and activate the environment using the Julia package manager:
```julia
using Pkg
Pkg.activate(".")
```
This will activate the local environment and make all required packages available.

## 3. Run a simulation
To run JuLattice, first include the main script:
```julia
include("JuLattice_main.jl")
```
Then call the run function with the path to a configuration file, e.g.:
```julia
JuLattice.run("config/cylinder.json")
```
This will execute the simulation defined in the config file and start real-time plotting.
You can replace "config/cylinder.json" with any other configuration file in the /config directory.

## Configuration
All simulation settings are defined in JSON files located in the config/ directory. Each configuration file allows you to specify:

- Domain size and resolution
- Boundary conditions
- Fluid properties (e.g., relaxation time)
- Obstacle geometry (if any)
- Simulation duration and output preferences

The user only needs to modify these configuration files to define new simulation setups—no code modification is required.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Other Branches
### [`JuLattice_for_teaching`](https://github.com/DomowojS/JuLattice/tree/JuLattice_for_teaching)

This branch contains a simplified, hard-coded version of JuLattice used for educational purposes.  
It is less modular but more straightforward, intended to demonstrate the core Lattice Boltzmann logic without configuration overhead.

> **Note**: For a full-featured version with configuration files and plotting, use the [`main`](https://github.com/DomowojS/JuLattice/tree/main) branch.
