# JuLattice

**JuLattice** is a 2D Lattice Boltzmann Method (LBM) solver written in Julia. It implements the D2Q9 lattice scheme and uses the Bhatnagar-Gross-Krook (BGK) approximation for collision, where a single relaxation time (defined by the user) governs the simulation dynamics.

## Features

- D2Q9 LBM implementation
- BGK collision model with user-defined relaxation time
- Configurable simulation parameters via JSON configuration files
- Real-time plotting of simulation results
- Written in pure Julia

## Requirements

- Julia ≥ **1.9.0**

## Getting Started

### 1. Clone the repository

Open a terminal and run:

```bash
git clone https://github.com/yourusername/JuLattice.git
cd JuLattice
