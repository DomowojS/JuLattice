# Run from: C:\...\JuLattice\

Write-Host "Starting simulation 1: Grid Study D/dx=30"
julia run_gridStudy_Ddeltax30.jl

Write-Host "Starting simulation 2: Ma Study D/dx=40 Ma=0.05"
julia run_MaStudy_Ddeltax40.jl

Write-Host "All simulations complete."