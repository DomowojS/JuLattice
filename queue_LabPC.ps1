# Run from: C:\...\JuLattice\

Write-Host "Starting simulation 1: Grid Study D/dx=25"
julia run_gridStudy_Ddeltax25.jl

Write-Host "Starting simulation 2: Domain Extension D/dx=20"
julia run_domainExtension_Ddeltax20.jl

Write-Host "All simulations complete."
