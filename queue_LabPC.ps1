# run_gridStudy_Windows.ps1
# Run with: Start-Process powershell -ArgumentList "-File run_gridStudy_Windows.ps1" -WindowStyle Normal
# or just: .\run_gridStudy_Windows.ps1

Set-Location $PSScriptRoot
$THREADS = [Environment]::ProcessorCount
New-Item -ItemType Directory -Force -Path "logs" | Out-Null

Write-Host "================================"
Write-Host "  JuLattice Simulation Queue"
Write-Host "================================"

Write-Host "[1/2] Starting: Grid Study D/dx=15 | VF-outflow"
julia -t $THREADS run_gridStudy_Ddeltax15_VF-outflow.jl `
    1>"logs\stdout_Ddx15.log" 2>"logs\stderr_Ddx15.log"
Write-Host "[1/2] Done (exit code $LASTEXITCODE): $(Get-Date)"

Write-Host "[2/2] Starting: Grid Study D/dx=20 | VF-outflow"
julia -t $THREADS run_gridStudy_Ddeltax20_VF-outflow.jl `
    1>"logs\stdout_Ddx20.log" 2>"logs\stderr_Ddx20.log"
Write-Host "[2/2] Done (exit code $LASTEXITCODE): $(Date)"

Write-Host "================================="
Write-Host "  All simulations complete"
Write-Host "================================="