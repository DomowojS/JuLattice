# ============================================================
# Grid study — Actuator disc
# Runs all runfiles back to back on ALL available CPU cores
# ============================================================

$julia   = "julia"
$project = "."
$threads = 16

$runfiles = @(
    "runfile_disc_gridStudyDdx=8.jl",
    "runfile_disc_gridStudyDdx=10.jl",
    "runfile_disc_gridStudyDdx=14.jl",
    "runfile_disc_gridStudyDdx=20.jl",
    "runfile_disc_gridStudyDdx=27.jl",
    "runfile_disc_gridStudyDdx=40.jl",
    "runfile_disc_gridStudyDdx=54.jl"
)

Write-Host "Using $threads threads (all logical cores)" -ForegroundColor Yellow

$total = $runfiles.Count
$run   = 1

foreach ($file in $runfiles) {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "  Run $run / $total  —  $file" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan

    $start = Get-Date
    & $julia --project=$project -t $threads $file
    $elapsed = (Get-Date) - $start

    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERROR: $file failed (exit code $LASTEXITCODE)" -ForegroundColor Red
        exit 1
    }

    Write-Host "  Done in $($elapsed.ToString('hh\:mm\:ss'))" -ForegroundColor Green
    $run++
}

Write-Host ""
Write-Host "All $total runs completed." -ForegroundColor Green