# Startup: 
#   Start-Process powershell -ArgumentList "-File .\run_gridStudy_disc.ps1" -WindowStyle Normal


$mutex = [System.Threading.Mutex]::new($false, "JuLattice_queue")
if (-not $mutex.WaitOne(0)) {
    Write-Host "ERROR: Another queue instance is already running. Exiting." -ForegroundColor Red
    exit 1
}

$julia   = "julia"
$project = "."
$threads = 16

# JULIA_EXCLUSIVE: no CPU yielding to other OS processes
$env:JULIA_EXCLUSIVE = "1"

$runfiles = @(
    "gridStudy2_Ddx20_Re2760_Ma0.1_Ct0.65.jl",
    "gridStudy2_Ddx27_Re2760_Ma0.1_Ct0.65.jl"
)

Write-Host "================================" -ForegroundColor Cyan
Write-Host "  JuLattice Simulation Queue" -ForegroundColor Cyan
Write-Host "  Threads: $threads  |  JULIA_EXCLUSIVE=1" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan

$total = $runfiles.Count
$run   = 1

foreach ($file in $runfiles) {
    $logFile = $file -replace '\.jl$', '.log'

    Start-Transcript -Path $logFile -Append

    Write-Host ""
    Write-Host "[$run/$total] Starting: $file  |  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan

    $start = Get-Date
    & $julia --project=$project -t $threads $file
    $exitCode = $LASTEXITCODE
    $elapsed  = (Get-Date) - $start

    Write-Host "[$run/$total] Done (exit code $exitCode) in $($elapsed.ToString('hh\:mm\:ss'))  |  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Green

    Stop-Transcript

    if ($exitCode -ne 0) {
        Write-Host "ERROR: $file failed  stopping queue." -ForegroundColor Red
        $mutex.ReleaseMutex()
        exit 1
    }

    $run++
}

Write-Host ""
Write-Host "=================================" -ForegroundColor Cyan
Write-Host "  All $total simulations complete" -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan

$mutex.ReleaseMutex()