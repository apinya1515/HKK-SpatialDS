# run_step23_nosleep.ps1
# Wrapper to run step23_spatial_compare.R with sleep/hibernate disabled
# Restores original power settings on completion (or Ctrl+C)

Write-Host "=== Saving current power settings ===" -ForegroundColor Cyan

# Save current sleep/hibernate timeouts (in seconds, from registry)
$acSleep = powercfg /query SCHEME_CURRENT SUB_SLEEP STANDBYIDLE | Select-String "Current AC Power Setting Index"
$dcSleep = powercfg /query SCHEME_CURRENT SUB_SLEEP STANDBYIDLE | Select-String "Current DC Power Setting Index"

# Extract hex values and convert to minutes
$acSleepMin = if ($acSleep -match '0x([0-9a-fA-F]+)') { [convert]::ToInt32($matches[1], 16) / 60 } else { 30 }
$dcSleepMin = if ($dcSleep -match '0x([0-9a-fA-F]+)') { [convert]::ToInt32($matches[1], 16) / 60 } else { 15 }

Write-Host "  AC sleep timeout: $acSleepMin min"
Write-Host "  DC sleep timeout: $dcSleepMin min"

# Function to restore settings
function Restore-PowerSettings {
    Write-Host "`n=== Restoring power settings ===" -ForegroundColor Yellow
    powercfg /change standby-timeout-ac $acSleepMin
    powercfg /change standby-timeout-dc $dcSleepMin
    powercfg /change hibernate-timeout-ac 60
    powercfg /change hibernate-timeout-dc 30
    Write-Host "  Restored AC sleep: $acSleepMin min, DC sleep: $dcSleepMin min"
}

# Register cleanup on script exit
Register-EngineEvent -SourceIdentifier PowerShell.Exiting -Action { Restore-PowerSettings } | Out-Null

Write-Host "`n=== Disabling sleep and hibernate ===" -ForegroundColor Green
powercfg /change standby-timeout-ac 0
powercfg /change standby-timeout-dc 0
powercfg /change hibernate-timeout-ac 0
powercfg /change hibernate-timeout-dc 0
Write-Host "  Sleep and hibernate DISABLED"

Write-Host "`n=== Starting step23_spatial_compare.R ===" -ForegroundColor Green
Write-Host "  Start time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Write-Host "  Log file: step23_run.log"
Write-Host ""

# Run R script with direct file redirection (more robust than Tee-Object pipes)
# Use Start-Process to avoid broken pipe issues
$logFile = Join-Path $PSScriptRoot "step23_run.log"
$proc = Start-Process -FilePath "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" `
    -ArgumentList "step23_spatial_compare.R" `
    -WorkingDirectory $PSScriptRoot `
    -RedirectStandardOutput $logFile `
    -RedirectStandardError (Join-Path $PSScriptRoot "step23_run_err.log") `
    -NoNewWindow -PassThru

Write-Host "  R process started (PID: $($proc.Id))"
Write-Host "  Monitor progress with: Get-Content step23_run.log -Tail 5 -Wait"
Write-Host ""

# Wait for completion
$proc.WaitForExit()
$exitCode = $proc.ExitCode

Write-Host "`n=== R script finished ===" -ForegroundColor $(if ($exitCode -eq 0) { "Green" } else { "Red" })
Write-Host "  Exit code: $exitCode"
Write-Host "  End time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"

# Show last few lines of log
Write-Host "`n=== Last 10 lines of log ===" -ForegroundColor Cyan
Get-Content $logFile -Tail 10

# Check for errors
$errLog = Join-Path $PSScriptRoot "step23_run_err.log"
if ((Test-Path $errLog) -and (Get-Item $errLog).Length -gt 0) {
    Write-Host "`n=== Error log (last 10 lines) ===" -ForegroundColor Red
    Get-Content $errLog -Tail 10
}

# Restore power settings
Restore-PowerSettings

Write-Host "`nDone!" -ForegroundColor Green
