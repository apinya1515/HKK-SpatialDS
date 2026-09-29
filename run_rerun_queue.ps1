# run_rerun_queue.ps1
# Re-runs every model with Delta_wAIC < 2 in Results/Final_Model_Comparison.csv through
# rerun_models.R ($Chains chain processes per model), keeping at most $MaxModels models running.
# Rank 1 models go first. Progress: Results/MCMC_rerun/queue.log and <SP>_Rank<N>/progress.log
#
# Start detached:
#   Start-Process powershell -ArgumentList '-ExecutionPolicy Bypass -File run_rerun_queue.ps1' -WindowStyle Hidden

# Each chain process commits ~2.9 GB (spatial); 3 chains x 4 models = 12 processes = ~35 GB.
# -AllCandidates: every model in Final_Model_Comparison.csv (for re-scoring WAIC), not only Delta_wAIC < 2
param([int]$MaxModels = 4, [int]$Chains = 3, [switch]$AllCandidates)

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$R = 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
$outDir = 'Results\MCMC_rerun'
New-Item -ItemType Directory -Force $outDir | Out-Null
$log = Join-Path $outDir 'queue.log'
function Log($m) { "[{0}] {1}" -f (Get-Date -Format 'MM-dd HH:mm:ss'), $m | Add-Content $log }

$codes = @{ 'Banteng' = 'BTG'; 'Sambar deer' = 'SBR'; 'Gaur' = 'GAR'; 'Muntjac' = 'MJK'; 'Wild boar' = 'PIG' }
$queue = Import-Csv 'Results\Final_Model_Comparison.csv' |
  Where-Object { $AllCandidates -or [double]$_.Delta_wAIC -lt 2 } |
  ForEach-Object { [pscustomobject]@{ Sp = $codes[$_.Species]; Rank = [int]$_.Rank; Type = $_.Type } } |
  Sort-Object @{ Expression = { $_.Type -ne 'Spatial' } }, Rank, Sp  # long spatial runs first
Log ("Queue ({0} models): {1}" -f $queue.Count, (($queue | ForEach-Object { "$($_.Sp)_R$($_.Rank)" }) -join ', '))

function Running-Models {
  Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
    Where-Object { $_.CommandLine -match 'rerun_models\.R"?\s+"?(\w+)"?\s+"?(\d+)"?' } |
    ForEach-Object { "{0}_Rank{1}" -f $Matches[1], $Matches[2] } | Sort-Object -Unique
}

# Safe to restart: models already running or finished (final RDS present) are skipped
$running0 = @(Running-Models)
# Iteration cap for models started while the cap in rerun_models.R was higher (150,000): at the
# 100,000-iteration check (chunk 10 of 10,000) such a model writes CONTINUE instead of STOP. Stop its
# chains and save the capped result (logged NOT converged) with check-only mode + FORCE_SAVE.
$capChunk = 10
function Enforce-Cap {
  foreach ($tag in @(Running-Models)) {
    $dec = "$outDir\$tag\decision$('{0:D2}' -f $capChunk).txt"
    if ((Test-Path $dec) -and ((Get-Content $dec -TotalCount 1) -eq 'CONTINUE') -and -not (Test-Path "$outDir\MCMC_Samples_$tag.rds")) {
      $sp, $rk = $tag -replace 'Rank', '' -split '_'
      Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
        Where-Object { $_.CommandLine -match "rerun_models\.R`"?\s+`"?$sp`"?\s+`"?$rk`"?\s" } |
        ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
      $env:FORCE_SAVE = '1'
      $p = Start-Process -FilePath $R -ArgumentList "rerun_models.R $sp $rk check" -NoNewWindow -Wait -PassThru `
        -RedirectStandardOutput "$outDir\$tag-capcheck.out" -RedirectStandardError "$outDir\$tag-capcheck.err"
      Remove-Item Env:FORCE_SAVE
      Log "$tag reached 100,000 iterations without passing: stopped, capped result saved (exit $($p.ExitCode))"
    }
  }
}

# Rhat threshold for models started while rerun_models.R used a stricter one (1.03): once a model's
# latest check has >= 30,000 iterations and max Rhat < 1.045, re-check its chunks with check-only
# mode (current rule). Only if that saves the final RDS are the model's chains stopped.
$rhatMax = 1.045
$rhatTried = @{}
function Enforce-Rhat {
  foreach ($tag in @(Running-Models)) {
    $lg = "$outDir\$tag\convergence_log.csv"
    if (-not (Test-Path $lg) -or (Test-Path "$outDir\MCMC_Samples_$tag.rds")) { continue }
    $last = (Import-Csv $lg)[-1]
    $car = if ($last.max_rhat_car -and $last.max_rhat_car -ne 'NA') { [double]$last.max_rhat_car } else { 0 }
    $worst = [math]::Max([double]$last.max_rhat_noncar, $car)
    $key = "$tag@$($last.chunk)"
    if ([int]$last.iterations_per_chain -ge 30000 -and $worst -lt $rhatMax -and -not $rhatTried.ContainsKey($key)) {
      $rhatTried[$key] = $true
      $sp, $rk = $tag -replace 'Rank', '' -split '_'
      $p = Start-Process -FilePath $R -ArgumentList "rerun_models.R $sp $rk check" -NoNewWindow -Wait -PassThru `
        -RedirectStandardOutput "$outDir\$tag-rhatcheck.out" -RedirectStandardError "$outDir\$tag-rhatcheck.err"
      if (Test-Path "$outDir\MCMC_Samples_$tag.rds") {
        Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" |
          Where-Object { $_.CommandLine -match "rerun_models\.R`"?\s+`"?$sp`"?\s+`"?$rk`"?\s" } |
          ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
        Log "$tag passed Rhat < $rhatMax at $($last.iterations_per_chain) iterations: final saved, chains stopped"
      } else {
        Log "$tag re-check under Rhat < $rhatMax did not pass (exit $($p.ExitCode)); chains keep running"
      }
    }
  }
}

$pending = [System.Collections.ArrayList]@($queue | Where-Object {
  $tag = "$($_.Sp)_Rank$($_.Rank)"
  -not ($running0 -contains $tag) -and -not (Test-Path "$outDir\MCMC_Samples_$tag.rds")
})
Log ("Pending after skipping running/finished: {0}" -f (($pending | ForEach-Object { "$($_.Sp)_R$($_.Rank)" }) -join ', '))
while ($pending.Count -gt 0 -or @(Running-Models).Count -gt 0) {
  Enforce-Cap
  Enforce-Rhat
  while ($pending.Count -gt 0 -and @(Running-Models).Count -lt $MaxModels) {
    $m = $pending[0]; $pending.RemoveAt(0)
    $tag = "$($m.Sp)_Rank$($m.Rank)"
    # Clear chunk/decision files left by an interrupted run of this model
    if (Test-Path "$outDir\$tag") { Remove-Item -Recurse -Force "$outDir\$tag" }
    foreach ($c in 1..$Chains) {
      Start-Process -FilePath $R -ArgumentList "rerun_models.R $($m.Sp) $($m.Rank) $c" -WindowStyle Hidden `
        -RedirectStandardOutput "$outDir\$tag-chain$c.out" -RedirectStandardError "$outDir\$tag-chain$c.err" | Out-Null
    }
    Log "started $tag ($($m.Type))"
    Start-Sleep -Seconds 20
  }
  Start-Sleep -Seconds 60
}
Log 'all models finished'
