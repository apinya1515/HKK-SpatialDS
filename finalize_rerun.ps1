# finalize_rerun.ps1
# Post-run steps after run_rerun_queue.ps1 finishes. Stops at the first failing step.
#   1. Verify every Delta_wAIC < 2 model has a final RDS and converged
#   2. Archive Results/MCMC/*.rds -> Results/MCMC_archive_<date>/, install the new RDS files
#      (+ MCMC_Samples_<SP>.rds alias = Rank 1, as run_compiled_posteriors.R did)
#   3. Regenerate every Results/ folder from the new samples, in dependency order
#   4. Post-hoc WAIC, MCMC settings record, old-vs-new estimate comparison
# Log: Results/MCMC_rerun/finalize.log

$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$R = 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
$rerun = 'Results\MCMC_rerun'
$log = "$rerun\finalize.log"
function Log($m) { $l = "[{0}] {1}" -f (Get-Date -Format 'MM-dd HH:mm:ss'), $m; $l | Add-Content $log; Write-Host $l }

# ---- 1. Verify
$codes = @{ 'Banteng' = 'BTG'; 'Sambar deer' = 'SBR'; 'Gaur' = 'GAR'; 'Muntjac' = 'MJK'; 'Wild boar' = 'PIG' }
$models = Import-Csv 'Results\Final_Model_Comparison.csv' | Where-Object { [double]$_.Delta_wAIC -lt 2 } |
  ForEach-Object { "{0}_Rank{1}" -f $codes[$_.Species], $_.Rank }
$bad = @()
foreach ($t in $models) {
  if (-not (Test-Path "$rerun\MCMC_Samples_$t.rds")) { $bad += "$t (no final RDS)"; continue }
  $last = (Import-Csv "$rerun\$t\convergence_log.csv")[-1]
  if ($last.pass -ne 'TRUE') { $bad += "$t (not converged)" }
}
if ($bad.Count) { Log ("STOP - not ready: " + ($bad -join '; ')); exit 1 }
Log "All $($models.Count) models converged."

# ---- 2. Archive old samples, install new ones
$archive = 'Results\MCMC_archive_' + (Get-Date -Format 'yyyy-MM-dd')
New-Item -ItemType Directory -Force $archive | Out-Null
Get-ChildItem Results\MCMC -Filter *.rds | Move-Item -Destination $archive
foreach ($t in $models) { Copy-Item "$rerun\MCMC_Samples_$t.rds" Results\MCMC\ }
foreach ($sp in 'BTG', 'GAR', 'MJK', 'PIG', 'SBR') {
  Copy-Item "Results\MCMC\MCMC_Samples_${sp}_Rank1.rds" "Results\MCMC\MCMC_Samples_$sp.rds"
}
Log "Old samples archived to $archive; new samples installed in Results\MCMC."

# ---- 3 & 4. Regenerate outputs (each script in a fresh R session)
$steps = @(
  @('generate_model_summaries.R'),          # model_summary (xlsx, Master csv used below)
  @('create_tables.R'),                     # tables
  @('generate_detection_outputs.R'),        # detection
  @('generate_car_maps.R'),                 # CAR_map
  @('generate_covariate_boxplots.R'),       # Covariates
  @('update_maps_covariates_averaging.R'),  # Maps + Density + Averaging
  @('generate_presentation_outputs.R'),     # Community density, heatmap, synthesis tables
  @('generate_traceplots.R'),               # traceplots
  @('compute_waic_posthoc.R'),              # tables/WAIC_Rerun_Comparison.csv
  @('summarize_rerun.R', $archive)          # model_summary/Rerun_*.csv
)
foreach ($s in $steps) {
  $name = $s[0]
  Log "Running $name ..."
  $out = "$rerun\finalize_$($name -replace '\.R$','')"
  # Start-Process: R writes normal messages to stderr, which PowerShell 5.1 would treat as errors
  $p = Start-Process -FilePath $R -ArgumentList $s -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput "$out.out" -RedirectStandardError "$out.err"
  if ($p.ExitCode -ne 0) { Log "FAILED: $name (exit $($p.ExitCode)) - see $out.err"; exit 1 }
  Log "  done: $name"
}
Log 'Finalize complete.'
