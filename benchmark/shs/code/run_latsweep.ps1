# run_latsweep.ps1 -- latitude stress test (341-point meridian, data/points_latsweep.csv).
#
# Thin wrapper: the latsweep set is run by run_all.ps1 with the same engines, repeats and outputs
# (results/raw, then results/csv and the figures). Same parameters as run_all.ps1, minus -Sets.
# Tool folders come from the environment variables GRAFLAB_DIR, MATLAB_EXE and MSYS2_UCRT64_BIN
# (see benchmark/shs/README.md), or from the parameters below.
#     pwsh benchmark/shs/code/run_latsweep.ps1 -GrafLabDir <folder with GrafLab.m>
[CmdletBinding()]
param(
  [int]$Repeat = 3,
  [string[]]$SkipEngines = @(),
  [switch]$OnlyGsh,
  [string]$GrafLabDir,
  [string]$Matlab,
  [string]$RawDir,
  [string]$ImgDir,
  [int]$Nmax = 2190,
  [int]$MaxPoints = 0,
  [string]$Msys2Bin,
  [switch]$Append,
  [switch]$DryRun
)
$ErrorActionPreference = 'Stop'
# Only the parameters the caller actually gave are passed on; run_all.ps1 applies its own defaults.
$forward = @{ Sets = @('latsweep') }
foreach ($k in $PSBoundParameters.Keys) { $forward[$k] = $PSBoundParameters[$k] }
& "$PSScriptRoot\run_all.ps1" @forward
