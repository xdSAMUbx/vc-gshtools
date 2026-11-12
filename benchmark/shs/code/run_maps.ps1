# run_maps.ps1 -- spatial maps for the final comparison (not part of the timed benchmark).
# Synthesises V, T and dg on data/points_grid.csv with GrafLab (sfc, msfc, era) and gsh (fast), one repetition,
# into results/raw/maps, then draws img/map_colorado_grid_graflab_vs_gsh.png.
#   pwsh benchmark/shs/code/run_maps.ps1 -GrafLabDir <folder with GrafLab.m>
[CmdletBinding()]
param([string]$GrafLabDir = $env:GRAFLAB_DIR, [string]$Matlab = $env:MATLAB_EXE, [int]$Nmax = 2190)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path "$PSScriptRoot\..\..\..").Path
$bench = "$root\benchmark\shs"
$gfc = "$bench\data\XGM2019e_2190_tide_free.gfc"
$pts = "$bench\data\points_grid.csv"
$raw = "$bench\results\raw\maps"
$py = "$root\.venv\Scripts\python.exe"
if (-not $Matlab) { $Matlab = 'matlab' }
if (-not (Test-Path $pts)) { & $py "$bench\code\python\make_map_grid.py" }
New-Item -ItemType Directory -Force $raw | Out-Null
& "$root\build\clang\bin\gsh_shs_bench.exe" --gfc $gfc --points $pts --recurrence fast --nmax $Nmax `
  --out "$raw\gsh_fast_grid_t12.csv" --timings "$raw\timings.csv" --set grid --threads 12 --repeat 1
if ($LASTEXITCODE) { throw "gsh failed" }
foreach ($rec in 'sfc', 'msfc', 'era') {
  $cmd = "addpath('$bench\code\matlab'); run_graflab('$gfc','$pts',$Nmax,'$rec','grid','$raw',1,'$GrafLabDir')"
  & $Matlab -batch $cmd
  if ($LASTEXITCODE) { throw "GrafLab $rec failed" }
}
& $py "$bench\code\python\plot_maps.py" --raw $raw --img "$bench\img"
