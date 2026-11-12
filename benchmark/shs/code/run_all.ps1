# run_all.ps1 -- reproduces the spherical-harmonic synthesis benchmark, sequentially, on one machine.
#
# Run from anywhere (PowerShell 7). The final, timed run is ONE command (see benchmark/shs/README.md):
#     pwsh benchmark/shs/code/run_all.ps1 -Sets colorado,latsweep -GrafLabDir <folder with GrafLab.m>
#
# Principle: every third-party tool is run in its official, as-distributed form. No parallelism is added by
# this benchmark to a tool that does not have it. Run matrix per point set (3 repetitions everywhere):
#   gsh        fast at 12 and 1 threads; msfc, sfc, era at 12 threads        parallel_mode library (OpenMP)
#   GrafLab    sfc, msfc, era, MATLAB default threads                        parallel_mode library (MATLAB)
#   pyshtools  1 thread (serial library)                                     parallel_mode none
#   shtools    1 thread (serial driver loop)                                 parallel_mode none
#   pyharm     12 and 1 threads (one batched CHarm call)                     parallel_mode library
#   charm      --parallel charm at 12 and 1 threads (one batched call)       parallel_mode library
# Each run computes V, T and dg together. Outputs: raw\<engine>_<rec>_<set>_t<N>.csv (+ .md), raw\timings.csv.
# Then consolidate.py writes csv\<folder>\res_value.csv, res_time.csv; compare.py writes
# csv\<folder>\comp_value.csv, raw\summary_tables.md and the figures in img\.
#
# Parameters:
#   -Sets          colorado (default), latsweep, or colorado,latsweep
#   -Repeat        timed repetitions per run (default 3)
#   -SkipEngines   any of gsh, graflab, pyshtools, pyharm, shtools, charm
#   -OnlyGsh       run only the gsh stages; keep the frozen references already in raw\ (GrafLab, pyshtools,
#                  pyharm, shtools, charm), consolidate, compare, and check gsh against the stored comparison
#   -GrafLabDir    folder with GrafLab.m. Default: environment variable GRAFLAB_DIR. Required for GrafLab.
#   -Matlab        matlab executable. Default: environment variable MATLAB_EXE, else the one on PATH.
#   -Msys2Bin      MSYS2 UCRT64 bin folder (runtime DLLs of the CHarm and SHTOOLS drivers). Default:
#                  environment variable MSYS2_UCRT64_BIN. Required for charm and shtools.
#   -RawDir        raw outputs folder (default benchmark/shs/results/raw); csv\ is written next to it
#   -ImgDir        figure folder (default benchmark/shs/img)
#   -Nmax          maximum degree (default 2190). Smoke tests: nmax <= 60.
#   -MaxPoints     smoke tests only: first N points of each set (default 0 = all points)
#   -Append        keep an existing raw\timings.csv and add rows (needed to resume)
#   -DryRun        print the plan and every command line, run nothing
#
# Prerequisites: benchmark/shs/README.md (gsh build, Python venv, MATLAB + GrafLab, CHarm and SHTOOLS drivers,
# decompressed model). The machine must be otherwise idle: wall times are only meaningful on a quiet system.
[CmdletBinding()]
param(
  [string[]]$Sets = @('colorado'),
  [int]$Repeat = 3,
  [string[]]$SkipEngines = @(),
  [switch]$OnlyGsh,
  [string]$GrafLabDir = $env:GRAFLAB_DIR,
  [string]$Matlab = $env:MATLAB_EXE,
  [string]$RawDir = '',
  [string]$ImgDir = '',
  [int]$Nmax = 2190,
  [int]$MaxPoints = 0,
  [string]$Msys2Bin = $env:MSYS2_UCRT64_BIN,
  [switch]$Append,
  [switch]$DryRun
)
$ErrorActionPreference = 'Stop'

$root = (Resolve-Path "$PSScriptRoot\..\..\..").Path
$bench = "$root\benchmark\shs"
$data = "$bench\data"
$raw = if ($RawDir) { [IO.Path]::GetFullPath($RawDir) } else { "$bench\results\raw" }
$img = if ($ImgDir) { [IO.Path]::GetFullPath($ImgDir) } else { "$bench\img" }
$resultsDir = Split-Path -Parent $raw
$csvDir = Join-Path $resultsDir 'csv'
$gfc = "$data\XGM2019e_2190_tide_free.gfc"
$gsh = "$root\build\clang\bin\gsh_shs_bench.exe"
$charmBin = "$root\build\external\charm-driver"
$charmExe = "$charmBin\charm_shs_bench.exe"
$shtoolsBin = "$root\build\external\shtools_bench"
$shtoolsExe = "$shtoolsBin\shtools_shs_bench.exe"
$py = "$root\.venv\Scripts\python.exe"
$threads = [Environment]::ProcessorCount
$timings = "$raw\timings.csv"
$knownEngines = @('gsh', 'graflab', 'pyshtools', 'pyharm', 'shtools', 'charm')
$knownSets = @('colorado', 'latsweep')

# Accept "-Sets a,b" and "-Sets a, b" under -File (where the comma list arrives as one string).
$Sets = @($Sets | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
$SkipEngines = @($SkipEngines | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
foreach ($e in $SkipEngines) { if ($knownEngines -notcontains $e) { throw "Unknown engine in -SkipEngines: $e" } }
foreach ($s in $Sets) { if ($knownSets -notcontains $s) { throw "Unknown set in -Sets: $s" } }
if ($Repeat -lt 1) { throw '-Repeat must be >= 1' }
if ($Nmax -lt 2) { throw '-Nmax must be >= 2' }
$skip = @{}
foreach ($e in $knownEngines) { $skip[$e] = $SkipEngines -contains $e }
if ($OnlyGsh) { foreach ($e in $knownEngines) { if ($e -ne 'gsh') { $skip[$e] = $true } } ; $skip['gsh'] = $false }

# Tool locations that are not given explicitly: never hard-coded to a user's machine.
if (-not $skip['graflab']) {
  if (-not $GrafLabDir) {
    $msg = 'GrafLab folder not set: pass -GrafLabDir <folder with GrafLab.m> or set GRAFLAB_DIR. GrafLab is not redistributed (see benchmark/shs/README.md).'
    if ($DryRun) { Write-Host "NOTE: $msg" -ForegroundColor Yellow } else { throw $msg }
  }
  if (-not $Matlab) {
    $found = Get-Command matlab -ErrorAction SilentlyContinue
    if ($found) { $Matlab = $found.Source }
    elseif (-not $DryRun) { throw 'MATLAB not found: pass -Matlab <matlab.exe> or set MATLAB_EXE.' }
  }
}
if ((-not $skip['charm'] -or -not $skip['shtools']) -and -not $Msys2Bin) {
  $msg = 'MSYS2 UCRT64 bin folder not set: pass -Msys2Bin <...>\ucrt64\bin or set MSYS2_UCRT64_BIN (needed by the CHarm and SHTOOLS drivers).'
  if ($DryRun) { Write-Host "NOTE: $msg" -ForegroundColor Yellow } else { throw $msg }
}

$total = [Diagnostics.Stopwatch]::StartNew()
function Stamp { Get-Date -Format 'HH:mm:ss' }
function Fmt([TimeSpan]$t) { '{0:hh\:mm\:ss}' -f $t }
function Step([string]$msg) { Write-Host "[$(Stamp)] $msg" }

# Paths handed to the drivers are relative to the repository root when inside it (no local paths in outputs).
function Rel([string]$p) {
  $full = [IO.Path]::GetFullPath($p)
  $base = $root.TrimEnd('\') + '\'
  if ($full.StartsWith($base, [StringComparison]::OrdinalIgnoreCase)) { return $full.Substring($base.Length) }
  return $full
}

# Publication check: no local paths, user folder names or personal names in the outputs.
$privatePattern = 'C:\\Users|C:/Users|Downloads|samue|Calderon|Documents\\programacion'
function Find-PrivateText {
  $files = @()
  if (Test-Path $csvDir) { $files += Get-ChildItem $csvDir -Recurse -File }
  $files += Get-ChildItem $resultsDir -Filter '*.md' -File -ErrorAction SilentlyContinue
  if (Test-Path $img) { $files += Get-ChildItem $img -File }
  if (-not $files) { return @() }
  return @(Select-String -Path $files.FullName -Pattern $privatePattern -ErrorAction SilentlyContinue)
}

# ---------------------------------------------------------------------------------------------
# 0) Machine information and idle warning (before anything runs)
# ---------------------------------------------------------------------------------------------
Write-Host '============================================================'
Write-Host ' SHS benchmark: run'
Write-Host '============================================================'
try {
  $cpu = Get-CimInstance Win32_Processor | Select-Object -First 1
  $os = Get-CimInstance Win32_OperatingSystem
  $cs = Get-CimInstance Win32_ComputerSystem
  Write-Host ("Machine : {0} ({1} cores, {2} logical processors)" -f $cpu.Name.Trim(), $cpu.NumberOfCores,
    $cpu.NumberOfLogicalProcessors)
  Write-Host ("RAM     : {0:N1} GiB   OS: {1} {2}" -f ($cs.TotalPhysicalMemory / 1GB), $os.Caption, $os.Version)
  Write-Host ("CPU load now: {0} %   threads used by the parallel runs: {1}" -f $cpu.LoadPercentage, $threads)
  $busy = Get-CimInstance Win32_PerfFormattedData_PerfProc_Process -ErrorAction SilentlyContinue |
    Where-Object { $_.Name -ne '_Total' -and $_.Name -ne 'Idle' } |
    Sort-Object PercentProcessorTime -Descending | Select-Object -First 5
  if ($busy) {
    Write-Host 'Busiest processes right now:'
    foreach ($b in $busy) { Write-Host ("  {0,-30} {1,6} %" -f $b.Name, $b.PercentProcessorTime) }
  }
} catch {
  Write-Host "Machine information not available: $($_.Exception.Message)"
}
$mode = if ($OnlyGsh) { 'OnlyGsh (gsh stages only, frozen references kept)' } else { 'full' }
Write-Host "PowerShell $($PSVersionTable.PSVersion); mode: $mode; sets: $($Sets -join ', '); repeat: $Repeat; nmax: $Nmax; points: $(if ($MaxPoints -gt 0) { "first $MaxPoints (smoke)" } else { 'all' })"
Write-Host "Raw outputs : $raw"
Write-Host "CSV outputs : $csvDir"
Write-Host "Figures     : $img"
Write-Host ''
Write-Host 'WARNING: the machine must be otherwise idle. Close browsers, IDEs, sync clients, MATLAB GUI and any build.' -ForegroundColor Yellow
Write-Host 'WARNING: wall times of this run are only comparable if nothing else competes for the CPU.' -ForegroundColor Yellow
$pre = Find-PrivateText
if ($pre) {
  Write-Host "WARNING: $($pre.Count) line(s) with local paths or personal names already exist in the outputs (checked at the end; remove the files): " -ForegroundColor Yellow
  foreach ($m in ($pre | Select-Object -First 5)) { Write-Host "  $($m.Path):$($m.LineNumber)" }
}
if (-not $DryRun) { Write-Host 'Starting in 10 seconds (Ctrl+C to abort)...'; Start-Sleep -Seconds 10 }

# ---------------------------------------------------------------------------------------------
# 1) Preflight: every required file must exist before the first run starts
# ---------------------------------------------------------------------------------------------
$need = @()
if (-not (Test-Path $gfc)) { $need += $gfc }
if (-not $skip['gsh'] -and -not (Test-Path $gsh)) { $need += $gsh }
if (-not $skip['charm']) {
  if (-not (Test-Path $charmExe)) { $need += $charmExe }
  if (-not (Test-Path "$charmBin\libcharm.dll")) { $need += "$charmBin\libcharm.dll" }
}
if (-not $skip['shtools'] -and -not (Test-Path $shtoolsExe)) { $need += $shtoolsExe }
if ((-not $skip['charm'] -or -not $skip['shtools']) -and $Msys2Bin -and -not (Test-Path "$Msys2Bin\libgomp-1.dll")) {
  $need += "libgomp-1.dll in the MSYS2 bin folder ($Msys2Bin)"
}
if (-not (Test-Path $py)) { $need += $py }
if (-not $skip['graflab']) {
  if ($GrafLabDir -and -not (Test-Path "$GrafLabDir\GrafLab.m")) { $need += "$GrafLabDir\GrafLab.m" }
  if ($Matlab -and -not (Test-Path $Matlab)) { $need += $Matlab }
}
if ($need.Count -gt 0) {
  $msg = "Missing files:`n  " + ($need -join "`n  ")
  if ($DryRun) { Write-Host $msg -ForegroundColor Yellow } else { throw $msg }
}

# ---------------------------------------------------------------------------------------------
# 2) Raw folder: a final run starts from an empty timings file; OnlyGsh replaces the gsh rows only
# ---------------------------------------------------------------------------------------------
$baseline = @{}
$compGsh = Join-Path $csvDir 'gshtools\comp_value.csv'
if (Test-Path $compGsh) {
  foreach ($r in (Import-Csv $compGsh | Where-Object { $_.pair -eq 'gsh_vs_graflab' })) {
    $baseline["$($r.recurrence_a)|$($r.recurrence_b)|$($r.set)|$($r.quantity)"] = [double]$r.max_abs_d
  }
}
if ($OnlyGsh) {
  if (-not (Test-Path $timings)) { throw "-OnlyGsh needs the frozen raw outputs: $timings does not exist." }
  foreach ($set in $Sets) {
    if (-not (Get-ChildItem $raw -Filter "graflab_*_${set}.csv" -ErrorAction SilentlyContinue)) {
      throw "-OnlyGsh: no frozen GrafLab reference for set $set in $raw."
    }
  }
  if (-not $DryRun) {
    Copy-Item $timings "$raw\timings.before_onlygsh.txt" -Force
    $kept = Get-Content $timings | Where-Object { $_ -notmatch '^gsh,' }
    Set-Content -Path $timings -Value $kept -Encoding utf8
    Write-Host "OnlyGsh: gsh rows removed from timings.csv (copy kept as timings.before_onlygsh.txt)"
  }
} elseif ((Test-Path $timings) -and -not $Append) {
  if ($DryRun) { Write-Host "NOTE: $timings exists; a real run would stop here (use -Append or move it)." -ForegroundColor Yellow }
  else { throw "$timings already exists. Move it away or pass -Append (this run would add rows to it)." }
}
if (-not $DryRun) { New-Item -ItemType Directory -Force $raw | Out-Null }

# Points: the full sets, or the first -MaxPoints rows of each set for smoke tests.
function Get-Points([string]$set) {
  $src = "$data\points_$set.csv"
  if ($MaxPoints -le 0) { return Rel $src }
  $dir = "$raw\points_subset"
  $dst = "$dir\points_${set}_first$MaxPoints.csv"
  if (-not $DryRun) {
    New-Item -ItemType Directory -Force $dir | Out-Null
    Get-Content $src -TotalCount ($MaxPoints + 1) | Set-Content -Encoding utf8 $dst
  }
  return Rel $dst
}

# ---------------------------------------------------------------------------------------------
# 3) Helpers: run one tool, fail on non-zero exit; one progress line per run
# ---------------------------------------------------------------------------------------------
function Invoke-Tool([string]$exe, [string[]]$arguments, [string[]]$extraPath = @()) {
  if ($DryRun) { Write-Host "    DRY RUN: $exe $($arguments -join ' ')"; return }
  $old = $env:PATH
  if ($extraPath.Count -gt 0) { $env:PATH = (($extraPath | Where-Object { $_ }) -join ';') + ';' + $old }
  try {
    & $exe @arguments
    $code = $LASTEXITCODE
  } finally {
    $env:PATH = $old
  }
  if ($code -ne 0) { throw "exit code $code from: $exe $($arguments -join ' ')" }
}

function Run-Step([string]$label, [scriptblock]$body) {
  Step "START $label"
  $sw = [Diagnostics.Stopwatch]::StartNew()
  try {
    & $body
  } catch {
    Write-Host "[$(Stamp)] FAILED $label : $($_.Exception.Message)" -ForegroundColor Red
    throw
  }
  $sw.Stop()
  Step ("DONE  {0} in {1} | total elapsed {2}" -f $label, (Fmt $sw.Elapsed), (Fmt $total.Elapsed))
}

$rawRel = Rel $raw
$timingsRel = Rel $timings
$gfcRel = Rel $gfc
$msysPath = if ($Msys2Bin) { $Msys2Bin } else { '' }   # DLL search folders: absolute on purpose

Push-Location $root
if (-not $DryRun) { Start-Transcript -Path "$raw\run_all.log" -Append | Out-Null }
try {
  # -------------------------------------------------------------------------------------------
  # 4) The runs, set by set
  # -------------------------------------------------------------------------------------------
  foreach ($set in $Sets) {
    $pts = Get-Points $set
    Step "=== set ${set}: points $pts ==="

    if (-not $skip['gsh']) {
      foreach ($rec in 'fast', 'msfc', 'sfc', 'era') {
        $cfg = if ($rec -eq 'fast') { @($threads, 1) } else { @($threads) }
        foreach ($th in $cfg) {
          Run-Step "gsh $rec $set t$th" {
            Invoke-Tool (Rel $gsh) @('--gfc', $gfcRel, '--points', $pts, '--recurrence', $rec, '--nmax', $Nmax,
              '--out', "$rawRel\gsh_${rec}_${set}_t$th.csv", '--timings', $timingsRel, '--set', $set,
              '--threads', $th, '--repeat', $Repeat)
          }
        }
      }
    }

    if (-not $skip['graflab']) {
      foreach ($rec in 'sfc', 'msfc', 'era') {
        Run-Step "GrafLab $rec $set (MATLAB threads, repeat $Repeat)" {
          # MATLAB starts in its own folder: absolute paths here; run_graflab.m writes repo-relative names.
          $cmd = "addpath('$bench\code\matlab'); run_graflab('$gfc','$(Resolve-Path $pts)',$Nmax,'$rec','$set','$raw',$Repeat,'$GrafLabDir')"
          Invoke-Tool $Matlab @('-batch', $cmd)
        }
      }
    }

    # pyshtools: serial library, one thread (official form, no parallelism added)
    if (-not $skip['pyshtools']) {
      Run-Step "pyshtools $set t1" {
        Invoke-Tool (Rel $py) @("benchmark\shs\code\python\run_pyshtools.py", '--gfc', $gfcRel, '--points', $pts, '--nmax', $Nmax,
          '--out', "$rawRel\pyshtools_native_${set}_t1.csv", '--timings', $timingsRel, '--set', $set,
          '--repeat', $Repeat, '--threads', 1)
      }
    }

    if (-not $skip['pyharm']) {
      foreach ($th in $threads, 1) {
        Run-Step "pyharm $set t$th" {
          Invoke-Tool (Rel $py) @("benchmark\shs\code\python\run_pyharm.py", '--gfc', $gfcRel, '--points', $pts, '--nmax', $Nmax,
            '--out', "$rawRel\pyharm_native_${set}_t$th.csv", '--timings', $timingsRel, '--set', $set,
            '--repeat', $Repeat, '--threads', $th)
        }
      }
    }

    # SHTOOLS driver: serial loop (threads 1); its OpenMP is not used in the official runs.
    if (-not $skip['shtools']) {
      Run-Step "shtools driver $set t1" {
        Invoke-Tool (Rel $shtoolsExe) @('--gfc', $gfcRel, '--points', $pts, '--nmax', $Nmax,
          '--out', "$rawRel\shtools_native_${set}_t1.csv", '--timings', $timingsRel, '--set', $set,
          '--threads', 1, '--repeat', $Repeat) @($shtoolsBin, $msysPath)
      }
    }

    # CHarm: one batched library call over all points (--parallel charm), at 12 and 1 threads.
    if (-not $skip['charm']) {
      foreach ($th in $threads, 1) {
        Run-Step "charm --parallel charm $set t$th" {
          Invoke-Tool (Rel $charmExe) @('--gfc', $gfcRel, '--points', $pts, '--nmax', $Nmax,
            '--out', "$rawRel\charm_native_charm_${set}_t$th.csv", '--timings', $timingsRel, '--set', $set,
            '--threads', $th, '--repeat', $Repeat, '--parallel', 'charm') @($charmBin, $msysPath)
        }
      }
    }
  }

  # -------------------------------------------------------------------------------------------
  # 5) Consolidation, comparison, figures, regression check (OnlyGsh), publication check
  # -------------------------------------------------------------------------------------------
  Run-Step 'consolidate.py (results\csv)' {
    Invoke-Tool (Rel $py) @("benchmark\shs\code\python\consolidate.py", '--raw', $rawRel, '--csv', (Rel $csvDir))
  }
  Run-Step 'compare.py (comp_value.csv, summary_tables.md, figures)' {
    Invoke-Tool (Rel $py) @("benchmark\shs\code\python\compare.py", '--raw', $rawRel, '--csv', (Rel $csvDir), '--img', (Rel $img))
  }

  if ($OnlyGsh -and -not $DryRun) {
    Step 'Regression check: gsh against the frozen GrafLab, max abs(d) of V, T, dg'
    if ($baseline.Count -eq 0) {
      Write-Host 'No stored gsh comparison before this run (results\csv\gshtools\comp_value.csv): nothing to compare. Baseline recorded from this run.' -ForegroundColor Yellow
    } else {
      $failures = @()
      $rowsOut = @()
      foreach ($r in (Import-Csv $compGsh | Where-Object { $_.pair -eq 'gsh_vs_graflab' })) {
        $key = "$($r.recurrence_a)|$($r.recurrence_b)|$($r.set)|$($r.quantity)"
        $new = [double]$r.max_abs_d
        if (-not $baseline.ContainsKey($key)) { $rowsOut += "  NEW     $key  max|d| = $new"; continue }
        $old = $baseline[$key]
        $tol = [Math]::Max(1e-9 * [Math]::Abs($old), 1e-9)   # relative 1e-9, floor 1e-9 in the quantity unit (round-off of the fast kernel)
        $ok = [Math]::Abs($new - $old) -le $tol
        $status = if ($ok) { 'OK     ' } else { 'CHANGED' }
        $rowsOut += "  $status $key  before = $old  now = $new"
        if (-not $ok) { $failures += $key }
      }
      $rowsOut | ForEach-Object { Write-Host $_ }
      if ($failures.Count -gt 0) { throw "Regression check failed for: $($failures -join '; ')" }
      Write-Host 'Regression check passed (tolerance: relative 1e-9 or absolute 1e-9 in the quantity unit, whichever is larger).'
    }
  }

  $post = Find-PrivateText
  if ($post) {
    Write-Host 'Publication check FAILED: local paths or personal names in the outputs:' -ForegroundColor Red
    foreach ($m in ($post | Select-Object -First 20)) { Write-Host "  $($m.Path):$($m.LineNumber): $($m.Line.Trim())" }
    throw "Publication check: $($post.Count) line(s) with local paths or personal names. Fix the outputs before publishing."
  }
  Write-Host 'Publication check passed: no local paths or personal names in csv\, results\*.md or img\.'
  Step ("ALL DONE in {0}" -f (Fmt $total.Elapsed))
} finally {
  Pop-Location
  if (-not $DryRun) { Stop-Transcript | Out-Null }
}
