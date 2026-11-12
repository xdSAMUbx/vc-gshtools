<#
.SYNOPSIS
  Builds vc-gshtools from scratch on Windows 11 x64: checks the tools, bootstraps vcpkg if needed,
  configures, builds and runs the tests.

.DESCRIPTION
  Steps, in order:
    1. Checks the platform (Windows, x64) and that the script runs inside the repository.
    2. Finds Visual Studio C++ tools (vswhere), git, CMake, Ninja and, for the clang preset, clang-cl.
       Missing tools are reported with the winget command that installs them.
    3. With -InstallPrereqs only: installs the missing tools with winget (asks first unless -Yes).
       winget may ask for administrator rights (UAC); this script never elevates itself.
    4. Loads the compiler environment of Visual Studio into this process only (Launch-VsDevShell.ps1).
    5. Finds vcpkg at -VcpkgRoot, $env:VCPKG_ROOT or C:\vcpkg. If it is missing, clones
       https://github.com/microsoft/vcpkg there and runs bootstrap-vcpkg.bat (announced before it runs).
    6. cmake --preset <Preset>, cmake --build --preset <Preset>, ctest --preset <Preset> (unless -SkipTests).
    7. Prints where the binaries are and the next step.

  Without -InstallPrereqs the script installs nothing. Without -DryRun it changes only the build folder
  build\<Preset> and, if vcpkg is missing, the vcpkg folder. Re-running it is safe.

  Exit codes: 0 success, 1 failure, 2 missing prerequisites (nothing was installed).

.PARAMETER Preset
  msvc (cl.exe) or clang (clang-cl.exe). Default: clang.

.PARAMETER VcpkgRoot
  vcpkg folder. Default: $env:VCPKG_ROOT if set, otherwise C:\vcpkg.

.PARAMETER InstallPrereqs
  Allows winget to install the missing tools (Git, CMake, Ninja, LLVM for the clang preset, and the
  Visual Studio Build Tools C++ workload). Off by default.

.PARAMETER Yes
  With -InstallPrereqs, does not ask for confirmation.

.PARAMETER SkipTests
  Builds but does not run ctest.

.PARAMETER Clean
  Deletes build\<Preset> before configuring.

.PARAMETER DryRun
  Prints what would be done. Changes nothing (no installs, no clone, no build folder, no ctest).

.EXAMPLE
  .\config\config.ps1 -DryRun
.EXAMPLE
  .\config\config.ps1 -Preset clang -InstallPrereqs
#>
[CmdletBinding()]
param(
    [ValidateSet('msvc', 'clang')][string]$Preset = 'clang',
    [string]$VcpkgRoot,
    [switch]$InstallPrereqs,
    [switch]$Yes,
    [switch]$SkipTests,
    [switch]$Clean,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
$BuildDir = Join-Path (Join-Path $RepoRoot 'build') $Preset
$VcpkgUrl = 'https://github.com/microsoft/vcpkg.git'
$ExitMissing = 2

function Write-Step([string]$Text) { Write-Host ''; Write-Host "==> $Text" -ForegroundColor Cyan }
function Write-Note([string]$Text) { Write-Host "    $Text" }
function Write-Plan([string]$Text) { Write-Host "    [dry run] $Text" -ForegroundColor Yellow }

function Find-Tool([string]$Name, [string[]]$Fallbacks) {
    $cmd = Get-Command -Name $Name -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($cmd) { return $cmd.Source }
    foreach ($f in $Fallbacks) {
        if ($f -and (Test-Path -LiteralPath $f -PathType Leaf)) { return $f }
    }
    return $null
}

function Invoke-Checked([string]$File, [string[]]$Arguments) {
    if ($DryRun) { Write-Plan ("would run: {0} {1}" -f $File, ($Arguments -join ' ')); return }
    & $File @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Command failed (exit $LASTEXITCODE): $File $($Arguments -join ' ')" }
}

function Get-VsInstallPath {
    $pf86 = ${env:ProgramFiles(x86)}
    if (-not $pf86) { return $null }
    $vswhere = Join-Path $pf86 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere -PathType Leaf)) { return $null }
    $path = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
        -property installationPath 2>$null | Select-Object -First 1
    if ($path -and (Test-Path -LiteralPath $path)) { return [string]$path }
    return $null
}

function Get-Prerequisites([string]$VsPath) {
    $pf = $env:ProgramFiles
    $vsCMake = if ($VsPath) { Join-Path $VsPath 'Common7\IDE\CommonExtensions\Microsoft\CMake' } else { $null }
    $list = @()
    $list += [pscustomobject]@{ Name = 'Visual Studio C++ tools'; Found = [bool]$VsPath; WingetId = $null }
    $list += [pscustomobject]@{ Name = 'git'; Found = [bool](Find-Tool 'git' @((Join-Path $pf 'Git\cmd\git.exe'))); WingetId = 'Git.Git' }
    $cmake = Find-Tool 'cmake' @((Join-Path $pf 'CMake\bin\cmake.exe'), $(if ($vsCMake) { Join-Path $vsCMake 'CMake\bin\cmake.exe' }))
    $list += [pscustomobject]@{ Name = 'cmake'; Found = [bool]$cmake; WingetId = 'Kitware.CMake' }
    $ninja = Find-Tool 'ninja' @($(if ($vsCMake) { Join-Path $vsCMake 'Ninja\ninja.exe' }))
    $list += [pscustomobject]@{ Name = 'ninja'; Found = [bool]$ninja; WingetId = 'Ninja-build.Ninja' }
    if ($Preset -eq 'clang') {
        $clang = Find-Tool 'clang-cl' @((Join-Path $pf 'LLVM\bin\clang-cl.exe'))
        $list += [pscustomobject]@{ Name = 'clang-cl (LLVM)'; Found = [bool]$clang; WingetId = 'LLVM.LLVM' }
    }
    return $list
}

# ---------------------------------------------------------------- 1. platform and repository
Write-Step '1/7 Platform and repository'
if ($env:OS -ne 'Windows_NT') { throw 'This script supports Windows only.' }
$isX64 = ($env:PROCESSOR_ARCHITECTURE -eq 'AMD64') -or ($env:PROCESSOR_ARCHITEW6432 -eq 'AMD64')
if (-not $isX64) { throw 'This script supports Windows x64 only.' }
Write-Note ("Windows: {0}; PowerShell {1} ({2})" -f [Environment]::OSVersion.VersionString, $PSVersionTable.PSVersion, $PSVersionTable.PSEdition)
Write-Note "Repository: $RepoRoot"
Write-Note "Preset: $Preset (compiler: $(if ($Preset -eq 'msvc') { 'cl.exe' } else { 'clang-cl.exe' }))"
if ($DryRun) { Write-Note 'DRY RUN: nothing will be installed, cloned, built or changed.' }
if (-not (Test-Path -LiteralPath (Join-Path $RepoRoot 'CMakePresets.json') -PathType Leaf)) {
    throw "CMakePresets.json not found in $RepoRoot. Run this script from a vc-gshtools checkout."
}

# ---------------------------------------------------------------- 2. tools
Write-Step '2/7 Tools (Visual Studio C++, git, CMake, Ninja)'
$vsPath = Get-VsInstallPath
$prereqs = Get-Prerequisites $vsPath
foreach ($p in $prereqs) {
    $state = if ($p.Found) { 'found' } else { 'MISSING' }
    Write-Note ("{0,-26} {1}" -f $p.Name, $state)
}
$missing = @($prereqs | Where-Object { -not $_.Found })

if ($missing.Count -gt 0) {
    Write-Host ''
    Write-Note 'Install commands (run them yourself, or rerun this script with -InstallPrereqs):'
    foreach ($m in $missing) {
        if ($m.WingetId) {
            Write-Note ("  winget install --id {0} -e --accept-package-agreements --accept-source-agreements" -f $m.WingetId)
        } else {
            Write-Note '  winget install --id Microsoft.VisualStudio.BuildTools -e --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"'
        }
    }
    if (-not $InstallPrereqs) {
        if ($DryRun) {
            Write-Plan 'a real run would stop here with exit code 2 (nothing installed); rerun with -InstallPrereqs to let winget install the tools'
        } else {
            Write-Host ''
            Write-Host 'Missing prerequisites. Nothing was installed. Use -InstallPrereqs to let winget install them.' -ForegroundColor Red
            exit $ExitMissing
        }
    }
    if ($InstallPrereqs -and -not (Find-Tool 'winget' @())) {
        throw 'winget was not found. Install the missing tools manually (see the commands above) and rerun.'
    }
    if ($InstallPrereqs) {
    Write-Host ''
    Write-Note 'Installing with winget. winget may ask for administrator rights (UAC); this script does not elevate itself.'
    Write-Note 'The Visual Studio Build Tools download is large (several GB) and can take a long time.'
    if (-not $DryRun -and -not $Yes) {
        if (-not [Environment]::UserInteractive -or [Console]::IsInputRedirected) {
            throw 'Confirmation needed: rerun with -Yes, or run interactively.'
        }
        $answer = Read-Host 'Install the missing tools now? [y/N]'
        if ($answer -notmatch '^(y|yes)$') { Write-Note 'Cancelled. Nothing was installed.'; exit $ExitMissing }
    }
    foreach ($m in $missing) {
        if ($m.WingetId) {
            Invoke-Checked 'winget' @('install', '--id', $m.WingetId, '-e', '--accept-package-agreements', '--accept-source-agreements')
        } else {
            Invoke-Checked 'winget' @('install', '--id', 'Microsoft.VisualStudio.BuildTools', '-e', '--accept-package-agreements',
                '--accept-source-agreements', '--override', '--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended')
        }
    }
    if (-not $DryRun) {
        # Re-read PATH so the tools installed just now are visible in this session.
        $env:Path = [Environment]::GetEnvironmentVariable('Path', 'Machine') + ';' + [Environment]::GetEnvironmentVariable('Path', 'User')
        $vsPath = Get-VsInstallPath
        $prereqs = Get-Prerequisites $vsPath
        $still = @($prereqs | Where-Object { -not $_.Found })
        if ($still.Count -gt 0) {
            throw ("Still missing after installation: {0}. Open a new terminal and rerun this script." -f (($still | ForEach-Object { $_.Name }) -join ', '))
        }
    }
    }
}

# ---------------------------------------------------------------- 3. compiler environment
Write-Step '3/7 Compiler environment'
if ($DryRun) {
    Write-Plan "would load $vsPath\Common7\Tools\Launch-VsDevShell.ps1 into this process (-Arch amd64 -HostArch amd64)"
    if ($Preset -eq 'clang') { Write-Plan 'would add the LLVM bin folder to PATH for this process' }
} else {
    $launch = Join-Path $vsPath 'Common7\Tools\Launch-VsDevShell.ps1'
    if (-not (Test-Path -LiteralPath $launch -PathType Leaf)) { throw "Launch-VsDevShell.ps1 not found under $vsPath." }
    & $launch -Arch amd64 -HostArch amd64 -SkipAutomaticLocation | Out-Null
    Write-Note 'Visual Studio environment loaded for this process only.'
    if ($Preset -eq 'clang') {
        $llvmBin = Join-Path $env:ProgramFiles 'LLVM\bin'
        if (-not (Get-Command clang-cl -ErrorAction SilentlyContinue) -and (Test-Path -LiteralPath $llvmBin)) {
            $env:Path = $llvmBin + ';' + $env:Path
        }
        if (-not (Get-Command clang-cl -ErrorAction SilentlyContinue)) { throw 'clang-cl not found on PATH.' }
    }
    $cl = Get-Command cl.exe -ErrorAction SilentlyContinue
    if ($cl) { Write-Note "cl.exe: $($cl.Source)" }
}

# ---------------------------------------------------------------- 4. vcpkg
Write-Step '4/7 vcpkg'
$vcpkgRootResolved = if ($VcpkgRoot) { $VcpkgRoot } elseif ($env:VCPKG_ROOT) { $env:VCPKG_ROOT } else { 'C:\vcpkg' }
$vcpkgRootResolved = [IO.Path]::GetFullPath($vcpkgRootResolved).TrimEnd('\')
$toolchain = Join-Path $vcpkgRootResolved 'scripts\buildsystems\vcpkg.cmake'
$bootstrap = Join-Path $vcpkgRootResolved 'bootstrap-vcpkg.bat'
Write-Note "vcpkg root: $vcpkgRootResolved"
if (-not (Test-Path -LiteralPath $toolchain -PathType Leaf)) {
    if ((Test-Path -LiteralPath $vcpkgRootResolved) -and @(Get-ChildItem -LiteralPath $vcpkgRootResolved -Force).Count -gt 0) {
        throw "vcpkg root '$vcpkgRootResolved' exists, is not empty, and is not a vcpkg checkout. Use -VcpkgRoot to point elsewhere."
    }
    Write-Note "vcpkg is not installed there. It will be cloned from $VcpkgUrl (full clone: the manifest pins a baseline commit) and bootstrapped."
    Write-Note 'This downloads from GitHub.'
    if (-not $DryRun) {
        Invoke-Checked 'git' @('clone', $VcpkgUrl, $vcpkgRootResolved)
    } else {
        Write-Plan "would run: git clone $VcpkgUrl $vcpkgRootResolved"
    }
}
if (-not (Test-Path -LiteralPath (Join-Path $vcpkgRootResolved 'vcpkg.exe') -PathType Leaf)) {
    Write-Note 'Bootstrapping vcpkg (bootstrap-vcpkg.bat -disableMetrics).'
    if (-not $DryRun) {
        if (-not (Test-Path -LiteralPath $bootstrap -PathType Leaf)) { throw "bootstrap-vcpkg.bat not found in $vcpkgRootResolved." }
        & cmd.exe /c "`"$bootstrap`" -disableMetrics"
        if ($LASTEXITCODE -ne 0) { throw "bootstrap-vcpkg.bat failed (exit $LASTEXITCODE)." }
    } else {
        Write-Plan 'would run: bootstrap-vcpkg.bat -disableMetrics'
    }
} else {
    Write-Note 'vcpkg is already bootstrapped.'
}
$env:VCPKG_ROOT = $vcpkgRootResolved
Write-Note "VCPKG_ROOT set for this process: $vcpkgRootResolved"

# ---------------------------------------------------------------- 5. configure
Write-Step "5/7 Configure (cmake --preset $Preset)"
Write-Note 'The first configure builds the dependencies (FFTW, Eigen, GoogleTest) with vcpkg. This takes a long time and needs internet access.'
if ($Clean -and (Test-Path -LiteralPath $BuildDir)) {
    $buildRoot = Join-Path $RepoRoot 'build'
    if (-not $BuildDir.StartsWith($buildRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Refusing to delete $BuildDir." }
    if ($DryRun) { Write-Plan "would delete $BuildDir" } else { Remove-Item -LiteralPath $BuildDir -Recurse -Force }
}
Push-Location $RepoRoot
try {
    Invoke-Checked 'cmake' @('--preset', $Preset)

    Write-Step "6/7 Build (cmake --build --preset $Preset)"
    Invoke-Checked 'cmake' @('--build', '--preset', $Preset)

    if ($SkipTests) {
        Write-Step '7/7 Tests skipped (-SkipTests)'
    } else {
        Write-Step "7/7 Tests (ctest --preset $Preset)"
        Invoke-Checked 'ctest' @('--preset', $Preset)
    }
} finally {
    Pop-Location
}

# ---------------------------------------------------------------- summary
Write-Host ''
if ($DryRun) {
    Write-Host 'Dry run finished. Nothing was changed.' -ForegroundColor Green
    exit 0
}
$binDir = Join-Path $BuildDir 'bin'
Write-Host 'Build finished.' -ForegroundColor Green
Write-Note "Binaries: $binDir"
if (Test-Path -LiteralPath $binDir) {
    Get-ChildItem -LiteralPath $binDir -Filter '*.exe' | ForEach-Object { Write-Note ("  {0}" -f $_.Name) }
}
Write-Host ''
Write-Note 'Next step: the SHS benchmark, see benchmark\shs\results\Summary.md (section "How each engine was run")'
Write-Note 'and benchmark\shs\data\DATA.md. Download the model from the Release data-v1, then decompress it:'
Write-Note '  xz -dk benchmark\shs\data\XGM2019e_2190_tide_free.gfc.xz   (needs xz; or any xz-compatible tool)'
exit 0
