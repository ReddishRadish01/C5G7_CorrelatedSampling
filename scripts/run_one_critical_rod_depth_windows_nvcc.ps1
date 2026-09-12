param(
    [Parameter(Mandatory = $true, Position = 0)]
    [double]$RodDepthCm,

    [int]$Particles = 500000,
    [int]$Inactive = 300,
    [int]$Active = 200,
    [string]$CudaArch = "sm_120",
    [string]$NvccPath = "nvcc",
    [string]$ExePath = "build\windows_nvcc\c5g7.exe",
    [string[]]$CaseArgs = @("--critical-search-binary"),
    [switch]$NoBuild,
    [switch]$NoVsDevEnv,
    [switch]$AllowUnsupportedCompiler,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ExtraProgramArgs
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RootDir = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RootDir

if (-not $NoBuild) {
    $buildArgs = @(
        "-NvccPath", $NvccPath,
        "-CudaArch", $CudaArch,
        "-Output", $ExePath
    )
    if ($NoVsDevEnv) { $buildArgs += "-NoVsDevEnv" }
    if ($AllowUnsupportedCompiler) { $buildArgs += "-AllowUnsupportedCompiler" }

    & (Join-Path $PSScriptRoot "build_windows_nvcc.ps1") @buildArgs
}

$exe = Resolve-Path -LiteralPath (Join-Path $RootDir $ExePath)
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$depth = [string]::Format($culture, "{0:F6}", $RodDepthCm)

$programArgs = @()
$programArgs += $CaseArgs
$programArgs += "--rod-depth=$depth"
$programArgs += "--particles=$Particles"
$programArgs += "--inactive=$Inactive"
$programArgs += "--active=$Active"
$programArgs += $ExtraProgramArgs

& $exe @programArgs
exit $LASTEXITCODE
