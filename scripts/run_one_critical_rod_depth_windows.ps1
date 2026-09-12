param(
    [Parameter(Mandatory = $true, Position = 0)]
    [double]$RodDepthCm,

    [int]$Particles = 500000,
    [int]$Inactive = 300,
    [int]$Active = 200,
    [string]$Configuration = "Release",
    [string]$ExePath = "",
    [string[]]$CaseArgs = @("--critical-search-binary")
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RootDir = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RootDir

if (-not $ExePath) {
    $ExePath = Join-Path $RootDir "x64\$Configuration\C5G7_3-DExtension.exe"
}

if (-not (Test-Path -LiteralPath $ExePath)) {
    throw "Executable not found: $ExePath. Build Release x64 in Visual Studio first, or pass -ExePath."
}

$culture = [System.Globalization.CultureInfo]::InvariantCulture
$depth = [string]::Format($culture, "{0:F6}", $RodDepthCm)

& $ExePath @CaseArgs "--rod-depth=$depth" "--particles=$Particles" "--inactive=$Inactive" "--active=$Active"
exit $LASTEXITCODE
