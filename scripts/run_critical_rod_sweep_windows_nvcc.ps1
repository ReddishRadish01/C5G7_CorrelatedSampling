param(
    [double]$StartCm = 32.0,
    [double]$EndCm = 37.0,
    [double]$CellHeightCm = 0.21,
    [int]$Particles = 500000,
    [int]$Inactive = 300,
    [int]$Active = 200,
    [string]$CudaArch = "sm_120",
    [string]$NvccPath = "nvcc",
    [string]$ExePath = "build\windows_nvcc\c5g7.exe",
    [string]$LogDir = "results\critical_rod_sweep_windows_nvcc",
    [string[]]$CaseArgs = @("--critical-search-binary"),
    [switch]$NoBuild,
    [switch]$NoVsDevEnv,
    [switch]$AllowUnsupportedCompiler
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
if ($PSVersionTable.PSVersion.Major -ge 7) {
    $PSNativeCommandUseErrorActionPreference = $false
}

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
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

$summary = Join-Path $LogDir "summary.csv"
"depth_cm,actual_depth_cm,k_mean,k_stderr_mean,status,log_file" | Set-Content -LiteralPath $summary

$culture = [System.Globalization.CultureInfo]::InvariantCulture
$startIdx = [math]::Ceiling(($StartCm - 1.0e-9) / $CellHeightCm)
$endIdx = [math]::Floor(($EndCm + 1.0e-9) / $CellHeightCm)

if ($startIdx -gt $endIdx) {
    throw "No quantized rod depths found in [$StartCm, $EndCm] with cell height $CellHeightCm cm."
}

$firstDepth = [string]::Format($culture, "{0:F6}", $startIdx * $CellHeightCm)
$lastDepth = [string]::Format($culture, "{0:F6}", $endIdx * $CellHeightCm)

Write-Host "Running direct-nvcc CriticalSearch_Binary rod-depth sweep"
Write-Host "  exe: $exe"
Write-Host "  actual depth range: $firstDepth to $lastDepth cm"
Write-Host "  cell height: $CellHeightCm cm"
Write-Host "  particles/inactive/active: $Particles/$Inactive/$Active"
Write-Host "  logs: $LogDir"

for ($idx = $startIdx; $idx -le $endIdx; ++$idx) {
    $depth = [string]::Format($culture, "{0:F6}", $idx * $CellHeightCm)
    $depthTag = $depth.Replace(".", "p")
    $logFile = Join-Path $LogDir "rod_${depthTag}cm.log"

    Write-Host "==> rod depth $depth cm"

    $programArgs = @()
    $programArgs += $CaseArgs
    $programArgs += "--rod-depth=$depth"
    $programArgs += "--particles=$Particles"
    $programArgs += "--inactive=$Inactive"
    $programArgs += "--active=$Active"

    $statusCode = 0
    try {
        & $exe @programArgs > $logFile 2>&1
        $statusCode = $LASTEXITCODE
        if ($null -eq $statusCode) {
            $statusCode = 0
        }
    }
    catch {
        $_ | Out-File -LiteralPath $logFile -Append
        $statusCode = 1
    }

    $actual = ""
    $kMean = ""
    $kStderr = ""
    if (Test-Path -LiteralPath $logFile) {
        $lines = Get-Content -LiteralPath $logFile
        $actualLine = $lines | Select-String -Pattern "Actual (?:modeled|quantized) rod depth:" | Select-Object -Last 1
        $kMeanLine = $lines | Select-String -Pattern "^k_mean:" | Select-Object -Last 1
        $kStderrLine = $lines | Select-String -Pattern "^k_stderr\(mean\):" | Select-Object -Last 1

        if ($actualLine) { $actual = ($actualLine.Line -replace ".*Actual (?:modeled|quantized) rod depth:\s*", "" -replace "\s*cm.*", "").Trim() }
        if ($kMeanLine) { $kMean = ($kMeanLine.Line -replace "^k_mean:\s*", "").Trim() }
        if ($kStderrLine) { $kStderr = ($kStderrLine.Line -replace "^k_stderr\(mean\):\s*", "").Trim() }
    }

    if ($statusCode -eq 0) {
        "$depth,$actual,$kMean,$kStderr,OK,$logFile" | Add-Content -LiteralPath $summary
    }
    else {
        Write-Warning "Depth $depth failed with exit code $statusCode; see $logFile"
        "$depth,$actual,$kMean,$kStderr,FAILED:$statusCode,$logFile" | Add-Content -LiteralPath $summary
    }
}

Write-Host "Done. Summary written to $summary"
