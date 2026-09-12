[CmdletBinding()]
param(
    [string]$OutputPath = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot

if (-not (Get-Command tar -ErrorAction SilentlyContinue)) {
    throw "Required command 'tar' was not found in PATH."
}

if ([string]::IsNullOrWhiteSpace($OutputPath)) {
    $OutputPath = Join-Path $projectRoot "dist\C5G7_CS_Newton_SLURM_current.tar.gz"
}
elseif (-not [System.IO.Path]::IsPathRooted($OutputPath)) {
    $OutputPath = Join-Path $projectRoot $OutputPath
}

$OutputPath = [System.IO.Path]::GetFullPath($OutputPath)
$outputDirectory = Split-Path -Parent $OutputPath
New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null

$requiredFiles = @(
    "Makefile",
    "main.cu",
    "Core.cu",
    "Cycle.cu",
    "Neutron.cu",
    "RNG.cu",
    "Tally.cu",
    "XSParser.cu",
    "C5txt/CR.txt",
    "Geometry_ExtensionRod/C5G7CoreGeometry.txt",
    "scripts/critical_rod_newton_cs.env",
    "scripts/run_critical_rod_newton_cs.sh",
    "scripts/run_critical_rod_newton_cs_slurm.sbatch"
)

foreach ($relativePath in $requiredFiles) {
    $localPath = Join-Path $projectRoot ($relativePath -replace "/", "\")
    if (-not (Test-Path -LiteralPath $localPath -PathType Leaf)) {
        throw "Required project file is missing: $relativePath"
    }
}

Push-Location $projectRoot
try {
    & tar -czf $OutputPath `
        --exclude="./.git" `
        --exclude="./.vs" `
        --exclude="./build" `
        --exclude="./x64" `
        --exclude="./results" `
        --exclude="./resultPlot" `
        --exclude="./kTrend" `
        --exclude="./tmp" `
        --exclude="./dist" `
        --exclude="./.codex-build" `
        --exclude="./C5G7_Cor.9b61fcbf" `
        --exclude="./*.log" `
        --exclude="./newton_history.csv" `
        --exclude="./cs_history_*.txt" `
        .
    if ($LASTEXITCODE -ne 0) {
        throw "tar creation failed with exit code $LASTEXITCODE."
    }
}
finally {
    Pop-Location
}

$archiveEntries = @(& tar -tzf $OutputPath)
if ($LASTEXITCODE -ne 0) {
    throw "tar verification failed with exit code $LASTEXITCODE."
}
$normalizedEntries = @($archiveEntries | ForEach-Object { $_ -replace '^\./', '' })

$missingEntries = @($requiredFiles | Where-Object { $normalizedEntries -notcontains $_ })
if ($missingEntries.Count -gt 0) {
    throw "Archive verification failed; missing: $($missingEntries -join ', ')"
}

$archive = Get-Item -LiteralPath $OutputPath
$hash = Get-FileHash -LiteralPath $OutputPath -Algorithm SHA256

Write-Host "SLURM transfer package is ready."
Write-Host "Path: $($archive.FullName)"
Write-Host ("Size: {0:N2} MiB" -f ($archive.Length / 1MB))
Write-Host "SHA256: $($hash.Hash)"
Write-Host "Verified required files: $($requiredFiles.Count)"
