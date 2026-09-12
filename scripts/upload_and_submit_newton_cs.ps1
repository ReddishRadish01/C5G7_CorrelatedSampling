[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[A-Za-z0-9._@-]+$')]
    [string]$Server,

    [ValidatePattern('^[A-Za-z0-9._/~\-]+$')]
    [string]$RemoteDir = "C5G7_runs",

    [ValidateRange(1, 65535)]
    [int]$Port = 22,

    [ValidatePattern('^[A-Za-z0-9._-]*$')]
    [string]$Partition = "",

    [ValidatePattern('^[A-Za-z0-9._-]*$')]
    [string]$Account = "",

    [ValidatePattern('^[A-Za-z0-9._/+\-]*$')]
    [string]$CudaModule = "",

    [double]$InitialDepthCm = 64.24,
    [double]$LowerCm = 21.42,
    [double]$UpperCm = 128.52,
    [double]$ProbeStepCm = 0.25,
    [double]$MinProbeStepCm = 0.03125,
    [int]$MaxProbeRetries = 3,
    [double]$Damping = 1.0,
    [double]$MaxNewtonStepCm = 30.0,
    [double]$KTolerancePcm = 10.0,
    [double]$MaxKStderrPcm = 10.0,
    [double]$ConfidenceSigma = 2.0,
    [double]$SourceTailAbsPcm = 1.0,
    [double]$SourceTailRelFraction = 0.02,
    [double]$MaxSlopeRelativeDifference = 0.05,
    [int]$Particles = 2000000,
    [int]$Inactive = 100,
    [int]$Active = 200,
    [int]$MaxIterations = 8,
    [int]$CsGenerations = 10,
    [UInt64]$BaseSeed = 92235922383,

    [switch]$UploadOnly
)

$ErrorActionPreference = "Stop"
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$projectRoot = Split-Path -Parent $PSScriptRoot
$runName = "C5G7_" + (Get-Date -Format "yyyyMMdd_HHmmss")
$remoteRunDir = "$($RemoteDir.TrimEnd('/'))/$runName"
$archivePath = Join-Path ([System.IO.Path]::GetTempPath()) ("c5g7_" + [guid]::NewGuid().ToString("N") + ".tar.gz")
$remoteArchive = "$remoteRunDir/source.tar.gz"

foreach ($commandName in @("ssh", "scp", "tar")) {
    if (-not (Get-Command $commandName -ErrorAction SilentlyContinue)) {
        throw "Required command '$commandName' was not found in PATH."
    }
}

if ($LowerCm -ge $UpperCm) { throw "LowerCm must be less than UpperCm." }
if ($InitialDepthCm -lt $LowerCm -or $InitialDepthCm -gt $UpperCm) {
    throw "InitialDepthCm must be inside [LowerCm, UpperCm]."
}
if ($ProbeStepCm -le 0 -or $MinProbeStepCm -le 0 -or $MinProbeStepCm -gt $ProbeStepCm) {
    throw "Require 0 < MinProbeStepCm <= ProbeStepCm."
}
if ($MaxProbeRetries -lt 0) { throw "MaxProbeRetries cannot be negative." }
if ($Damping -le 0 -or $Damping -gt 1) { throw "Damping must be in (0, 1]." }
if ($MaxNewtonStepCm -le 0) { throw "MaxNewtonStepCm must be positive." }
if ($KTolerancePcm -lt 0 -or $MaxKStderrPcm -lt 0 -or $ConfidenceSigma -le 0 -or
    $SourceTailAbsPcm -lt 0 -or $SourceTailRelFraction -lt 0 -or
    $MaxSlopeRelativeDifference -lt 0 -or
    ($SourceTailAbsPcm -eq 0 -and $SourceTailRelFraction -eq 0)) {
    throw "Invalid statistical or source-tail tolerance."
}
if ($Particles -le 0 -or $Inactive -lt 0 -or $Active -lt 2 -or $MaxIterations -le 0) {
    throw "Invalid particle, cycle, or iteration count."
}
if ($CsGenerations -lt 1 -or $CsGenerations -gt 10) {
    throw "CsGenerations must be between 1 and 10 for this build."
}
if ($BaseSeed -gt [UInt64]281474976710655) {
    throw "BaseSeed must fit in the RNG's 48-bit state (0..281474976710655)."
}

function Format-Invariant([double]$Value) {
    return $Value.ToString("R", $culture)
}

try {
    Push-Location $projectRoot
    try {
        & tar -czf $archivePath `
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
        if ($LASTEXITCODE -ne 0) { throw "tar failed with exit code $LASTEXITCODE." }
    }
    finally {
        Pop-Location
    }

    Write-Host "Creating remote run directory: $Server`:$remoteRunDir"
    & ssh -p $Port $Server "mkdir -p $remoteRunDir"
    if ($LASTEXITCODE -ne 0) { throw "Remote mkdir failed with exit code $LASTEXITCODE." }

    Write-Host "Uploading source archive..."
    & scp -P $Port $archivePath "${Server}:$remoteArchive"
    if ($LASTEXITCODE -ne 0) { throw "scp failed with exit code $LASTEXITCODE." }

    & ssh -p $Port $Server "cd $remoteRunDir && tar -xzf source.tar.gz && chmod +x scripts/*.sh && mkdir -p results"
    if ($LASTEXITCODE -ne 0) { throw "Remote extraction failed with exit code $LASTEXITCODE." }

    if ($UploadOnly) {
        Write-Host "Upload complete: $Server`:$remoteRunDir"
        return
    }

    $sbatchOptions = @("--parsable")
    if ($Partition) { $sbatchOptions += "--partition=$Partition" }
    if ($Account) { $sbatchOptions += "--account=$Account" }

    $settings = @(
        "INITIAL_DEPTH_CM=$(Format-Invariant $InitialDepthCm)",
        "LOWER_CM=$(Format-Invariant $LowerCm)",
        "UPPER_CM=$(Format-Invariant $UpperCm)",
        "PROBE_STEP_CM=$(Format-Invariant $ProbeStepCm)",
        "MIN_PROBE_STEP_CM=$(Format-Invariant $MinProbeStepCm)",
        "MAX_PROBE_RETRIES=$MaxProbeRetries",
        "DAMPING=$(Format-Invariant $Damping)",
        "MAX_NEWTON_STEP_CM=$(Format-Invariant $MaxNewtonStepCm)",
        "K_TOL_PCM=$(Format-Invariant $KTolerancePcm)",
        "MAX_K_STDERR_PCM=$(Format-Invariant $MaxKStderrPcm)",
        "CONFIDENCE_SIGMA=$(Format-Invariant $ConfidenceSigma)",
        "SOURCE_TAIL_ABS_PCM=$(Format-Invariant $SourceTailAbsPcm)",
        "SOURCE_TAIL_REL_FRACTION=$(Format-Invariant $SourceTailRelFraction)",
        "MAX_SLOPE_RELATIVE_DIFFERENCE=$(Format-Invariant $MaxSlopeRelativeDifference)",
        "PARTICLES=$Particles",
        "INACTIVE=$Inactive",
        "ACTIVE=$Active",
        "MAX_ITER=$MaxIterations",
        "CS_GENERATIONS=$CsGenerations",
        "BASE_SEED=$BaseSeed"
    )
    if ($CudaModule) { $settings += "CUDA_MODULE=$CudaModule" }

    $submitCommand = @(
        "cd", $remoteRunDir, "&&", "sbatch"
    ) + $sbatchOptions + @(
        "scripts/run_critical_rod_newton_cs_slurm.sbatch"
    ) + $settings
    $remoteCommand = $submitCommand -join " "

    Write-Host "Submitting Slurm job..."
    $jobOutput = & ssh -p $Port $Server $remoteCommand
    if ($LASTEXITCODE -ne 0) { throw "sbatch failed with exit code $LASTEXITCODE." }
    $jobId = ([string]($jobOutput | Select-Object -Last 1)).Trim()
    if ($jobId -notmatch '^[0-9]+(?:;[A-Za-z0-9._-]+)?$') {
        throw "Could not parse the Slurm job id from: $($jobOutput -join ' | ')"
    }
    $jobNumber = $jobId.Split(';')[0]

    Write-Host "Submitted job: $jobId"
    Write-Host "Remote directory: $Server`:$remoteRunDir"
    Write-Host "Monitor: ssh -p $Port $Server squeue -j $jobNumber"
    Write-Host "Slurm log: $remoteRunDir/results/slurm-c5g7-newton-cs-$jobNumber.out"
}
finally {
    if (Test-Path -LiteralPath $archivePath) {
        Remove-Item -LiteralPath $archivePath -Force
    }
}
