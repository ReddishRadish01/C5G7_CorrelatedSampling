[CmdletBinding()]
param(
    [string]$ConfigPath = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrWhiteSpace($ConfigPath)) {
    $ConfigPath = Join-Path $scriptDir 'critical_rod_newton_cs_windows_config.psd1'
}
$culture = [System.Globalization.CultureInfo]::InvariantCulture
$numberPattern = '[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?'
$rngModulus = [UInt64]281474976710656
$rngMask = [UInt64]281474976710655

$rootDir = (Resolve-Path (Join-Path $scriptDir '..')).Path
Set-Location $rootDir

$resolvedConfig = (Resolve-Path -LiteralPath $ConfigPath).Path
$config = Import-PowerShellDataFile -LiteralPath $resolvedConfig
$requiredKeys = @(
    'Particles', 'Inactive', 'Active', 'BaseSeed',
    'InitialDepthCm', 'LowerCm', 'UpperCm', 'TargetK',
    'ProbeStepCm', 'MinProbeStepCm', 'MaxProbeRetries', 'CsGenerations',
    'Damping', 'MaxNewtonStepCm', 'MaxIterations',
    'DepthToleranceCm', 'KTolerancePcm', 'MaxKStderrPcm',
    'MinAbsSlope', 'MinSlopeSigma', 'ConfidenceSigma',
    'SourceTailAbsPcm', 'SourceTailRelFraction', 'MaxSlopeRelativeDifference',
    'RequireSourceConvergence', 'FinalDirectCheck',
    'ExePath', 'NvccPath', 'CudaArch', 'ThreadsPerBlock',
    'DeviceLto', 'FastMath', 'ForceBuild', 'NoBuild'
)
foreach ($key in $requiredKeys) {
    if (-not $config.ContainsKey($key)) {
        throw "Missing '$key' in configuration file: $resolvedConfig"
    }
}

$particles = [int]$config.Particles
$inactive = [int]$config.Inactive
$active = [int]$config.Active
$baseSeed = [UInt64]$config.BaseSeed
$initialDepth = [double]$config.InitialDepthCm
$lower = [double]$config.LowerCm
$upper = [double]$config.UpperCm
$targetK = [double]$config.TargetK
$probeStepDefault = [double]$config.ProbeStepCm
$minProbeStep = [double]$config.MinProbeStepCm
$maxProbeRetries = [int]$config.MaxProbeRetries
$csGenerations = [int]$config.CsGenerations
$damping = [double]$config.Damping
$maxNewtonStep = [double]$config.MaxNewtonStepCm
$maxIterations = [int]$config.MaxIterations
$depthTolerance = [double]$config.DepthToleranceCm
$kTolerancePcm = [double]$config.KTolerancePcm
$maxKStderrPcm = [double]$config.MaxKStderrPcm
$minAbsSlope = [double]$config.MinAbsSlope
$minSlopeSigma = [double]$config.MinSlopeSigma
$confidenceSigma = [double]$config.ConfidenceSigma
$sourceTailAbsPcm = [double]$config.SourceTailAbsPcm
$sourceTailRelFraction = [double]$config.SourceTailRelFraction
$maxSlopeRelativeDifference = [double]$config.MaxSlopeRelativeDifference
$requireSourceConvergence = [bool]$config.RequireSourceConvergence
$finalDirectCheck = [bool]$config.FinalDirectCheck
$exeSetting = [string]$config.ExePath

if ($particles -le 0 -or $inactive -lt 0 -or $active -lt 2) {
    throw 'Particles must be positive, Inactive nonnegative, and Active at least 2.'
}
if ($lower -ge $upper -or $initialDepth -lt $lower -or $initialDepth -gt $upper) {
    throw 'Require LowerCm < UpperCm and LowerCm <= InitialDepthCm <= UpperCm.'
}
if ($probeStepDefault -le 0 -or $minProbeStep -le 0 -or $minProbeStep -gt $probeStepDefault) {
    throw 'Require 0 < MinProbeStepCm <= ProbeStepCm.'
}
if ($maxProbeRetries -lt 0 -or $csGenerations -lt 1 -or $csGenerations -gt 10 -or $maxIterations -le 0) {
    throw 'Invalid retry, CS-generation, or iteration count.'
}
if ($damping -le 0 -or $damping -gt 1 -or $maxNewtonStep -le 0 -or $depthTolerance -le 0) {
    throw 'Damping must be in (0,1], and Newton/depth limits must be positive.'
}
if ($kTolerancePcm -lt 0 -or $maxKStderrPcm -lt 0 -or $minAbsSlope -le 0 -or
    $minSlopeSigma -le 0 -or $confidenceSigma -le 0 -or $sourceTailAbsPcm -lt 0 -or
    $sourceTailRelFraction -lt 0 -or $maxSlopeRelativeDifference -lt 0 -or
    ($sourceTailAbsPcm -eq 0 -and $sourceTailRelFraction -eq 0)) {
    throw 'Invalid statistical tolerance or slope threshold.'
}
if ($baseSeed -gt $rngMask) {
    throw 'BaseSeed must fit in the 48-bit RNG state.'
}

function Format-Number([double]$Value, [string]$Format = 'G12') {
    return $Value.ToString($Format, $culture)
}

function Get-ProjectPath([string]$PathValue) {
    if ([IO.Path]::IsPathRooted($PathValue)) {
        return [IO.Path]::GetFullPath($PathValue)
    }
    return [IO.Path]::GetFullPath((Join-Path $rootDir $PathValue))
}

function Get-OddSeed([UInt64]$RawSeed) {
    [UInt64]$value = $RawSeed % $rngModulus
    if (($value % 2) -eq 0) { $value++ }
    return $value
}

function Get-LastDouble([string]$Text, [string]$Pattern, [string]$GroupName = 'value') {
    $found = [regex]::Matches($Text, $Pattern, [Text.RegularExpressions.RegexOptions]::Multiline)
    if ($found.Count -eq 0) { throw "Could not parse pattern: $Pattern" }
    return [double]::Parse($found[$found.Count - 1].Groups[$GroupName].Value, $culture)
}

function Get-LastUInt64([string]$Text, [string]$Pattern, [string]$GroupName = 'value') {
    $found = [regex]::Matches($Text, $Pattern, [Text.RegularExpressions.RegexOptions]::Multiline)
    if ($found.Count -eq 0) { throw "Could not parse pattern: $Pattern" }
    return [UInt64]::Parse($found[$found.Count - 1].Groups[$GroupName].Value, $culture)
}

function Invoke-C5G7 {
    param(
        [double]$ReferenceDepth,
        [AllowNull()][Nullable[double]]$PerturbedDepth,
        [UInt64]$Seed,
        [string]$LogPath
    )

    $programArgs = @(
        '--critical-search-binary',
        "--rod-depth=$(Format-Number $ReferenceDepth)",
        "--particles=$particles",
        "--inactive=$inactive",
        "--active=$active",
        "--seed=$Seed"
    )
    if ($null -ne $PerturbedDepth) {
        $programArgs += "--cs-depth=$(Format-Number ([double]$PerturbedDepth))"
        $programArgs += "--cs-generations=$csGenerations"
    }

    $captured = [System.Collections.Generic.List[string]]::new()
    $utf8NoBom = [Text.UTF8Encoding]::new($false)
    $writer = [IO.StreamWriter]::new($LogPath, $false, $utf8NoBom)
    try {
        & $script:exePath @programArgs 2>&1 | ForEach-Object {
            $line = $_.ToString()
            [void]$captured.Add($line)
            $writer.WriteLine($line)
            $writer.Flush()
            Write-Host $line
        }
        $status = $LASTEXITCODE
    }
    finally {
        $writer.Dispose()
    }
    if ($status -ne 0) {
        throw "C5G7 exited with code $status. See $LogPath"
    }
    return ($captured -join [Environment]::NewLine)
}

$exePath = Get-ProjectPath $exeSetting
$buildScript = Join-Path $scriptDir 'build_windows_nvcc.ps1'
$forceBuild = [bool]$config.ForceBuild
$noBuild = [bool]$config.NoBuild
$needsBuild = -not (Test-Path -LiteralPath $exePath)
if (-not $needsBuild -and -not $noBuild) {
    $exeTime = (Get-Item -LiteralPath $exePath).LastWriteTimeUtc
    $sourceFiles = @(
        Get-ChildItem -Path (Join-Path $rootDir '*.cu') -File
        Get-ChildItem -Path (Join-Path $rootDir '*.cuh') -File
        Get-Item -LiteralPath $buildScript
    )
    $needsBuild = $null -ne ($sourceFiles | Where-Object { $_.LastWriteTimeUtc -gt $exeTime } | Select-Object -First 1)
}

if ($noBuild) {
    if (-not (Test-Path -LiteralPath $exePath)) {
        throw "NoBuild is enabled, but the executable does not exist: $exePath"
    }
    Write-Host "Using existing executable: $exePath"
}
elseif ($forceBuild -or $needsBuild) {
    Write-Host 'Building the local CUDA executable...'
    $buildParams = @{
        NvccPath       = [string]$config.NvccPath
        CudaArch       = [string]$config.CudaArch
        Output         = $exeSetting
        ThreadsPerBlock = [int]$config.ThreadsPerBlock
    }
    if (-not [bool]$config.DeviceLto) { $buildParams.NoDeviceLto = $true }
    if ([bool]$config.FastMath) { $buildParams.FastMath = $true }
    & $buildScript @buildParams
}
else {
    Write-Host "Build is up to date: $exePath"
}

if (-not (Test-Path -LiteralPath $exePath)) {
    throw "Executable was not produced: $exePath"
}

$runTag = Get-Date -Format 'yyyyMMdd_HHmmss'
$logDir = Join-Path $rootDir "results\critical_rod_newton_cs_local_$runTag"
$historyPath = Join-Path $logDir 'newton_history.csv'
$resultPath = Join-Path $logDir 'result.txt'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$depth = [Math]::Min($upper, [Math]::Max($lower, $initialDepth))
$converged = $false
$lastEvaluatedDepth = $depth
$lastK = [double]::NaN
$lastKStderr = [double]::NaN
$lastSlope = [double]::NaN
$lastSlopeStderr = [double]::NaN
$lastLog = ''

Write-Host ''
Write-Host 'Local correlated-sampling Newton search'
Write-Host "  config: $resolvedConfig"
Write-Host "  initial/bounds: $(Format-Number $depth) cm / [$(Format-Number $lower), $(Format-Number $upper)] cm"
Write-Host "  particles/cycles: $particles / $inactive inactive + $active active"
Write-Host "  base seed: $baseSeed (all run seeds are forced odd)"
Write-Host "  source-tail budget: max($(Format-Number $sourceTailAbsPcm) pcm, $(Format-Number (100.0 * $sourceTailRelFraction))% of |delta k|)"
Write-Host "  consecutive-probe slope tolerance: $(Format-Number (100.0 * $maxSlopeRelativeDifference))% or $(Format-Number $confidenceSigma)-sigma"
Write-Host "  output: $logDir"

for ($iteration = 1; $iteration -le $maxIterations; $iteration++) {
    $probeStep = $probeStepDefault
    $retry = 0
    $previousProbeSlope = [double]::NaN
    $previousProbeSlopeStderr = [double]::NaN
    $slopeConsistencyDifference = [double]::NaN
    $slopeConsistencyLimit = [double]::NaN
    $slopeConsistent = $true
    $probeAssessment = 'not-assessed'

    while ($true) {
        if ($depth + $probeStep -le $upper) {
            $referenceDepth = $depth
            $perturbedDepth = $depth + $probeStep
            $currentSide = 'reference'
        }
        else {
            $referenceDepth = [Math]::Max($lower, $depth - $probeStep)
            $perturbedDepth = $depth
            $currentSide = 'perturbed'
        }
        $actualProbeStep = $perturbedDepth - $referenceDepth
        if ($actualProbeStep -le 0) { throw "Cannot form a positive CS stencil at depth $depth." }

        [UInt64]$seedRaw = $baseSeed + [UInt64]($iteration * 104729 + $retry * 1009)
        $seed = Get-OddSeed $seedRaw
        $depthTag = (Format-Number $depth 'F8').Replace('.', 'p')
        $logPath = Join-Path $logDir ("iter_{0:D2}_retry_{1}_d{2}.log" -f $iteration, $retry, $depthTag)

        Write-Host ''
        Write-Host "Iteration $iteration, retry ${retry}: d=$(Format-Number $depth) cm, CS [$(Format-Number $referenceDepth) -> $(Format-Number $perturbedDepth)] cm, seed=$seed"
        $text = Invoke-C5G7 -ReferenceDepth $referenceDepth -PerturbedDepth $perturbedDepth -Seed $seed -LogPath $logPath

        $kReference = Get-LastDouble $text "^Paired reference k:\s*(?<value>$numberPattern)\s*$"
        $kReferenceStderr = Get-LastDouble $text "^k_stderr\(mean\):\s*(?<value>$numberPattern)\s*$"
        $kPerturbed = Get-LastDouble $text "^Perturbed k, NPS \+ fission-source propagation.*:\s*(?<value>$numberPattern)\s*$"
        $deltaK = Get-LastDouble $text "^Total delta k:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$"
        $deltaKStderr = Get-LastDouble $text "^Total delta k:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$" 'stderr'
        $sourceIncrement = Get-LastDouble $text "^Last propagation increment,.*:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$"
        $sourceIncrementStderr = Get-LastDouble $text "^Last propagation increment,.*:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$" 'stderr'
        $diagnosticCount = Get-LastUInt64 $text '^Support/non-finite diagnostic count:\s*(?<value>[0-9]+)\s*$'
        if ($diagnosticCount -ne 0) {
            throw "CS support/non-finite diagnostic count is $diagnosticCount. See $logPath"
        }

        $slope = $deltaK / $actualProbeStep
        $slopeStderr = $deltaKStderr / $actualProbeStep
        $sourceTailPcm = [Math]::Abs($sourceIncrement) * 1.0e5
        $sourceTailUpperPcm = ([Math]::Abs($sourceIncrement) +
            $confidenceSigma * $sourceIncrementStderr) * 1.0e5
        $sourceTailLimitPcm = [Math]::Max(
            $sourceTailAbsPcm,
            $sourceTailRelFraction * [Math]::Abs($deltaK) * 1.0e5)
        $sourceTailFraction = if ([Math]::Abs($deltaK) -gt 0) {
            [Math]::Abs($sourceIncrement / $deltaK)
        }
        else {
            [double]::PositiveInfinity
        }
        $sourceTailSignificant = [Math]::Abs($sourceIncrement) -gt
            $confidenceSigma * $sourceIncrementStderr
        $sourceTailPractical = $sourceTailUpperPcm -le $sourceTailLimitPcm

        $slopeConsistent = $true
        if ($retry -gt 0) {
            $slopeConsistencyDifference = [Math]::Abs($slope - $previousProbeSlope)
            $combinedSlopeStderr = [Math]::Sqrt(
                $slopeStderr * $slopeStderr +
                $previousProbeSlopeStderr * $previousProbeSlopeStderr)
            $slopeScale = [Math]::Max([Math]::Abs($slope), [Math]::Abs($previousProbeSlope))
            $slopeConsistencyLimit = [Math]::Max(
                $confidenceSigma * $combinedSlopeStderr,
                $maxSlopeRelativeDifference * $slopeScale)
            $slopeConsistent = $slopeConsistencyDifference -le $slopeConsistencyLimit
        }

        Write-Host "  source tail: $(Format-Number $sourceTailPcm 'F6') pcm; confidence upper=$(Format-Number $sourceTailUpperPcm 'F6') pcm; allowed=$(Format-Number $sourceTailLimitPcm 'F6') pcm"
        if ($retry -gt 0) {
            Write-Host "  consecutive-probe slope difference=$(Format-Number $slopeConsistencyDifference); allowed=$(Format-Number $slopeConsistencyLimit) per cm"
        }

        $needsSmallerProbe = $false
        $validationReason = ''
        if ($requireSourceConvergence) {
            if (-not $sourceTailPractical) {
                $needsSmallerProbe = $true
                $validationReason = 'the source-tail confidence bound exceeds its practical error budget'
            }
            elseif ($retry -eq 0 -and $sourceTailSignificant) {
                $needsSmallerProbe = $true
                $validationReason = 'the source tail is statistically resolved, so a half-step slope check is required'
            }
            elseif ($retry -gt 0 -and -not $slopeConsistent) {
                $needsSmallerProbe = $true
                $validationReason = 'the consecutive probe slopes are inconsistent'
            }
        }

        if ($needsSmallerProbe) {
            $nextProbe = $probeStep / 2.0
            if ($retry -lt $maxProbeRetries -and $nextProbe -ge $minProbeStep) {
                Write-Warning "$validationReason; validating with h=$(Format-Number $nextProbe) cm."
                $previousProbeSlope = $slope
                $previousProbeSlopeStderr = $slopeStderr
                $probeStep = $nextProbe
                $retry++
                continue
            }
            @(
                'status=CS_PROBE_VALIDATION_FAILED'
                "reason=$validationReason"
                "depth_cm=$(Format-Number $depth)"
                "probe_step_cm=$(Format-Number $actualProbeStep)"
                "source_tail_upper_pcm=$(Format-Number $sourceTailUpperPcm)"
                "source_tail_limit_pcm=$(Format-Number $sourceTailLimitPcm)"
                "slope_consistency_difference=$(Format-Number $slopeConsistencyDifference)"
                "slope_consistency_limit=$(Format-Number $slopeConsistencyLimit)"
                "last_log=$logPath"
            ) | Set-Content -LiteralPath $resultPath -Encoding UTF8
            throw "CS probe validation failed: $validationReason. See $logPath"
        }

        if (-not $requireSourceConvergence) {
            $probeAssessment = 'checks-disabled'
        }
        elseif ($retry -gt 0) {
            $probeAssessment = 'practical-tail-and-slope-consistent'
        }
        else {
            $probeAssessment = 'practical-tail-not-statistically-resolved'
        }
        break
    }

    $slope = $deltaK / $actualProbeStep
    $slopeStderr = $deltaKStderr / $actualProbeStep
    if ($currentSide -eq 'reference') {
        $kCurrent = $kReference
        $kCurrentStderr = $kReferenceStderr
    }
    else {
        $kCurrent = $kPerturbed
        $kCurrentStderr = [Math]::Sqrt($kReferenceStderr * $kReferenceStderr + $deltaKStderr * $deltaKStderr)
    }

    $residual = $kCurrent - $targetK
    $residualPcm = $residual * 1.0e5
    $absResidualPcm = [Math]::Abs($residualPcm)
    $kStderrPcm = $kCurrentStderr * 1.0e5
    $lastEvaluatedDepth = $depth
    $lastK = $kCurrent
    $lastKStderr = $kCurrentStderr
    $lastSlope = $slope
    $lastSlopeStderr = $slopeStderr
    $lastLog = $logPath

    Write-Host "  k(d)=$(Format-Number $kCurrent) +/- $(Format-Number $kCurrentStderr); residual=$(Format-Number $residualPcm 'F6') pcm"
    Write-Host "  dk/dd=$(Format-Number $slope) +/- $(Format-Number $slopeStderr) per cm"

    $stderrOkay = $maxKStderrPcm -le 0 -or $kStderrPcm -le $maxKStderrPcm
    if ($absResidualPcm -le $kTolerancePcm -and $stderrOkay) {
        $converged = $true
        [pscustomobject][ordered]@{
            Iteration = $iteration; Seed = $seed; RequestedDepthCm = $depth
            ReferenceDepthCm = $referenceDepth; PerturbedDepthCm = $perturbedDepth
            ProbeStepCm = $actualProbeStep; KCurrent = $kCurrent; KStderr = $kCurrentStderr
            KReference = $kReference; KPerturbed = $kPerturbed; DeltaK = $deltaK
            DeltaKStderr = $deltaKStderr; SlopePerCm = $slope; SlopeStderrPerCm = $slopeStderr
            ResidualPcm = $residualPcm; LastSourceIncrement = $sourceIncrement
            LastSourceIncrementStderr = $sourceIncrementStderr; SourceTailPcm = $sourceTailPcm
            SourceTailUpperPcm = $sourceTailUpperPcm; SourceTailLimitPcm = $sourceTailLimitPcm
            SourceTailFraction = $sourceTailFraction; ProbeRetry = $retry; ProbeAssessment = $probeAssessment
            SlopeConsistencyDifference = $slopeConsistencyDifference; SlopeConsistencyLimit = $slopeConsistencyLimit
            RawNewtonDepthCm = $null
            NextDepthCm = $null; AppliedStepCm = 0.0; Status = 'CONVERGED'; LogFile = $logPath
        } | Export-Csv -LiteralPath $historyPath -NoTypeInformation -Append
        Write-Host "Converged at d=$(Format-Number $depth) cm."
        break
    }

    if ([Math]::Abs($slope) -lt $minAbsSlope) {
        throw "|dk/dd| is below MinAbsSlope=$minAbsSlope."
    }
    if ($slope -ge 0) {
        throw "The rod-worth slope is nonnegative: $slope."
    }
    if ($slopeStderr -gt 0 -and [Math]::Abs($slope) -lt $minSlopeSigma * $slopeStderr) {
        throw 'The CS slope is not statistically significant. Increase ProbeStepCm, Particles, or Active.'
    }

    $rawNewtonDepth = $depth - $residual / $slope
    $rawStep = $rawNewtonDepth - $depth
    $limitedStep = [Math]::Min($maxNewtonStep, [Math]::Max(-$maxNewtonStep, $rawStep))
    $nextDepth = $depth + $damping * $limitedStep
    $nextDepth = [Math]::Min($upper, [Math]::Max($lower, $nextDepth))
    $appliedStep = $nextDepth - $depth

    [pscustomobject][ordered]@{
        Iteration = $iteration; Seed = $seed; RequestedDepthCm = $depth
        ReferenceDepthCm = $referenceDepth; PerturbedDepthCm = $perturbedDepth
        ProbeStepCm = $actualProbeStep; KCurrent = $kCurrent; KStderr = $kCurrentStderr
        KReference = $kReference; KPerturbed = $kPerturbed; DeltaK = $deltaK
        DeltaKStderr = $deltaKStderr; SlopePerCm = $slope; SlopeStderrPerCm = $slopeStderr
        ResidualPcm = $residualPcm; LastSourceIncrement = $sourceIncrement
        LastSourceIncrementStderr = $sourceIncrementStderr; SourceTailPcm = $sourceTailPcm
        SourceTailUpperPcm = $sourceTailUpperPcm; SourceTailLimitPcm = $sourceTailLimitPcm
        SourceTailFraction = $sourceTailFraction; ProbeRetry = $retry; ProbeAssessment = $probeAssessment
        SlopeConsistencyDifference = $slopeConsistencyDifference; SlopeConsistencyLimit = $slopeConsistencyLimit
        RawNewtonDepthCm = $rawNewtonDepth
        NextDepthCm = $nextDepth; AppliedStepCm = $appliedStep; Status = 'CONTINUE'; LogFile = $logPath
    } | Export-Csv -LiteralPath $historyPath -NoTypeInformation -Append

    Write-Host "  raw Newton depth=$(Format-Number $rawNewtonDepth) cm"
    Write-Host "  accepted next depth=$(Format-Number $nextDepth) cm"

    if ([Math]::Abs($appliedStep) -le $depthTolerance) {
        Write-Warning 'Newton step reached the depth tolerance before the k/statistical tolerance.'
        break
    }
    $depth = $nextDepth
}

if (-not $converged) {
    @(
        'status=NOT_CONVERGED'
        "last_evaluated_depth_cm=$(Format-Number $lastEvaluatedDepth)"
        "last_k=$(Format-Number $lastK)"
        "last_k_stderr=$(Format-Number $lastKStderr)"
        "last_slope_per_cm=$(Format-Number $lastSlope)"
        "source_tail_upper_pcm=$(Format-Number $sourceTailUpperPcm)"
        "source_tail_limit_pcm=$(Format-Number $sourceTailLimitPcm)"
        "probe_assessment=$probeAssessment"
        "history=$historyPath"
        "last_log=$lastLog"
    ) | Set-Content -LiteralPath $resultPath -Encoding UTF8
    Write-Error "Newton search did not converge. See $resultPath"
    exit 4
}

$verificationStatus = 'SKIPPED'
$verificationK = [double]::NaN
$verificationStderr = [double]::NaN
$verificationErrorPcm = [double]::NaN
$verificationLog = ''
if ($finalDirectCheck) {
    [UInt64]$verificationRaw = $baseSeed + [UInt64](($maxIterations + 1) * 104729)
    $verificationSeed = Get-OddSeed $verificationRaw
    $depthTag = (Format-Number $depth 'F8').Replace('.', 'p')
    $verificationLog = Join-Path $logDir "final_direct_d$depthTag.log"
    Write-Host ''
    Write-Host "Independent direct verification at d=$(Format-Number $depth) cm, seed=$verificationSeed"
    $verificationText = Invoke-C5G7 -ReferenceDepth $depth -PerturbedDepth $null -Seed $verificationSeed -LogPath $verificationLog
    $verificationK = Get-LastDouble $verificationText "^k_mean:\s*(?<value>$numberPattern)\s*$"
    $verificationStderr = Get-LastDouble $verificationText "^k_stderr\(mean\):\s*(?<value>$numberPattern)\s*$"
    $verificationErrorPcm = [Math]::Abs($verificationK - $targetK) * 1.0e5
    $verificationLimitPcm = $kTolerancePcm + $confidenceSigma * $verificationStderr * 1.0e5
    $verificationStatus = if ($verificationErrorPcm -le $verificationLimitPcm) { 'PASSED' } else { 'FAILED' }
    Write-Host "  direct k=$(Format-Number $verificationK) +/- $(Format-Number $verificationStderr)"
    Write-Host "  |k-target|=$(Format-Number $verificationErrorPcm 'F6') pcm; limit=$(Format-Number $verificationLimitPcm 'F6') pcm"
}

$depthUncertainty = $lastKStderr / [Math]::Abs($lastSlope)
@(
    'status=CONVERGED'
    "critical_depth_cm=$(Format-Number $depth)"
    "cs_k=$(Format-Number $lastK)"
    "cs_k_stderr=$(Format-Number $lastKStderr)"
    "cs_slope_per_cm=$(Format-Number $lastSlope)"
    "cs_slope_stderr_per_cm=$(Format-Number $lastSlopeStderr)"
    "source_tail_upper_pcm=$(Format-Number $sourceTailUpperPcm)"
    "source_tail_limit_pcm=$(Format-Number $sourceTailLimitPcm)"
    "probe_assessment=$probeAssessment"
    "approx_depth_stderr_cm=$(Format-Number $depthUncertainty)"
    "target_k=$(Format-Number $targetK)"
    "verification_status=$verificationStatus"
    "verification_k=$(Format-Number $verificationK)"
    "verification_k_stderr=$(Format-Number $verificationStderr)"
    "verification_abs_error_pcm=$(Format-Number $verificationErrorPcm)"
    "history=$historyPath"
    "last_cs_log=$lastLog"
    "verification_log=$verificationLog"
) | Set-Content -LiteralPath $resultPath -Encoding UTF8

Write-Host ''
Write-Host "Critical rod depth: $(Format-Number $depth) cm"
Write-Host "Approximate one-sigma depth uncertainty: $(Format-Number $depthUncertainty) cm"
Write-Host "Verification: $verificationStatus"
Write-Host "Result: $resultPath"

if ($verificationStatus -eq 'FAILED') { exit 5 }
