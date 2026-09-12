param(
    [double]$LowerCm = 82.0,
    [double]$UpperCm = 84.0,
    [double]$TargetK = 1.0,
    [double]$DepthToleranceCm = 0.01,
    [double]$KTolerancePcm = 0.0,
    [double]$AroCheckDepthCm = 21.42,
    [double]$AroExpectedK = 1.14390,
    [double]$AroMinTolerancePcm = 20.0,
    [int]$ReplicatesPerDepth = 3,
    [int]$MaxReplicateFailures = 2,
    [UInt64]$BaseSeed = 92235922383,
    [double]$ConfidenceSigma = 2.0,
    [double]$ScatterWarningRatio = 2.0,
    [int]$MaxIterations = 10,
    [int]$Particles = 1000000,
    [int]$Inactive = 150,
    [int]$Active = 200,
    [string]$CudaArch = "sm_120",
    [string]$NvccPath = "nvcc",
    [string]$ExePath = "build\windows_nvcc\c5g7.exe",
    [ValidateSet(32, 64, 128, 256, 512, 1024)]
    [int]$ThreadsPerBlock = 32,
    [string]$LogDir = "",
    [string[]]$CaseArgs = @("--critical-search-binary"),
    [switch]$NoBuild,
    [switch]$NoVsDevEnv,
    [switch]$AllowUnsupportedCompiler,
    [switch]$NoDeviceLto,
    [switch]$FastMath,
    [switch]$SkipAroCheck
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
if ($PSVersionTable.PSVersion.Major -ge 7) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$RootDir = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RootDir
$Culture = [System.Globalization.CultureInfo]::InvariantCulture

if ($LowerCm -ge $UpperCm) {
    throw "LowerCm must be less than UpperCm."
}
if ($DepthToleranceCm -le 0.0) {
    throw "DepthToleranceCm must be greater than zero."
}
if ($KTolerancePcm -lt 0.0) {
    throw "KTolerancePcm cannot be negative."
}
if ($AroCheckDepthCm -lt 0.0 -or $AroExpectedK -le 0.0 -or $AroMinTolerancePcm -lt 0.0) {
    throw "ARO check depth/tolerance cannot be negative, and AroExpectedK must be positive."
}
if ($ReplicatesPerDepth -lt 2 -or $MaxReplicateFailures -lt 0) {
    throw "ReplicatesPerDepth must be at least 2, and MaxReplicateFailures cannot be negative."
}
if ($ConfidenceSigma -le 0.0 -or $ScatterWarningRatio -le 0.0) {
    throw "ConfidenceSigma and ScatterWarningRatio must be greater than zero."
}
if ($MaxIterations -lt 1) {
    throw "MaxIterations must be at least 1."
}
if ($Particles -lt 1 -or $Inactive -lt 0 -or $Active -lt 2) {
    throw "Particles must be positive, Inactive cannot be negative, and Active must be at least 2."
}

if (-not $NoBuild) {
    $buildArgs = @{
        NvccPath = $NvccPath
        CudaArch = $CudaArch
        Output = $ExePath
        ThreadsPerBlock = $ThreadsPerBlock
    }
    if ($NoVsDevEnv) { $buildArgs.NoVsDevEnv = $true }
    if ($AllowUnsupportedCompiler) { $buildArgs.AllowUnsupportedCompiler = $true }
    if ($NoDeviceLto) { $buildArgs.NoDeviceLto = $true }
    if ($FastMath) { $buildArgs.FastMath = $true }

    & (Join-Path $PSScriptRoot "build_windows_nvcc.ps1") @buildArgs
}

$exeCandidate = if ([System.IO.Path]::IsPathRooted($ExePath)) {
    $ExePath
}
else {
    Join-Path $RootDir $ExePath
}
if (-not (Test-Path -LiteralPath $exeCandidate)) {
    throw "Executable not found: $exeCandidate. Run without -NoBuild to compile it first."
}
$Exe = (Resolve-Path -LiteralPath $exeCandidate).Path

if ([string]::IsNullOrWhiteSpace($LogDir)) {
    $runTag = Get-Date -Format "yyyyMMdd_HHmmss"
    $LogDir = "results\critical_rod_binary_search_windows_nvcc_$runTag"
}
$logDirPath = if ([System.IO.Path]::IsPathRooted($LogDir)) {
    $LogDir
}
else {
    Join-Path $RootDir $LogDir
}
New-Item -ItemType Directory -Force -Path $logDirPath | Out-Null

$Summary = Join-Path $logDirPath "summary.csv"
$ReplicateSummary = Join-Path $logDirPath "replicates.csv"
$Result = Join-Path $logDirPath "result.txt"
"evaluation,requested_depth_cm,actual_depth_cm,k_mean,k_stderr_conservative,k_stderr_within,k_stderr_between,replicates,scatter_ratio,status,seeds" | Set-Content -LiteralPath $Summary
"evaluation,attempt,seed,requested_depth_cm,actual_depth_cm,k_mean,k_stderr_mean,status,log_file" | Set-Content -LiteralPath $ReplicateSummary

$script:Evaluation = 0
$script:BestRun = $null

function Format-Number {
    param(
        [double]$Value,
        [string]$Format = "F8"
    )
    return $Value.ToString($Format, $Culture)
}

function Read-LastValue {
    param(
        [string[]]$Lines,
        [string]$Pattern,
        [string]$Replacement
    )

    $match = $Lines | Select-String -Pattern $Pattern | Select-Object -Last 1
    if (-not $match) {
        return ""
    }
    return ($match.Line -replace $Replacement, "").Trim()
}

function Update-BestRun {
    param([pscustomobject]$Run)

    if (-not $Run.Valid) {
        return
    }
    $delta = [math]::Abs($Run.KMean - $TargetK)
    if ($null -eq $script:BestRun -or $delta -lt $script:BestRun.Delta) {
        $script:BestRun = [pscustomobject]@{
            Depth = $Run.Depth
            Actual = $Run.Actual
            KMean = $Run.KMean
            KStderr = $Run.KStderr
            WithinStderr = $Run.WithinStderr
            BetweenStderr = $Run.BetweenStderr
            ScatterRatio = $Run.ScatterRatio
            ReplicateCount = $Run.ReplicateCount
            Delta = $delta
        }
    }
}

function Get-IndependentSeed {
    param(
        [double]$Depth,
        [int]$Attempt
    )

    [UInt64]$modulus = 281474976710656
    [UInt64]$depthKey = [math]::Round($Depth * 1000000.0)
    [UInt64]$depthOffset = ($depthKey * [UInt64]1000003) % $modulus
    [UInt64]$attemptOffset = ([UInt64]$Attempt * [UInt64]104729) % $modulus
    [UInt64]$seed = (($BaseSeed % $modulus) + $depthOffset + $attemptOffset) % $modulus
    if (($seed % 2) -eq 0) { $seed++ }
    return $seed
}

function Invoke-SingleRodCase {
    param(
        [double]$Depth,
        [int]$Evaluation,
        [int]$Attempt,
        [UInt64]$Seed
    )

    $depthText = Format-Number $Depth
    $depthTag = $depthText.Replace(".", "p")
    $logFile = Join-Path $logDirPath ("eval_{0:D2}_attempt_{1:D2}_seed_{2}_rod_{3}cm.log" -f $Evaluation, $Attempt, $Seed, $depthTag)

    Write-Host "    replicate attempt $Attempt, seed $Seed"
    $programArgs = @()
    $programArgs += $CaseArgs
    $programArgs += "--rod-depth=$depthText"
    $programArgs += "--particles=$Particles"
    $programArgs += "--inactive=$Inactive"
    $programArgs += "--active=$Active"
    $programArgs += "--seed=$Seed"

    $statusCode = 0
    $previousErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & $Exe @programArgs > $logFile 2>&1
        $statusCode = $LASTEXITCODE
        if ($null -eq $statusCode) { $statusCode = 0 }
    }
    catch {
        $_ | Out-File -LiteralPath $logFile -Append
        $statusCode = 1
    }
    finally {
        $ErrorActionPreference = $previousErrorActionPreference
    }

    $lines = if (Test-Path -LiteralPath $logFile) {
        @(Get-Content -LiteralPath $logFile)
    }
    else {
        @()
    }
    $actualText = Read-LastValue $lines "Actual (?:modeled|quantized) rod depth:" ".*Actual (?:modeled|quantized) rod depth:\s*|\s*cm.*"
    $kMeanText = Read-LastValue $lines "^k_mean:" "^k_mean:\s*"
    $kStderrText = Read-LastValue $lines "^k_stderr\(mean\):" "^k_stderr\(mean\):\s*"

    $kMean = 0.0
    $kStderr = [double]::NaN
    $validK = [double]::TryParse($kMeanText, [System.Globalization.NumberStyles]::Float, $Culture, [ref]$kMean)
    if (-not [string]::IsNullOrWhiteSpace($kStderrText)) {
        [void][double]::TryParse($kStderrText, [System.Globalization.NumberStyles]::Float, $Culture, [ref]$kStderr)
    }

    $valid = ($statusCode -eq 0 -and $validK)
    $statusText = if ($valid) { "OK" } else { "FAILED:$statusCode" }
    "{0},{1},{2},{3},{4},{5},{6},{7},{8}" -f $Evaluation, $Attempt, $Seed, $depthText, $actualText, $kMeanText, $kStderrText, $statusText, $logFile | Add-Content -LiteralPath $ReplicateSummary

    $run = [pscustomobject]@{
        Valid = $valid
        Depth = [double]$depthText
        Actual = $actualText
        KMean = $kMean
        KStderr = $kStderr
        LogFile = $logFile
        StatusCode = $statusCode
        Seed = $Seed
    }

    if ($valid) {
        Write-Host "        k = $kMeanText, stderr = $kStderrText"
    }
    else {
        Write-Warning "Depth $depthText seed $Seed failed with exit code $statusCode; keeping $logFile."
    }
    return $run
}

function Invoke-RodCase {
    param([double]$Depth)

    $script:Evaluation++
    $evaluation = $script:Evaluation
    $depthText = Format-Number $Depth
    Write-Host "==> evaluation ${evaluation}: rod depth $depthText cm, averaging $ReplicatesPerDepth independent runs"

    $runs = New-Object System.Collections.Generic.List[object]
    $maxAttempts = $ReplicatesPerDepth + $MaxReplicateFailures
    for ($attempt = 1; $attempt -le $maxAttempts -and $runs.Count -lt $ReplicatesPerDepth; ++$attempt) {
        $seed = Get-IndependentSeed $Depth $attempt
        $run = Invoke-SingleRodCase $Depth $evaluation $attempt $seed
        if ($run.Valid) {
            $runs.Add($run)
        }
    }

    if ($runs.Count -lt $ReplicatesPerDepth) {
        $seeds = (@($runs | ForEach-Object { $_.Seed }) -join ";")
        "{0},{1},,NaN,NaN,NaN,NaN,{2},NaN,FAILED,{3}" -f $evaluation, $depthText, $runs.Count, $seeds | Add-Content -LiteralPath $Summary
        return [pscustomobject]@{
            Valid = $false
            Depth = $Depth
            Actual = ""
            KMean = [double]::NaN
            KStderr = [double]::NaN
            WithinStderr = [double]::NaN
            BetweenStderr = [double]::NaN
            ScatterRatio = [double]::NaN
            ReplicateCount = $runs.Count
            LogFile = (@($runs | ForEach-Object { $_.LogFile }) -join ";")
        }
    }

    $count = $runs.Count
    $kMean = [double](($runs | Measure-Object -Property KMean -Average).Average)
    $sumSquaredDeviation = 0.0
    $sumStderrSquared = 0.0
    $finiteStderrCount = 0
    foreach ($run in $runs) {
        $difference = $run.KMean - $kMean
        $sumSquaredDeviation += $difference * $difference
        if (-not [double]::IsNaN($run.KStderr)) {
            $sumStderrSquared += $run.KStderr * $run.KStderr
            $finiteStderrCount++
        }
    }

    $sampleStdDev = [math]::Sqrt($sumSquaredDeviation / ($count - 1))
    $betweenStderr = $sampleStdDev / [math]::Sqrt($count)
    $withinStderr = if ($finiteStderrCount -gt 0) {
        [math]::Sqrt($sumStderrSquared) / $count
    }
    else {
        [double]::NaN
    }
    $individualRmsStderr = if ($finiteStderrCount -gt 0) {
        [math]::Sqrt($sumStderrSquared / $finiteStderrCount)
    }
    else {
        [double]::NaN
    }
    $scatterRatio = if (-not [double]::IsNaN($individualRmsStderr) -and $individualRmsStderr -gt 0.0) {
        $sampleStdDev / $individualRmsStderr
    }
    else {
        [double]::NaN
    }
    $conservativeStderr = if ([double]::IsNaN($withinStderr)) {
        $betweenStderr
    }
    else {
        [math]::Max($withinStderr, $betweenStderr)
    }

    $actual = $runs[0].Actual
    $seeds = (@($runs | ForEach-Object { $_.Seed }) -join ";")
    "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}" -f $evaluation, $depthText, $actual, (Format-Number $kMean F7), (Format-Number $conservativeStderr G), (Format-Number $withinStderr G), (Format-Number $betweenStderr G), $count, (Format-Number $scatterRatio F3), "OK", $seeds | Add-Content -LiteralPath $Summary

    $aggregate = [pscustomobject]@{
        Valid = $true
        Depth = [double]$depthText
        Actual = $actual
        KMean = $kMean
        KStderr = $conservativeStderr
        WithinStderr = $withinStderr
        BetweenStderr = $betweenStderr
        ScatterRatio = $scatterRatio
        ReplicateCount = $count
        LogFile = (@($runs | ForEach-Object { $_.LogFile }) -join ";")
        Seeds = $seeds
    }

    Write-Host "    averaged k = $(Format-Number $kMean F7) +/- $(Format-Number $conservativeStderr G)"
    Write-Host "    within-run SE = $(Format-Number $withinStderr G); between-run SE = $(Format-Number $betweenStderr G); scatter ratio = $(Format-Number $scatterRatio F2)"
    if (-not [double]::IsNaN($scatterRatio) -and $scatterRatio -gt $ScatterWarningRatio) {
        Write-Warning "Replicate scatter is $(Format-Number $scatterRatio F2)x the runs' reported stderr. This may indicate GPU nondeterminism or underestimated uncertainty."
    }

    Update-BestRun $aggregate
    return $aggregate
}

function Get-EndpointRun {
    param(
        [double]$Depth,
        [bool]$IsLower
    )

    $run = Invoke-RodCase $Depth
    if ($run.Valid) {
        return $run
    }

    $fallback = if ($IsLower) {
        $LowerCm + 0.05 * ($UpperCm - $LowerCm)
    }
    else {
        $UpperCm - 0.05 * ($UpperCm - $LowerCm)
    }
    Write-Warning "Endpoint failed; trying nearby depth $(Format-Number $fallback) cm."
    return Invoke-RodCase $fallback
}

function Test-Bracket {
    param(
        [double]$FirstK,
        [double]$SecondK
    )
    $first = $FirstK - $TargetK
    $second = $SecondK - $TargetK
    return (($first -le 0.0 -and $second -ge 0.0) -or ($first -ge 0.0 -and $second -le 0.0))
}

Write-Host "CriticalSearch_Binary rod-depth search (Windows direct NVCC)"
Write-Host "  exe: $Exe"
Write-Host "  initial bracket: $(Format-Number $LowerCm F2) to $(Format-Number $UpperCm F2) cm"
Write-Host "  fuel insertion below top moderator: $(Format-Number ($LowerCm - 21.42) F2) to $(Format-Number ($UpperCm - 21.42) F2) cm"
Write-Host "  target k: $(Format-Number $TargetK F7)"
Write-Host "  depth tolerance: $DepthToleranceCm cm; k tolerance: $KTolerancePcm pcm"
Write-Host "  particles/inactive/active: $Particles/$Inactive/$Active"
Write-Host "  independent replicates/depth: $ReplicatesPerDepth; base seed: $BaseSeed"
Write-Host "  statistical decision interval: +/- $ConfidenceSigma aggregate stderr"
Write-Host "  output: $logDirPath"

$aroCheckStatus = "SKIPPED"
$aroCheckK = ""
$aroCheckStderr = ""
$aroCheckDeltaPcm = ""
$aroCheckLimitPcm = ""
if (-not $SkipAroCheck) {
    Write-Host ""
    Write-Host "Running ARO equivalence check at $(Format-Number $AroCheckDepthCm F2) cm..."
    $aroRun = Invoke-RodCase $AroCheckDepthCm
    if (-not $aroRun.Valid) {
        @(
            "STATUS=ARO_CHECK_FAILED"
            "reason=simulation failed"
            "aro_check_depth_cm=$(Format-Number $AroCheckDepthCm)"
            "aro_expected_k=$(Format-Number $AroExpectedK F7)"
            "log=$($aroRun.LogFile)"
        ) | Tee-Object -FilePath $Result
        exit 5
    }

    $aroCheckK = Format-Number $aroRun.KMean F7
    $aroCheckStderr = Format-Number $aroRun.KStderr G
    $aroCheckDeltaPcmValue = [math]::Abs((($aroRun.KMean - 1.0) / $aroRun.KMean) - (($AroExpectedK - 1.0) / $AroExpectedK)) * 100000.0
    $aroSigmaPcm = if ([double]::IsNaN($aroRun.KStderr)) {
        0.0
    }
    else {
        $aroRun.KStderr / ($aroRun.KMean * $aroRun.KMean) * 100000.0
    }
    $aroCheckLimitPcmValue = [math]::Max($AroMinTolerancePcm, 3.0 * $aroSigmaPcm)
    $aroCheckDeltaPcm = Format-Number $aroCheckDeltaPcmValue F3
    $aroCheckLimitPcm = Format-Number $aroCheckLimitPcmValue F3

    Write-Host "    expected k = $(Format-Number $AroExpectedK F7)"
    Write-Host "    observed k = $aroCheckK +/- $aroCheckStderr"
    Write-Host "    difference = $aroCheckDeltaPcm pcm; allowed = $aroCheckLimitPcm pcm"

    if ($aroCheckDeltaPcmValue -gt $aroCheckLimitPcmValue) {
        @(
            "STATUS=ARO_CHECK_FAILED"
            "reason=multiplication factor mismatch"
            "aro_check_depth_cm=$(Format-Number $AroCheckDepthCm)"
            "aro_expected_k=$(Format-Number $AroExpectedK F7)"
            "aro_observed_k=$aroCheckK"
            "aro_k_stderr=$aroCheckStderr"
            "aro_difference_pcm=$aroCheckDeltaPcm"
            "aro_allowed_pcm=$aroCheckLimitPcm"
            "log=$($aroRun.LogFile)"
        ) | Tee-Object -FilePath $Result
        exit 5
    }

    $aroCheckStatus = "PASSED"
    Write-Host "    ARO equivalence check passed."
    Write-Host ""
}

$lowRun = Get-EndpointRun $LowerCm $true
if (-not $lowRun.Valid) {
    "STATUS=FAILED`nCould not obtain a valid lower-endpoint result." | Tee-Object -FilePath $Result
    exit 3
}

$highRun = Get-EndpointRun $UpperCm $false
if (-not $highRun.Valid) {
    "STATUS=FAILED`nCould not obtain a valid upper-endpoint result." | Tee-Object -FilePath $Result
    exit 3
}

if (-not (Test-Bracket $lowRun.KMean $highRun.KMean)) {
    @(
        "STATUS=NOT_BRACKETED"
        "The initial interval does not bracket target k=$TargetK."
        "lower_depth_cm=$(Format-Number $lowRun.Depth) lower_k=$(Format-Number $lowRun.KMean F7)"
        "upper_depth_cm=$(Format-Number $highRun.Depth) upper_k=$(Format-Number $highRun.KMean F7)"
        "Expand -LowerCm/-UpperCm and run again."
    ) | Tee-Object -FilePath $Result
    exit 4
}

$stopReason = "maximum iterations reached"
$searchStatus = "COMPLETE"
for ($iteration = 1; $iteration -le $MaxIterations; ++$iteration) {
    $width = $highRun.Depth - $lowRun.Depth
    if ($width -le $DepthToleranceCm) {
        $stopReason = "depth bracket tolerance reached"
        break
    }

    $midpoint = ($lowRun.Depth + $highRun.Depth) / 2.0
    $midRun = Invoke-RodCase $midpoint
    if (-not $midRun.Valid) {
        $leftCandidate = $lowRun.Depth + 0.45 * $width
        $rightCandidate = $lowRun.Depth + 0.55 * $width
        Write-Warning "Midpoint failed; trying $(Format-Number $leftCandidate) cm and then $(Format-Number $rightCandidate) cm."
        $midRun = Invoke-RodCase $leftCandidate
        if (-not $midRun.Valid) {
            $midRun = Invoke-RodCase $rightCandidate
        }
    }

    if (-not $midRun.Valid) {
        $stopReason = "three candidate simulations failed in the current bracket"
        $searchStatus = "PARTIAL"
        break
    }

    $confidenceHalfWidth = $ConfidenceSigma * $midRun.KStderr
    if (-not [double]::IsNaN($confidenceHalfWidth) -and [math]::Abs($midRun.KMean - $TargetK) -le $confidenceHalfWidth) {
        $stopReason = "target k lies inside the averaged confidence interval"
        $lowRun = $midRun
        $highRun = $midRun
        break
    }

    $deltaPcm = [math]::Abs($midRun.KMean - $TargetK) * 100000.0
    if ($KTolerancePcm -gt 0.0 -and $deltaPcm -le $KTolerancePcm) {
        $stopReason = "k tolerance reached"
        $lowRun = $midRun
        $highRun = $midRun
        break
    }

    if (Test-Bracket $lowRun.KMean $midRun.KMean) {
        $highRun = $midRun
    }
    else {
        $lowRun = $midRun
    }

    Write-Host "    bracket is now $(Format-Number $lowRun.Depth) cm (k=$(Format-Number $lowRun.KMean F7)) to $(Format-Number $highRun.Depth) cm (k=$(Format-Number $highRun.KMean F7))"
}

if ($lowRun.KMean -ne $highRun.KMean) {
    $estimatedDepth = $lowRun.Depth + ($TargetK - $lowRun.KMean) * ($highRun.Depth - $lowRun.Depth) / ($highRun.KMean - $lowRun.KMean)
    $estimatedDepth = [math]::Max($lowRun.Depth, [math]::Min($highRun.Depth, $estimatedDepth))
}
else {
    $estimatedDepth = $lowRun.Depth
}

$resultLines = @(
    "STATUS=$searchStatus"
    "stop_reason=$stopReason"
    "target_k=$(Format-Number $TargetK F7)"
    "replicates_per_depth=$ReplicatesPerDepth"
    "base_seed=$BaseSeed"
    "confidence_sigma=$(Format-Number $ConfidenceSigma F2)"
    "aro_check_status=$aroCheckStatus"
    "aro_check_depth_cm=$(Format-Number $AroCheckDepthCm)"
    "aro_expected_k=$(Format-Number $AroExpectedK F7)"
    "aro_observed_k=$aroCheckK"
    "aro_k_stderr=$aroCheckStderr"
    "aro_difference_pcm=$aroCheckDeltaPcm"
    "aro_allowed_pcm=$aroCheckLimitPcm"
    "final_lower_depth_cm=$(Format-Number $lowRun.Depth) final_lower_k=$(Format-Number $lowRun.KMean F7)"
    "final_upper_depth_cm=$(Format-Number $highRun.Depth) final_upper_k=$(Format-Number $highRun.KMean F7)"
    "estimated_critical_depth_cm=$(Format-Number $estimatedDepth)"
    "estimated_fuel_insertion_below_top_moderator_cm=$(Format-Number ($estimatedDepth - 21.42))"
    "best_sampled_depth_cm=$(Format-Number $script:BestRun.Depth)"
    "best_sampled_actual_depth_cm=$($script:BestRun.Actual)"
    "best_sampled_k=$(Format-Number $script:BestRun.KMean F7)"
    "best_sampled_k_stderr=$(Format-Number $script:BestRun.KStderr G)"
    "best_sampled_within_stderr=$(Format-Number $script:BestRun.WithinStderr G)"
    "best_sampled_between_stderr=$(Format-Number $script:BestRun.BetweenStderr G)"
    "best_sampled_scatter_ratio=$(Format-Number $script:BestRun.ScatterRatio F3)"
    "summary=$Summary"
    "replicates=$ReplicateSummary"
)
$resultLines | Tee-Object -FilePath $Result
