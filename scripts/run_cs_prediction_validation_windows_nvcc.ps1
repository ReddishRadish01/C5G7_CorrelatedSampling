[CmdletBinding()]
param(
    [string]$ConfigPath = '',
    [string]$OutputDir = '',
    [switch]$SmokeTest,
    [switch]$NoBuild,
    [switch]$ForceBuild
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrWhiteSpace($ConfigPath)) {
    $ConfigPath = Join-Path $scriptDir 'cs_prediction_validation_windows_config.psd1'
}

$culture = [Globalization.CultureInfo]::InvariantCulture
$numberPattern = '[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?'
$rngModulus = [UInt64]281474976710656
$rngMask = [UInt64]281474976710655
$rootDir = (Resolve-Path (Join-Path $scriptDir '..')).Path
Set-Location $rootDir

$resolvedConfig = (Resolve-Path -LiteralPath $ConfigPath).Path
$config = Import-PowerShellDataFile -LiteralPath $resolvedConfig
$requiredKeys = @(
    'ReferenceDepthCm', 'ProbeStepsCm', 'ModeledUpperDepthCm',
    'Particles', 'Inactive', 'Active', 'CsGenerations', 'BaseSeed',
    'ConfidenceSigma', 'ExpectedDeltaKSign',
    'SourceTailAbsPcm', 'SourceTailRelFraction', 'ProgressEveryCycles',
    'ExePath', 'NvccPath', 'CudaArch', 'ThreadsPerBlock',
    'DeviceLto', 'FastMath', 'ForceBuild', 'NoBuild'
)
foreach ($key in $requiredKeys) {
    if (-not $config.ContainsKey($key)) {
        throw "Missing '$key' in configuration file: $resolvedConfig"
    }
}

$referenceDepth = [double]$config.ReferenceDepthCm
$probeSteps = [double[]]@($config.ProbeStepsCm)
$modeledUpperDepth = [double]$config.ModeledUpperDepthCm
$particles = [int]$config.Particles
$inactive = [int]$config.Inactive
$active = [int]$config.Active
$csGenerations = [int]$config.CsGenerations
$baseSeed = [UInt64]$config.BaseSeed
$confidenceSigma = [double]$config.ConfidenceSigma
$expectedDeltaKSign = [int]$config.ExpectedDeltaKSign
$sourceTailAbsPcm = [double]$config.SourceTailAbsPcm
$sourceTailRelFraction = [double]$config.SourceTailRelFraction
$progressEveryCycles = [int]$config.ProgressEveryCycles
$exeSetting = [string]$config.ExePath

if ($SmokeTest) {
    $probeSteps = [double[]]@(0.125)
    $particles = 20000
    $inactive = 10
    $active = 12
    $progressEveryCycles = 2
}

if ($particles -le 0 -or $inactive -lt 0 -or $active -lt 2) {
    throw 'Particles must be positive, Inactive nonnegative, and Active at least 2.'
}
if ($csGenerations -lt 1 -or $csGenerations -gt 10) {
    throw 'CsGenerations must be between 1 and 10.'
}
if ($active -lt 2 -or ($inactive + $active) -le $csGenerations) {
    throw 'There are not enough cycles to produce CS samples at the selected propagation level.'
}
if ($probeSteps.Count -eq 0) { throw 'ProbeStepsCm must contain at least one step.' }
if ($expectedDeltaKSign -ne -1 -and $expectedDeltaKSign -ne 1) {
    throw 'ExpectedDeltaKSign must be -1 or +1.'
}
if ($confidenceSigma -le 0 -or $sourceTailAbsPcm -lt 0 -or $sourceTailRelFraction -lt 0) {
    throw 'Invalid confidence or source-tail setting.'
}
if ($baseSeed -gt $rngMask) { throw 'BaseSeed must fit in the 48-bit RNG state.' }

$seenSteps = @{}
foreach ($step in $probeSteps) {
    if ($step -le 0) {
        throw 'This implementation validates positive insertion steps only; every ProbeStepsCm value must be > 0.'
    }
    if ($referenceDepth + $step -gt $modeledUpperDepth) {
        throw "Target depth $($referenceDepth + $step) cm exceeds ModeledUpperDepthCm=$modeledUpperDepth."
    }
    $key = $step.ToString('G17', $culture)
    if ($seenSteps.ContainsKey($key)) { throw "Duplicate probe step: $key cm" }
    $seenSteps[$key] = $true
}
$probeSteps = [double[]]@($probeSteps | Sort-Object)

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

function Get-LastInt64([string]$Text, [string]$Pattern, [string]$GroupName = 'value') {
    $found = [regex]::Matches($Text, $Pattern, [Text.RegularExpressions.RegexOptions]::Multiline)
    if ($found.Count -eq 0) { throw "Could not parse pattern: $Pattern" }
    return [Int64]::Parse($found[$found.Count - 1].Groups[$GroupName].Value, $culture)
}

function Read-DirectResult([string]$LogPath) {
    $text = Get-Content -LiteralPath $LogPath -Raw
    return [pscustomobject]@{
        RequestedDepthCm = Get-LastDouble $text "^Requested rod depth:\s*(?<value>$numberPattern)\s+cm\s*$"
        ActualDepthCm = Get-LastDouble $text "^Actual modeled rod depth:\s*(?<value>$numberPattern)\s+cm\s*$"
        K = Get-LastDouble $text "^k_mean:\s*(?<value>$numberPattern)\s*$"
        KStderr = Get-LastDouble $text "^k_stderr\(mean\):\s*(?<value>$numberPattern)\s*$"
        ActiveCycles = Get-LastInt64 $text '^Active cycles:\s*(?<value>[0-9]+)\s*$'
    }
}

function Read-CsResult([string]$LogPath) {
    $text = Get-Content -LiteralPath $LogPath -Raw
    return [pscustomobject]@{
        ReferenceDepthCm = Get-LastDouble $text "^Reference rod depth \(cm\):\s*(?<value>$numberPattern)\s*$"
        PerturbedDepthCm = Get-LastDouble $text "^Perturbed rod depth \(cm\):\s*(?<value>$numberPattern)\s*$"
        PairedReferenceK = Get-LastDouble $text "^Paired reference k:\s*(?<value>$numberPattern)\s*$"
        PerturbedK = Get-LastDouble $text "^Perturbed k, NPS \+ fission-source propagation.*:\s*(?<value>$numberPattern)\s*$"
        DeltaK = Get-LastDouble $text "^Total delta k:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$"
        DeltaKStderr = Get-LastDouble $text "^Total delta k:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$" 'stderr'
        DeltaRhoPcm = Get-LastDouble $text "^Total delta reactivity \(pcm\):\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$"
        DeltaRhoStderrPcm = Get-LastDouble $text "^Total delta reactivity \(pcm\):\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$" 'stderr'
        LastSourceIncrement = Get-LastDouble $text "^Last propagation increment,.*:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$"
        LastSourceIncrementStderr = Get-LastDouble $text "^Last propagation increment,.*:\s*(?<value>$numberPattern)\s+\+/-\s+(?<stderr>$numberPattern)\s*$" 'stderr'
        SupportDiagnosticCount = Get-LastInt64 $text '^Support/non-finite diagnostic count:\s*(?<value>[0-9]+)\s*$'
        UsableActiveCycles = Get-LastInt64 $text '^Usable active cycles:\s*(?<value>[0-9]+)\s*$'
    }
}

function Show-ProgressLine([string]$Line) {
    if ($Line -match '^Cycle\s+(?<cycle>[0-9]+),') {
        $cycle = [int]$matches.cycle
        if ($progressEveryCycles -le 1 -or ($cycle % $progressEveryCycles) -eq 0) {
            Write-Host $Line
        }
        return
    }
    if ($Line -match '^(Correlated sampling enabled:|Total Time:|k_mean:|Total delta k:|Total delta reactivity|Last propagation increment|Support/non-finite diagnostic count:|WARNING:|NOTICE:)') {
        Write-Host $Line
    }
}

function Invoke-C5G7Run {
    param(
        [double]$Depth,
        [AllowNull()][Nullable[double]]$CsDepth,
        [UInt64]$Seed,
        [string]$LogPath
    )

    $programArgs = @(
        '--critical-search-binary',
        "--rod-depth=$(Format-Number $Depth 'G17')",
        "--particles=$particles",
        "--inactive=$inactive",
        "--active=$active",
        "--seed=$Seed"
    )
    if ($null -ne $CsDepth) {
        $programArgs += "--cs-depth=$(Format-Number ([double]$CsDepth) 'G17')"
        $programArgs += "--cs-generations=$csGenerations"
    }

    $utf8NoBom = [Text.UTF8Encoding]::new($false)
    $writer = [IO.StreamWriter]::new($LogPath, $false, $utf8NoBom)
    try {
        & $script:exePath @programArgs 2>&1 | ForEach-Object {
            $line = $_.ToString()
            $writer.WriteLine($line)
            Show-ProgressLine $line
        }
        $status = $LASTEXITCODE
    }
    finally {
        $writer.Dispose()
    }
    if ($status -ne 0) { throw "C5G7 exited with code $status. See $LogPath" }
}

function Get-OrRunDirect {
    param([double]$Depth, [UInt64]$Seed, [string]$LogPath, [string]$Label)
    if (Test-Path -LiteralPath $LogPath) {
        try {
            $result = Read-DirectResult $LogPath
            Write-Host "Reusing completed $Label log: $LogPath"
            return $result
        }
        catch {
            Write-Warning "Existing $Label log is incomplete; rerunning it."
        }
    }
    Write-Host "Running $Label at d=$(Format-Number $Depth) cm, seed=$Seed"
    Invoke-C5G7Run -Depth $Depth -CsDepth $null -Seed $Seed -LogPath $LogPath
    return Read-DirectResult $LogPath
}

function Get-OrRunCs {
    param([double]$Depth, [double]$TargetDepth, [UInt64]$Seed, [string]$LogPath)
    if (Test-Path -LiteralPath $LogPath) {
        try {
            $result = Read-CsResult $LogPath
            Write-Host "Reusing completed CS log: $LogPath"
            return $result
        }
        catch {
            Write-Warning 'Existing CS log is incomplete; rerunning it.'
        }
    }
    Write-Host "Running CS d=$(Format-Number $Depth) -> $(Format-Number $TargetDepth) cm, seed=$Seed"
    Invoke-C5G7Run -Depth $Depth -CsDepth $TargetDepth -Seed $Seed -LogPath $LogPath
    return Read-CsResult $LogPath
}

function Get-SignAssessment([double]$Mean, [double]$Stderr) {
    $lower = $Mean - $confidenceSigma * $Stderr
    $upper = $Mean + $confidenceSigma * $Stderr
    if ($expectedDeltaKSign -lt 0) {
        if ($upper -lt 0) { return 'RESOLVED_EXPECTED' }
        if ($lower -gt 0) { return 'RESOLVED_WRONG' }
    }
    else {
        if ($lower -gt 0) { return 'RESOLVED_EXPECTED' }
        if ($upper -lt 0) { return 'RESOLVED_WRONG' }
    }
    return 'UNRESOLVED'
}

function Get-AgreementLabel([double]$Z) {
    $absoluteZ = [Math]::Abs($Z)
    if ($absoluteZ -le 1.0) { return 'AGREE_1SIGMA' }
    if ($absoluteZ -le 2.0) { return 'AGREE_2SIGMA' }
    return 'DISAGREE_GT_2SIGMA'
}

function Write-ValidationSvg {
    param([object[]]$Rows, [string]$Path)
    if ($Rows.Count -eq 0) { return }

    $width = 1100.0
    $height = 760.0
    $left = 105.0
    $right = 45.0
    $top1 = 75.0
    $bottom1 = 430.0
    $top2 = 535.0
    $bottom2 = 690.0
    $plotWidth = $width - $left - $right
    $xMin = 0.0
    $xMax = [double](($Rows | Measure-Object -Property ProbeStepCm -Maximum).Maximum)
    if ($xMax -le 0) { $xMax = 1.0 }

    $yValues = [Collections.Generic.List[double]]::new()
    foreach ($row in $Rows) {
        [void]$yValues.Add([double]$row.DirectDeltaK - [double]$row.DirectDeltaKStderr)
        [void]$yValues.Add([double]$row.DirectDeltaK + [double]$row.DirectDeltaKStderr)
        [void]$yValues.Add([double]$row.CsDeltaK - [double]$row.CsDeltaKStderr)
        [void]$yValues.Add([double]$row.CsDeltaK + [double]$row.CsDeltaKStderr)
        [void]$yValues.Add(0.0)
    }
    $yMin = [double](($yValues | Measure-Object -Minimum).Minimum)
    $yMax = [double](($yValues | Measure-Object -Maximum).Maximum)
    $yPad = [Math]::Max(1.0e-6, 0.08 * ($yMax - $yMin))
    $yMin -= $yPad
    $yMax += $yPad

    function XPixel([double]$x) { return $left + ($x - $xMin) / ($xMax - $xMin) * $plotWidth }
    function Y1Pixel([double]$y) { return $bottom1 - ($y - $yMin) / ($yMax - $yMin) * ($bottom1 - $top1) }
    function Y2Pixel([double]$z) {
        $zLimited = [Math]::Max(-3.0, [Math]::Min(3.0, $z))
        return $bottom2 - ($zLimited + 3.0) / 6.0 * ($bottom2 - $top2)
    }

    $svg = [Text.StringBuilder]::new()
    [void]$svg.AppendLine('<?xml version="1.0" encoding="UTF-8"?>')
    [void]$svg.AppendLine("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$([int]$width)`" height=`"$([int]$height)`" viewBox=`"0 0 $([int]$width) $([int]$height)`">")
    [void]$svg.AppendLine('<rect width="100%" height="100%" fill="white"/>')
    [void]$svg.AppendLine('<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#1f2937}.title{font-size:24px;font-weight:700}.label{font-size:16px}.tick{font-size:13px;fill:#4b5563}.grid{stroke:#d1d5db;stroke-width:1}.axis{stroke:#111827;stroke-width:1.5}.direct{stroke:#1769aa;fill:#1769aa}.cs{stroke:#d97706;fill:#d97706}</style>')
    [void]$svg.AppendLine('<text x="105" y="36" class="title">CS prediction versus independent direct Monte Carlo</text>')

    for ($i = 0; $i -le 5; $i++) {
        $xValue = $xMax * $i / 5.0
        $x = XPixel $xValue
        [void]$svg.AppendLine("<line x1=`"$x`" y1=`"$top1`" x2=`"$x`" y2=`"$bottom2`" class=`"grid`"/>")
        [void]$svg.AppendLine("<text x=`"$x`" y=`"$($bottom2 + 25)`" text-anchor=`"middle`" class=`"tick`">$(Format-Number $xValue '0.####')</text>")
    }
    for ($i = 0; $i -le 5; $i++) {
        $yValue = $yMin + ($yMax - $yMin) * $i / 5.0
        $y = Y1Pixel $yValue
        [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$y`" x2=`"$($width - $right)`" y2=`"$y`" class=`"grid`"/>")
        [void]$svg.AppendLine("<text x=`"$($left - 12)`" y=`"$($y + 5)`" text-anchor=`"end`" class=`"tick`">$(Format-Number $yValue '0.0#####')</text>")
    }
    foreach ($zLine in @(-2.0, -1.0, 0.0, 1.0, 2.0)) {
        $y = Y2Pixel $zLine
        $dash = if ([Math]::Abs($zLine) -eq 2.0) { '6,5' } elseif ([Math]::Abs($zLine) -eq 1.0) { '3,4' } else { '' }
        [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$y`" x2=`"$($width - $right)`" y2=`"$y`" stroke=`"#9ca3af`" stroke-dasharray=`"$dash`"/>")
        [void]$svg.AppendLine("<text x=`"$($left - 12)`" y=`"$($y + 5)`" text-anchor=`"end`" class=`"tick`">$(Format-Number $zLine '0')</text>")
    }
    [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$top1`" x2=`"$left`" y2=`"$bottom1`" class=`"axis`"/>")
    [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$bottom1`" x2=`"$($width - $right)`" y2=`"$bottom1`" class=`"axis`"/>")
    [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$top2`" x2=`"$left`" y2=`"$bottom2`" class=`"axis`"/>")
    [void]$svg.AppendLine("<line x1=`"$left`" y1=`"$bottom2`" x2=`"$($width - $right)`" y2=`"$bottom2`" class=`"axis`"/>")
    [void]$svg.AppendLine("<text x=`"$([int]($width / 2))`" y=`"$($height - 18)`" text-anchor=`"middle`" class=`"label`">Additional rod insertion, delta z (cm)</text>")
    [void]$svg.AppendLine('<text x="24" y="255" transform="rotate(-90 24 255)" text-anchor="middle" class="label">Delta k</text>')
    [void]$svg.AppendLine('<text x="24" y="615" transform="rotate(-90 24 615)" text-anchor="middle" class="label">Agreement z-score</text>')
    [void]$svg.AppendLine('<text x="735" y="62" class="label" fill="#1769aa">● Direct MC difference (±1σ)</text>')
    [void]$svg.AppendLine('<text x="735" y="86" class="label" fill="#d97706">■ CS estimate (±1σ)</text>')

    foreach ($row in $Rows) {
        $x = XPixel ([double]$row.ProbeStepCm)
        $directY = Y1Pixel ([double]$row.DirectDeltaK)
        $directLow = Y1Pixel ([double]$row.DirectDeltaK - [double]$row.DirectDeltaKStderr)
        $directHigh = Y1Pixel ([double]$row.DirectDeltaK + [double]$row.DirectDeltaKStderr)
        $csY = Y1Pixel ([double]$row.CsDeltaK)
        $csLow = Y1Pixel ([double]$row.CsDeltaK - [double]$row.CsDeltaKStderr)
        $csHigh = Y1Pixel ([double]$row.CsDeltaK + [double]$row.CsDeltaKStderr)
        $offset = 5.0
        [void]$svg.AppendLine("<line x1=`"$($x - $offset)`" y1=`"$directLow`" x2=`"$($x - $offset)`" y2=`"$directHigh`" class=`"direct`"/>")
        [void]$svg.AppendLine("<line x1=`"$($x - $offset - 4)`" y1=`"$directLow`" x2=`"$($x - $offset + 4)`" y2=`"$directLow`" class=`"direct`"/>")
        [void]$svg.AppendLine("<line x1=`"$($x - $offset - 4)`" y1=`"$directHigh`" x2=`"$($x - $offset + 4)`" y2=`"$directHigh`" class=`"direct`"/>")
        [void]$svg.AppendLine("<circle cx=`"$($x - $offset)`" cy=`"$directY`" r=`"5`" class=`"direct`"/>")
        [void]$svg.AppendLine("<line x1=`"$($x + $offset)`" y1=`"$csLow`" x2=`"$($x + $offset)`" y2=`"$csHigh`" class=`"cs`"/>")
        [void]$svg.AppendLine("<line x1=`"$($x + $offset - 4)`" y1=`"$csLow`" x2=`"$($x + $offset + 4)`" y2=`"$csLow`" class=`"cs`"/>")
        [void]$svg.AppendLine("<line x1=`"$($x + $offset - 4)`" y1=`"$csHigh`" x2=`"$($x + $offset + 4)`" y2=`"$csHigh`" class=`"cs`"/>")
        [void]$svg.AppendLine("<rect x=`"$($x + $offset - 4.5)`" y=`"$($csY - 4.5)`" width=`"9`" height=`"9`" class=`"cs`"/>")
        $zY = Y2Pixel ([double]$row.DeltaKAgreementZ)
        $zColor = if ([Math]::Abs([double]$row.DeltaKAgreementZ) -le 2.0) { '#15803d' } else { '#b91c1c' }
        [void]$svg.AppendLine("<circle cx=`"$x`" cy=`"$zY`" r=`"5.5`" fill=`"$zColor`"/>")
    }

    [void]$svg.AppendLine('<text x="105" y="505" class="label">z = (CS delta-k − direct delta-k) / combined standard error; |z| ≤ 2 indicates statistical agreement.</text>')
    [void]$svg.AppendLine('</svg>')
    [IO.File]::WriteAllText($Path, $svg.ToString(), [Text.UTF8Encoding]::new($false))
}

$exePath = Get-ProjectPath $exeSetting
$buildScript = Join-Path $scriptDir 'build_windows_nvcc.ps1'
$forceBuildRequested = [bool]$config.ForceBuild -or $ForceBuild.IsPresent
$noBuildRequested = [bool]$config.NoBuild -or $NoBuild.IsPresent
$needsBuild = -not (Test-Path -LiteralPath $exePath)
if (-not $needsBuild -and -not $noBuildRequested) {
    $exeTime = (Get-Item -LiteralPath $exePath).LastWriteTimeUtc
    $sourceFiles = @(
        Get-ChildItem -Path (Join-Path $rootDir '*.cu') -File
        Get-ChildItem -Path (Join-Path $rootDir '*.cuh') -File
        Get-Item -LiteralPath $buildScript
    )
    $needsBuild = $null -ne ($sourceFiles | Where-Object { $_.LastWriteTimeUtc -gt $exeTime } | Select-Object -First 1)
}
if ($noBuildRequested) {
    if (-not (Test-Path -LiteralPath $exePath)) { throw "NoBuild requested, but executable is missing: $exePath" }
    Write-Host "Using existing executable: $exePath"
}
elseif ($forceBuildRequested -or $needsBuild) {
    $buildParams = @{
        NvccPath = [string]$config.NvccPath
        CudaArch = [string]$config.CudaArch
        Output = $exeSetting
        ThreadsPerBlock = [int]$config.ThreadsPerBlock
    }
    if (-not [bool]$config.DeviceLto) { $buildParams.NoDeviceLto = $true }
    if ([bool]$config.FastMath) { $buildParams.FastMath = $true }
    Write-Host 'Building local CUDA validation executable...'
    & $buildScript @buildParams
}
else {
    Write-Host "Build is up to date: $exePath"
}
if (-not (Test-Path -LiteralPath $exePath)) { throw "Executable was not produced: $exePath" }

if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $runTag = Get-Date -Format 'yyyyMMdd_HHmmss'
    $suffix = if ($SmokeTest) { '_smoke' } else { '' }
    $outputPath = Join-Path $rootDir "results\cs_prediction_validation_local_${runTag}${suffix}"
}
else {
    $outputPath = Get-ProjectPath $OutputDir
}
New-Item -ItemType Directory -Force -Path $outputPath | Out-Null

Write-Host ''
Write-Host 'CS prediction validation against independent direct Monte Carlo'
Write-Host "  reference depth: $(Format-Number $referenceDepth) cm"
Write-Host "  probe steps: $(@($probeSteps | ForEach-Object { Format-Number $_ }) -join ', ') cm"
Write-Host "  particles/cycles: $particles / $inactive inactive + $active active"
Write-Host "  CS propagation level: L$csGenerations"
Write-Host "  all seeds are independent and explicitly forced odd"
Write-Host "  output: $outputPath"

$baseLog = Join-Path $outputPath ("direct_reference_d{0}.log" -f ((Format-Number $referenceDepth 'F8').Replace('.', 'p')))
$baseRunSeed = Get-OddSeed $baseSeed
$directReference = Get-OrRunDirect -Depth $referenceDepth -Seed $baseRunSeed -LogPath $baseLog -Label 'direct reference'

$rows = [Collections.Generic.List[object]]::new()
for ($index = 0; $index -lt $probeSteps.Count; $index++) {
    $step = [double]$probeSteps[$index]
    $targetDepth = $referenceDepth + $step
    $targetSeed = Get-OddSeed ($baseSeed + [UInt64](2 * (2 * $index + 1)))
    $csSeed = Get-OddSeed ($baseSeed + [UInt64](2 * (2 * $index + 2)))
    $stepTag = (Format-Number $step 'F8').Replace('.', 'p')
    $targetLog = Join-Path $outputPath ("direct_target_step_{0}.log" -f $stepTag)
    $csLog = Join-Path $outputPath ("cs_step_{0}.log" -f $stepTag)

    Write-Host ''
    Write-Host "Validation point $($index + 1)/$($probeSteps.Count): delta-z=$(Format-Number $step) cm"
    $directTarget = Get-OrRunDirect -Depth $targetDepth -Seed $targetSeed -LogPath $targetLog -Label 'direct target'
    $cs = Get-OrRunCs -Depth $referenceDepth -TargetDepth $targetDepth -Seed $csSeed -LogPath $csLog

    $directDeltaK = $directTarget.K - $directReference.K
    $directDeltaKStderr = [Math]::Sqrt(
        $directReference.KStderr * $directReference.KStderr +
        $directTarget.KStderr * $directTarget.KStderr)
    $deltaKDifference = $cs.DeltaK - $directDeltaK
    $deltaKAgreementStderr = [Math]::Sqrt(
        $cs.DeltaKStderr * $cs.DeltaKStderr +
        $directDeltaKStderr * $directDeltaKStderr)
    $deltaKAgreementZ = if ($deltaKAgreementStderr -gt 0) {
        $deltaKDifference / $deltaKAgreementStderr
    }
    else { [double]::NaN }

    $predictedTargetK = $directReference.K + $cs.DeltaK
    $predictedTargetKStderr = [Math]::Sqrt(
        $directReference.KStderr * $directReference.KStderr +
        $cs.DeltaKStderr * $cs.DeltaKStderr)

    $rhoReference = 1.0 - 1.0 / $directReference.K
    $rhoTarget = 1.0 - 1.0 / $directTarget.K
    $directDeltaRhoPcm = ($rhoTarget - $rhoReference) * 1.0e5
    $directDeltaRhoStderrPcm = [Math]::Sqrt(
        [Math]::Pow($directReference.KStderr / ($directReference.K * $directReference.K), 2) +
        [Math]::Pow($directTarget.KStderr / ($directTarget.K * $directTarget.K), 2)) * 1.0e5
    $rhoDifferencePcm = $cs.DeltaRhoPcm - $directDeltaRhoPcm
    $rhoAgreementStderrPcm = [Math]::Sqrt(
        $cs.DeltaRhoStderrPcm * $cs.DeltaRhoStderrPcm +
        $directDeltaRhoStderrPcm * $directDeltaRhoStderrPcm)
    $rhoAgreementZ = if ($rhoAgreementStderrPcm -gt 0) {
        $rhoDifferencePcm / $rhoAgreementStderrPcm
    }
    else { [double]::NaN }

    $sourceTailPcm = [Math]::Abs($cs.LastSourceIncrement) * 1.0e5
    $sourceTailUpperPcm = ([Math]::Abs($cs.LastSourceIncrement) +
        $confidenceSigma * $cs.LastSourceIncrementStderr) * 1.0e5
    $sourceTailLimitPcm = [Math]::Max(
        $sourceTailAbsPcm,
        $sourceTailRelFraction * [Math]::Abs($cs.DeltaK) * 1.0e5)
    $sourceConverged = $sourceTailUpperPcm -le $sourceTailLimitPcm
    $agreementLabel = Get-AgreementLabel $deltaKAgreementZ
    $directSign = Get-SignAssessment $directDeltaK $directDeltaKStderr
    $csSign = Get-SignAssessment $cs.DeltaK $cs.DeltaKStderr

    $status = 'PASS'
    if ($cs.SupportDiagnosticCount -ne 0) { $status = 'INVALID_CS_DIAGNOSTIC' }
    elseif (-not $sourceConverged) { $status = 'SOURCE_TAIL_NOT_CONVERGED' }
    elseif ($directSign -eq 'RESOLVED_WRONG' -or $csSign -eq 'RESOLVED_WRONG') { $status = 'MONOTONICITY_FAIL' }
    elseif ([Math]::Abs($deltaKAgreementZ) -gt $confidenceSigma) { $status = 'DIRECT_CS_DISAGREE' }

    $row = [pscustomobject][ordered]@{
        ProbeStepCm = $step
        ReferenceDepthRequestedCm = $referenceDepth
        ReferenceDepthReportedCm = $directReference.ActualDepthCm
        TargetDepthRequestedCm = $targetDepth
        TargetDepthReportedCm = $directTarget.ActualDepthCm
        DirectReferenceSeed = $baseRunSeed
        DirectTargetSeed = $targetSeed
        CsSeed = $csSeed
        DirectReferenceK = $directReference.K
        DirectReferenceKStderr = $directReference.KStderr
        DirectTargetK = $directTarget.K
        DirectTargetKStderr = $directTarget.KStderr
        DirectDeltaK = $directDeltaK
        DirectDeltaKStderr = $directDeltaKStderr
        CsDeltaK = $cs.DeltaK
        CsDeltaKStderr = $cs.DeltaKStderr
        CsPredictedTargetK = $predictedTargetK
        CsPredictedTargetKStderr = $predictedTargetKStderr
        PredictedMinusDirectTargetK = $predictedTargetK - $directTarget.K
        DeltaKAgreementStderr = $deltaKAgreementStderr
        DeltaKAgreementZ = $deltaKAgreementZ
        DeltaKAgreement = $agreementLabel
        DirectSlopePerCm = $directDeltaK / $step
        DirectSlopeStderrPerCm = $directDeltaKStderr / $step
        CsSlopePerCm = $cs.DeltaK / $step
        CsSlopeStderrPerCm = $cs.DeltaKStderr / $step
        DirectDeltaRhoPcm = $directDeltaRhoPcm
        DirectDeltaRhoStderrPcm = $directDeltaRhoStderrPcm
        CsDeltaRhoPcm = $cs.DeltaRhoPcm
        CsDeltaRhoStderrPcm = $cs.DeltaRhoStderrPcm
        DeltaRhoAgreementZ = $rhoAgreementZ
        DirectSignAssessment = $directSign
        CsSignAssessment = $csSign
        LastSourceIncrementPcm = $sourceTailPcm
        SourceTailConfidenceUpperPcm = $sourceTailUpperPcm
        SourceTailLimitPcm = $sourceTailLimitPcm
        SourceTailConverged = $sourceConverged
        CsSupportDiagnosticCount = $cs.SupportDiagnosticCount
        Status = $status
        DirectReferenceLog = [IO.Path]::GetFileName($baseLog)
        DirectTargetLog = [IO.Path]::GetFileName($targetLog)
        CsLog = [IO.Path]::GetFileName($csLog)
    }
    [void]$rows.Add($row)
    Write-Host "  direct delta-k=$(Format-Number $directDeltaK) +/- $(Format-Number $directDeltaKStderr)"
    Write-Host "  CS delta-k=$(Format-Number $cs.DeltaK) +/- $(Format-Number $cs.DeltaKStderr)"
    Write-Host "  agreement z=$(Format-Number $deltaKAgreementZ 'F3'); direct sign=$directSign; status=$status"
}

$csvPath = Join-Path $outputPath 'validation.csv'
$rows | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8
$svgPath = Join-Path $outputPath 'cs_vs_direct_delta_k.svg'
Write-ValidationSvg -Rows $rows.ToArray() -Path $svgPath

$passCount = @($rows | Where-Object { $_.Status -eq 'PASS' }).Count
$summaryPath = Join-Path $outputPath 'summary.md'
$summary = [Collections.Generic.List[string]]::new()
[void]$summary.Add('# CS prediction validation against direct Monte Carlo')
[void]$summary.Add('')
[void]$summary.Add("- Reference depth: $(Format-Number $referenceDepth) cm")
[void]$summary.Add("- Histories per cycle: $particles")
[void]$summary.Add("- Cycles: $inactive inactive + $active active")
[void]$summary.Add("- CS source propagation: L$csGenerations")
[void]$summary.Add("- Independent odd seeds: direct reference $baseRunSeed; every target and CS run uses another seed")
[void]$summary.Add("- Overall: $passCount / $($rows.Count) probe steps passed all enabled checks")
[void]$summary.Add('')
[void]$summary.Add('The primary residual is `CS delta-k - [direct k(target) - direct k(reference)]`. Its z-score uses the quadrature sum of the CS delta-k standard error and both independent direct-k standard errors.')
[void]$summary.Add('')
[void]$summary.Add('| step (cm) | direct delta-k | CS delta-k | agreement z | direct sign | source tail | status |')
[void]$summary.Add('|---:|---:|---:|---:|---|---|---|')
foreach ($row in $rows) {
    $tailText = if ($row.SourceTailConverged) { 'pass' } else { 'fail' }
    [void]$summary.Add("| $(Format-Number $row.ProbeStepCm '0.#####') | $(Format-Number $row.DirectDeltaK '0.000000000') ± $(Format-Number $row.DirectDeltaKStderr '0.000000000') | $(Format-Number $row.CsDeltaK '0.000000000') ± $(Format-Number $row.CsDeltaKStderr '0.000000000') | $(Format-Number $row.DeltaKAgreementZ '0.00') | $($row.DirectSignAssessment) | $tailText | $($row.Status) |")
}
[void]$summary.Add('')
[void]$summary.Add('Interpretation: `UNRESOLVED` direct sign means that two independent direct runs cannot distinguish the small rod-worth signal from Monte Carlo noise at the selected confidence level. A positive observed direct delta-k in that case is not evidence that insertion physically increased reactivity. Increase histories or active cycles before interpreting the sign.')
[void]$summary.Add('')
[void]$summary.Add('Rows share the same direct reference calculation. Per-row z-scores are valid, but the row residuals are correlated and must not be treated as independent samples in a naive global chi-square test.')
[IO.File]::WriteAllLines($summaryPath, $summary, [Text.UTF8Encoding]::new($false))

$protocolPath = Join-Path $outputPath 'run_settings.txt'
$settings = @(
    "config=$resolvedConfig",
    "executable=$exePath",
    "reference_depth_cm=$(Format-Number $referenceDepth 'G17')",
    "probe_steps_cm=$(@($probeSteps | ForEach-Object { Format-Number $_ 'G17' }) -join ',')",
    "particles=$particles",
    "inactive=$inactive",
    "active=$active",
    "cs_generations=$csGenerations",
    "base_seed=$baseSeed",
    "confidence_sigma=$(Format-Number $confidenceSigma)",
    "source_tail_abs_pcm=$(Format-Number $sourceTailAbsPcm)",
    "source_tail_relative_fraction=$(Format-Number $sourceTailRelFraction)"
)
[IO.File]::WriteAllLines($protocolPath, $settings, [Text.UTF8Encoding]::new($false))

Write-Host ''
Write-Host "Validation complete: $passCount / $($rows.Count) rows passed."
Write-Host "  table:   $csvPath"
Write-Host "  summary: $summaryPath"
Write-Host "  graph:   $svgPath"

