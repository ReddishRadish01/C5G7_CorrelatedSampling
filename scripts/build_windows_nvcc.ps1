param(
    [string]$NvccPath = "nvcc",
    [string]$CudaArch = "sm_120",
    [string]$Output = "build\windows_nvcc\c5g7.exe",
    [ValidateSet(32, 64, 128, 256, 512, 1024)]
    [int]$ThreadsPerBlock = 32,
    [string[]]$ExtraNvccFlags = @(),
    [switch]$Clean,
    [switch]$NoVsDevEnv,
    [switch]$AllowUnsupportedCompiler,
    [switch]$NoDeviceLto,
    [switch]$FastMath
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RootDir = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RootDir

function Find-VsDevCmd {
    $vswhereCandidates = @(
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\Installer\vswhere.exe"
    )

    foreach ($vswhere in $vswhereCandidates) {
        if (Test-Path -LiteralPath $vswhere) {
            $found = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -find "Common7\Tools\VsDevCmd.bat" | Select-Object -First 1
            if ($found) {
                return $found
            }
        }
    }

    $fallbacks = @(
        "${env:ProgramFiles}\Microsoft Visual Studio\2026\Community\Common7\Tools\VsDevCmd.bat",
        "${env:ProgramFiles}\Microsoft Visual Studio\2026\Professional\Common7\Tools\VsDevCmd.bat",
        "${env:ProgramFiles}\Microsoft Visual Studio\2026\Enterprise\Common7\Tools\VsDevCmd.bat",
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat",
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Professional\Common7\Tools\VsDevCmd.bat",
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\Enterprise\Common7\Tools\VsDevCmd.bat"
    )

    foreach ($candidate in $fallbacks) {
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
    }

    return ""
}

function Import-VsDevEnvironment {
    param([string]$VsDevCmd)

    Write-Host "Loading MSVC environment from:"
    Write-Host "  $VsDevCmd"

    $envLines = & cmd.exe /d /s /c "call `"$VsDevCmd`" -arch=x64 -host_arch=x64 >nul && set"
    if ($LASTEXITCODE -ne 0) {
        throw "VsDevCmd failed with exit code $LASTEXITCODE"
    }

    $pathLines = @($envLines | Where-Object { $_ -imatch "^PATH=" })
    $vsPathLine = $pathLines |
        Where-Object { $_ -imatch "\\VC\\Tools\\MSVC\\.*\\bin\\HostX64\\x64" } |
        Select-Object -First 1
    if (-not $vsPathLine) {
        $vsPathLine = $pathLines | Select-Object -First 1
    }
    if (-not $vsPathLine) {
        throw "VsDevCmd did not provide a PATH value."
    }

    $vsPathValue = $vsPathLine.Substring($vsPathLine.IndexOf("=") + 1)
    foreach ($line in $envLines) {
        if ($line -match "^(.*?)=(.*)$") {
            if ($matches[1] -ieq "PATH") {
                continue
            }
            Set-Item -Path "Env:$($matches[1])" -Value $matches[2]
        }
    }

    $env:PATH = $vsPathValue
    $hostCompiler = $vsPathValue.Split(";") |
        ForEach-Object { $_.Trim().Trim('"') } |
        Where-Object { $_ -and (Test-Path -LiteralPath (Join-Path $_ "cl.exe")) } |
        Select-Object -First 1
    if (-not $hostCompiler) {
        throw "VsDevCmd completed, but cl.exe was not found in its PATH."
    }

    return $hostCompiler
}

$hostCompilerDir = ""
if (-not $NoVsDevEnv) {
    $cl = Get-Command cl.exe -ErrorAction SilentlyContinue
    if ($cl) {
        $hostCompilerDir = Split-Path -Parent $cl.Source
    }
    else {
        $vsDevCmd = Find-VsDevCmd
        if ($vsDevCmd) {
            $hostCompilerDir = Import-VsDevEnvironment -VsDevCmd $vsDevCmd
        }
        else {
            throw "Could not find VsDevCmd.bat or cl.exe. Install the Visual Studio C++ build tools."
        }
    }
}
else {
    $cl = Get-Command cl.exe -ErrorAction SilentlyContinue
    if ($cl) {
        $hostCompilerDir = Split-Path -Parent $cl.Source
    }
}

if (-not $hostCompilerDir -or -not (Test-Path -LiteralPath (Join-Path $hostCompilerDir "cl.exe"))) {
    throw "cl.exe was not found. Run without -NoVsDevEnv, or launch from a Visual Studio Developer PowerShell."
}

$nvcc = Get-Command $NvccPath -ErrorAction SilentlyContinue
if (-not $nvcc -and $NvccPath -eq "nvcc" -and $env:CUDA_PATH) {
    $cudaNvcc = Join-Path $env:CUDA_PATH "bin\nvcc.exe"
    if (Test-Path -LiteralPath $cudaNvcc) {
        $nvcc = Get-Item -LiteralPath $cudaNvcc
    }
}
if (-not $nvcc) {
    throw "nvcc not found. Add CUDA bin to PATH, set CUDA_PATH, or pass -NvccPath."
}

$sources = @(
    "main.cu",
    "Core.cu",
    "Cycle.cu",
    "Neutron.cu",
    "RNG.cu",
    "Tally.cu",
    "XSParser.cu"
)

foreach ($source in $sources) {
    if (-not (Test-Path -LiteralPath $source)) {
        throw "Missing source file: $source"
    }
}

$outputPath = Join-Path $RootDir $Output
$outputDir = Split-Path -Parent $outputPath

if ($Clean -and (Test-Path -LiteralPath $outputDir)) {
    Remove-Item -LiteralPath $outputDir -Recurse -Force
}

New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$nvccArgs = @(
    "-O3",
    "-DNDEBUG",
    "-std=c++17",
    "-rdc=true",
    "-arch=$CudaArch",
    "-ccbin=$hostCompilerDir",
    "-DC5G7_THREADS_PER_BLOCK=$ThreadsPerBlock",
    "-split-compile=0",
    "-Xptxas",
    "-O3",
    "-Xcompiler=/O2",
    "-Xcompiler=/EHsc",
    "-Xcompiler=/wd4819"
)

if (-not $NoDeviceLto) {
    $nvccArgs += "-dlto"
}

if ($FastMath) {
    $nvccArgs += "--use_fast_math"
}

if ($AllowUnsupportedCompiler) {
    $nvccArgs += "-allow-unsupported-compiler"
}

$nvccArgs += $ExtraNvccFlags
$nvccArgs += $sources
$nvccArgs += @("-o", $outputPath)

Write-Host "Building with direct nvcc..."
Write-Host "  nvcc: $($nvcc.Source)"
Write-Host "  host compiler: $(Join-Path $hostCompilerDir 'cl.exe')"
Write-Host "  arch: $CudaArch"
Write-Host "  threads/block: $ThreadsPerBlock"
Write-Host "  device LTO: $(-not $NoDeviceLto)"
Write-Host "  fast math: $($FastMath.IsPresent)"
Write-Host "  output: $outputPath"

& $nvcc.Source @nvccArgs
if ($LASTEXITCODE -ne 0) {
    throw "nvcc failed with exit code $LASTEXITCODE"
}

Write-Host "Built $outputPath"
