@{
    # Monte Carlo calculation
    Particles                 = 2000000
    Inactive                  = 100
    Active                    = 200
    BaseSeed                  = 92235922383

    # Rod-depth search, measured from the top of the core [cm]
    InitialDepthCm            = 64.24
    LowerCm                   = 21.42
    UpperCm                   = 128.52
    TargetK                   = 1.0

    # Correlated-sampling derivative and safeguarded Newton update
    ProbeStepCm               = 0.25
    MinProbeStepCm            = 0.03125
    MaxProbeRetries           = 3
    CsGenerations             = 10
    Damping                   = 1.0
    MaxNewtonStepCm           = 30.0
    MaxIterations             = 8

    # Statistical and stopping criteria
    DepthToleranceCm          = 0.01
    KTolerancePcm             = 10.0
    MaxKStderrPcm             = 10.0
    MinAbsSlope               = 1.0e-5
    MinSlopeSigma             = 3.0
    ConfidenceSigma           = 2.0

    # The L10-L9 source tail is acceptable when its confidence upper bound is
    # below max(SourceTailAbsPcm, SourceTailRelFraction * |DeltaK|).
    SourceTailAbsPcm          = 1.0
    SourceTailRelFraction     = 0.02
    # A required h/2 check accepts slopes within this relative difference or
    # within ConfidenceSigma times their combined standard error.
    MaxSlopeRelativeDifference = 0.05
    RequireSourceConvergence  = $true
    FinalDirectCheck          = $true

    # Local CUDA build
    ExePath                   = 'build\windows_nvcc\c5g7_fast_math.exe'
    NvccPath                  = 'nvcc'
    CudaArch                  = 'sm_120'
    ThreadsPerBlock           = 32
    DeviceLto                 = $true
    FastMath                  = $true
    ForceBuild                = $false
    NoBuild                   = $false
}
