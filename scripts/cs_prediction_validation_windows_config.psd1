@{
    # Validation states. Depth is measured from the top of the core [cm].
    ReferenceDepthCm       = 83.6388953145
    ProbeStepsCm           = @(0.03125, 0.0625, 0.125, 0.25, 0.5)
    ModeledUpperDepthCm    = 128.52

    # Local pilot statistics. For a final study, try 2,000,000 / 100 / 200.
    Particles              = 500000
    Inactive               = 100
    Active                 = 150
    CsGenerations          = 10
    BaseSeed               = 92239000101

    # Interpretation. An inserted rod is expected to reduce k.
    ConfidenceSigma        = 2.0
    ExpectedDeltaKSign     = -1

    # The L10-L9 source tail passes when its confidence upper bound is below
    # max(SourceTailAbsPcm, SourceTailRelFraction * |CS delta-k|).
    SourceTailAbsPcm       = 1.0
    SourceTailRelFraction  = 0.02

    # Console and local CUDA build.
    ProgressEveryCycles    = 10
    ExePath                = 'build\windows_nvcc\c5g7_cs_validation.exe'
    NvccPath               = 'nvcc'
    CudaArch               = 'sm_120'
    ThreadsPerBlock        = 32
    DeviceLto              = $true
    FastMath               = $true
    ForceBuild             = $false
    NoBuild                = $false
}
