#pragma once

#include "cudaHeader.cuh"

#include <cmath>
#include <cfloat>

// Logarithmic histograms cover 10^-12 through 10^12.  Bin 0 and the last
// bin collect underflow and overflow, respectively.
constexpr int CS_WEIGHT_HIST_INTERIOR_BINS = 480;
constexpr int CS_WEIGHT_HIST_BINS = CS_WEIGHT_HIST_INTERIOR_BINS + 2;
constexpr double CS_WEIGHT_HIST_LOG10_MIN = -12.0;
constexpr double CS_WEIGHT_HIST_LOG10_MAX = 12.0;
constexpr double CS_WEIGHT_HIST_LOG10_WIDTH =
	(CS_WEIGHT_HIST_LOG10_MAX - CS_WEIGHT_HIST_LOG10_MIN) /
	static_cast<double>(CS_WEIGHT_HIST_INTERIOR_BINS);
constexpr int CS_WEIGHT_HIST_SAMPLE_STRIDE = 64;

// All quantities are observed at fission-scoring collisions.  The maxima and
// comparison counters use every observation.  Histograms use a deterministic
// 1/64 sample of particle indices to keep diagnostic overhead small.
struct CsWeightDiagnostics {
	unsigned long long observations;
	unsigned long long invalidObservations;
	unsigned long long histogramObservations;
	unsigned long long npsAboveReference;
	unsigned long long sourceAboveReference;
	unsigned long long fissionAboveReference;
	unsigned long long scoreAboveReference;

	double maxReferenceWeight;
	double maxNpsWeight;
	double maxSourceWeight;
	double maxFissionWeight;
	double maxScoreWeight;
	double maxNpsRatio;
	double maxSourceRatio;
	double maxFissionRatio;
	double maxScoreRatio;
	double maxNpsIncrease;
	double maxSourceIncrease;
	double maxFissionIncrease;
	double maxScoreIncrease;
	double maxNpsDecrease;
	double maxSourceDecrease;
	double maxFissionDecrease;
	double maxScoreDecrease;

	unsigned long long referenceHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long npsHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long sourceHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long fissionHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long scoreHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long npsRatioHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long sourceRatioHistogram[CS_WEIGHT_HIST_BINS];
	unsigned long long scoreRatioHistogram[CS_WEIGHT_HIST_BINS];
};

D inline void csAtomicMaxDouble(double* address, double value) {
	if (!isfinite(value) || value < 0.0) return;
	unsigned long long* bits = reinterpret_cast<unsigned long long*>(address);
	unsigned long long old = *bits;
	while (__longlong_as_double(static_cast<long long>(old)) < value) {
		const unsigned long long assumed = old;
		old = atomicCAS(bits, assumed,
			static_cast<unsigned long long>(__double_as_longlong(value)));
		if (old == assumed) break;
	}
}

D inline int csWeightHistogramBinFromNaturalLog(double naturalLogValue) {
	if (!isfinite(naturalLogValue)) {
		return naturalLogValue > 0.0 ? CS_WEIGHT_HIST_BINS - 1 : 0;
	}
	constexpr double inverseLn10 = 0.43429448190325182765;
	const double log10Value = naturalLogValue * inverseLn10;
	if (log10Value < CS_WEIGHT_HIST_LOG10_MIN) return 0;
	if (log10Value >= CS_WEIGHT_HIST_LOG10_MAX) return CS_WEIGHT_HIST_BINS - 1;
	const int interior = static_cast<int>(
		(log10Value - CS_WEIGHT_HIST_LOG10_MIN) / CS_WEIGHT_HIST_LOG10_WIDTH);
	return 1 + interior;
}

D inline void recordCsWeightDiagnostics(
	CsWeightDiagnostics* diagnostics,
	bool histogramSample,
	double referenceWeight,
	double npsLogRatio,
	double sourceRatio,
	double fissionXsRatio) {
	if (diagnostics == nullptr || !(referenceWeight > 0.0) ||
		!isfinite(referenceWeight) || !isfinite(npsLogRatio) ||
		!(sourceRatio >= 0.0) || !isfinite(sourceRatio) ||
		!(fissionXsRatio >= 0.0) || !isfinite(fissionXsRatio)) {
		if (diagnostics != nullptr) atomicAdd(&diagnostics->invalidObservations, 1ULL);
		return;
	}

	const double referenceLogWeight = log(referenceWeight);
	const double sourceLogRatio = sourceRatio > 0.0
		? npsLogRatio + log(sourceRatio) : -DBL_MAX;
	const double fissionLogRatio = fissionXsRatio > 0.0
		? npsLogRatio + log(fissionXsRatio) : -DBL_MAX;
	const double scoreLogRatio = (sourceRatio > 0.0 && fissionXsRatio > 0.0)
		? npsLogRatio + log(sourceRatio) + log(fissionXsRatio) : -DBL_MAX;

	const double npsRatio = npsLogRatio < 700.0 ? exp(npsLogRatio) : DBL_MAX;
	const double sourceCombinedRatio = sourceLogRatio < 700.0 ? exp(sourceLogRatio) : DBL_MAX;
	const double fissionCombinedRatio = fissionLogRatio < 700.0 ? exp(fissionLogRatio) : DBL_MAX;
	const double scoreCombinedRatio = scoreLogRatio < 700.0 ? exp(scoreLogRatio) : DBL_MAX;
	const double npsWeight = referenceWeight * npsRatio;
	const double sourceWeight = referenceWeight * sourceCombinedRatio;
	const double fissionWeight = referenceWeight * fissionCombinedRatio;
	const double scoreWeight = referenceWeight * scoreCombinedRatio;

	atomicAdd(&diagnostics->observations, 1ULL);
	csAtomicMaxDouble(&diagnostics->maxReferenceWeight, referenceWeight);
	csAtomicMaxDouble(&diagnostics->maxNpsWeight, npsWeight);
	csAtomicMaxDouble(&diagnostics->maxSourceWeight, sourceWeight);
	csAtomicMaxDouble(&diagnostics->maxFissionWeight, fissionWeight);
	csAtomicMaxDouble(&diagnostics->maxScoreWeight, scoreWeight);
	csAtomicMaxDouble(&diagnostics->maxNpsRatio, npsRatio);
	csAtomicMaxDouble(&diagnostics->maxSourceRatio, sourceCombinedRatio);
	csAtomicMaxDouble(&diagnostics->maxFissionRatio, fissionCombinedRatio);
	csAtomicMaxDouble(&diagnostics->maxScoreRatio, scoreCombinedRatio);

	if (npsWeight > referenceWeight) {
		atomicAdd(&diagnostics->npsAboveReference, 1ULL);
		csAtomicMaxDouble(&diagnostics->maxNpsIncrease, npsWeight - referenceWeight);
	}
	else csAtomicMaxDouble(&diagnostics->maxNpsDecrease, referenceWeight - npsWeight);
	if (sourceWeight > referenceWeight) {
		atomicAdd(&diagnostics->sourceAboveReference, 1ULL);
		csAtomicMaxDouble(&diagnostics->maxSourceIncrease, sourceWeight - referenceWeight);
	}
	else csAtomicMaxDouble(&diagnostics->maxSourceDecrease, referenceWeight - sourceWeight);
	if (fissionWeight > referenceWeight) {
		atomicAdd(&diagnostics->fissionAboveReference, 1ULL);
		csAtomicMaxDouble(&diagnostics->maxFissionIncrease, fissionWeight - referenceWeight);
	}
	else csAtomicMaxDouble(&diagnostics->maxFissionDecrease, referenceWeight - fissionWeight);
	if (scoreWeight > referenceWeight) {
		atomicAdd(&diagnostics->scoreAboveReference, 1ULL);
		csAtomicMaxDouble(&diagnostics->maxScoreIncrease, scoreWeight - referenceWeight);
	}
	else csAtomicMaxDouble(&diagnostics->maxScoreDecrease, referenceWeight - scoreWeight);

	if (!histogramSample) return;
	atomicAdd(&diagnostics->histogramObservations, 1ULL);
	atomicAdd(&diagnostics->referenceHistogram[
		csWeightHistogramBinFromNaturalLog(referenceLogWeight)], 1ULL);
	atomicAdd(&diagnostics->npsHistogram[
		csWeightHistogramBinFromNaturalLog(referenceLogWeight + npsLogRatio)], 1ULL);
	atomicAdd(&diagnostics->sourceHistogram[
		csWeightHistogramBinFromNaturalLog(referenceLogWeight + sourceLogRatio)], 1ULL);
	atomicAdd(&diagnostics->fissionHistogram[
		csWeightHistogramBinFromNaturalLog(referenceLogWeight + fissionLogRatio)], 1ULL);
	atomicAdd(&diagnostics->scoreHistogram[
		csWeightHistogramBinFromNaturalLog(referenceLogWeight + scoreLogRatio)], 1ULL);
	atomicAdd(&diagnostics->npsRatioHistogram[
		csWeightHistogramBinFromNaturalLog(npsLogRatio)], 1ULL);
	atomicAdd(&diagnostics->sourceRatioHistogram[
		csWeightHistogramBinFromNaturalLog(sourceLogRatio)], 1ULL);
	atomicAdd(&diagnostics->scoreRatioHistogram[
		csWeightHistogramBinFromNaturalLog(scoreLogRatio)], 1ULL);
}
