
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "GpuManager.cuh"
#include "CoreManager.cuh"
#include "XSManager.cuh"
#include "Interaction.cuh"
#include "NeutronBankManager.cuh"
#include "Tally.cuh"
#include "Debug.cuh"
#include "Cycle.cuh"

//#include "MCdrawer.h"

#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <array>
#include <cmath>
#include <vector>
#include <limits>

//#define TALLY

//#define XSRESULTDEBUGind

#define TALLYLASTCYCLE

#ifndef C5G7_THREADS_PER_BLOCK
#define C5G7_THREADS_PER_BLOCK 32
#endif

enum class RunCase {
	Unrodded,
	RoddedA,
	RoddedB,
	Reference2003,
	CriticalSearch_Binary,
	ARO,
	TEST
};

// Edit these defaults when running from Visual Studio without command arguments.
//constexpr RunCase DEFAULT_RUN_CASE = RunCase::Reference2003;
constexpr RunCase DEFAULT_RUN_CASE = RunCase::RoddedB;
constexpr double DEFAULT_ROD_INSERTION_DEPTH_CM = 0.0; // Edit this for Visual Studio runs without command arguments.
constexpr int DEFAULT_PARTICLE_COUNT = 500000;
constexpr int DEFAULT_INACTIVE_CYCLES = 200;
constexpr int DEFAULT_ACTIVE_CYCLES = 201;
constexpr unsigned long long DEFAULT_RNG_SEED = 92235922383ULL;
constexpr double DEFAULT_TALLY_Z_POSITION_CM = 0.63;

struct RunConfig {
	RunCase runCase = DEFAULT_RUN_CASE;
	double rodInsertionDepth = DEFAULT_ROD_INSERTION_DEPTH_CM;
	int initialNum = DEFAULT_PARTICLE_COUNT;
	int inactiveCycle = DEFAULT_INACTIVE_CYCLES;
	int activeCycle = DEFAULT_ACTIVE_CYCLES;
	unsigned long long seed = DEFAULT_RNG_SEED;
	double tallyZPositionCm = DEFAULT_TALLY_Z_POSITION_CM;
	bool csEnabled = false;
	double csPerturbedRodDepth = 0.0;
	int csPropagationGenerations = CS_MAX_PROPAGATION_GENERATIONS;
	bool csTrackWeights = false;
};

struct RunningStatistics {
	int count = 0;
	double mean = 0.0;
	double M2 = 0.0;

	void add(double value) {
		++count;
		const double delta = value - mean;
		mean += delta / static_cast<double>(count);
		M2 += delta * (value - mean);
	}

	double standardDeviation() const {
		return count >= 2 ? std::sqrt(M2 / static_cast<double>(count - 1)) : 0.0;
	}

	double standardError() const {
		return count >= 2 ? standardDeviation() / std::sqrt(static_cast<double>(count)) : 0.0;
	}
};

struct StandardErrorEstimate {
	double iidCycle = 0.0;
	double batchMeans = 0.0;
	double conservative = 0.0;
	int batchSize = 0;
	int batchCount = 0;
};

// Power-iteration cycles are generally autocorrelated.  The ordinary
// s/sqrt(N) error assumes independent cycles, so also estimate the error from
// contiguous non-overlapping batches.  Reporting the larger of the two keeps
// short/noisy runs from appearing more precise merely because one estimator
// happened to fluctuate downward.
static StandardErrorEstimate estimateStandardError(
	const std::vector<double>& samples,
	double iidCycleStandardError) {
	StandardErrorEstimate estimate{};
	estimate.iidCycle = iidCycleStandardError;
	estimate.batchMeans = iidCycleStandardError;
	estimate.conservative = iidCycleStandardError;

	const int sampleCount = static_cast<int>(samples.size());
	if (sampleCount < 4) {
		return estimate;
	}

	estimate.batchSize = std::max(2, static_cast<int>(std::floor(std::sqrt(
		static_cast<double>(sampleCount)))));
	estimate.batchCount = sampleCount / estimate.batchSize;
	if (estimate.batchCount < 2) {
		estimate.batchSize = 0;
		estimate.batchCount = 0;
		return estimate;
	}

	RunningStatistics batchStatistics{};
	for (int batch = 0; batch < estimate.batchCount; ++batch) {
		double batchSum = 0.0;
		const int first = batch * estimate.batchSize;
		for (int offset = 0; offset < estimate.batchSize; ++offset) {
			batchSum += samples[first + offset];
		}
		batchStatistics.add(batchSum / static_cast<double>(estimate.batchSize));
	}

	estimate.batchMeans = batchStatistics.standardError();
	estimate.conservative = std::max(estimate.iidCycle, estimate.batchMeans);
	return estimate;
}

static RunCase parseRunCaseName(std::string name) {
	std::transform(name.begin(), name.end(), name.begin(),
		[](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
	name.erase(std::remove_if(name.begin(), name.end(),
		[](char ch) { return ch == '-' || ch == '_'; }), name.end());

	if (name == "unrodded") {
		return RunCase::Unrodded;
	}
	if (name == "roddeda") {
		return RunCase::RoddedA;
	}
	if (name == "roddedb") {
		return RunCase::RoddedB;
	}
	if (name == "reference2003" || name == "reference") {
		return RunCase::Reference2003;
	}
	if (name == "criticalsearchbinary" || name == "critical") {
		return RunCase::CriticalSearch_Binary;
	}
	if (name == "aro") {
		return RunCase::ARO;
	}
	if (name == "test") {
		return RunCase::TEST;
	}

	throw std::invalid_argument(
		"Unknown runCase '" + name +
		"'. Use unrodded, rodded-a, rodded-b, reference2003, critical, aro, or test.");
}

static RunConfig parseRunConfig(int argc, char** argv) {
	RunConfig config{};
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "--case") {
			if (i + 1 >= argc) {
				throw std::invalid_argument("--case requires a runCase value.");
			}
			config.runCase = parseRunCaseName(argv[++i]);
		}
		else if (arg.rfind("--case=", 0) == 0) {
			config.runCase = parseRunCaseName(arg.substr(7));
		}
		else if (arg == "--rodded-a" || arg == "--roddedA") {
			config.runCase = RunCase::RoddedA;
		}
		else if (arg == "--rodded-b" || arg == "--roddedB") {
			config.runCase = RunCase::RoddedB;
		}
		else if (arg == "--reference2003") {
			config.runCase = RunCase::Reference2003;
			config.rodInsertionDepth = 0.0;
		}
		else if (arg == "--unrodded") {
			config.runCase = RunCase::Unrodded;
		}
		else if (arg == "--critical-search-binary" || arg == "--critical") {
			config.runCase = RunCase::CriticalSearch_Binary;
		}
		else if (arg == "--rod-depth" && i + 1 < argc) {
			config.rodInsertionDepth = std::atof(argv[++i]);
		}
		else if (arg.rfind("--rod-depth=", 0) == 0) {
			config.rodInsertionDepth = std::atof(arg.substr(12).c_str());
		}
		else if (arg == "--particles" && i + 1 < argc) {
			config.initialNum = std::atoi(argv[++i]);
		}
		else if (arg.rfind("--particles=", 0) == 0) {
			config.initialNum = std::atoi(arg.substr(12).c_str());
		}
		else if (arg == "--inactive" && i + 1 < argc) {
			config.inactiveCycle = std::atoi(argv[++i]);
		}
		else if (arg.rfind("--inactive=", 0) == 0) {
			config.inactiveCycle = std::atoi(arg.substr(11).c_str());
		}
		else if (arg == "--active" && i + 1 < argc) {
			config.activeCycle = std::atoi(argv[++i]);
		}
		else if (arg.rfind("--active=", 0) == 0) {
			config.activeCycle = std::atoi(arg.substr(9).c_str());
		}
		else if (arg == "--seed" && i + 1 < argc) {
			config.seed = std::strtoull(argv[++i], nullptr, 10);
		}
		else if (arg.rfind("--seed=", 0) == 0) {
			config.seed = std::strtoull(arg.substr(7).c_str(), nullptr, 10);
		}
		else if (arg == "--tally-z" && i + 1 < argc) {
			config.tallyZPositionCm = std::atof(argv[++i]);
		}
		else if (arg.rfind("--tally-z=", 0) == 0) {
			config.tallyZPositionCm = std::atof(arg.substr(10).c_str());
		}
		else if (arg == "--cs-depth" && i + 1 < argc) {
			config.csPerturbedRodDepth = std::atof(argv[++i]);
			config.csEnabled = true;
		}
		else if (arg.rfind("--cs-depth=", 0) == 0) {
			config.csPerturbedRodDepth = std::atof(arg.substr(11).c_str());
			config.csEnabled = true;
		}
		else if (arg == "--cs-generations" && i + 1 < argc) {
			config.csPropagationGenerations = std::atoi(argv[++i]);
		}
		else if (arg.rfind("--cs-generations=", 0) == 0) {
			config.csPropagationGenerations = std::atoi(arg.substr(17).c_str());
		}
		else if (arg == "--cs-track-weights") {
			config.csTrackWeights = true;
		}
	}
	config.initialNum = std::max(1, config.initialNum);
	config.inactiveCycle = std::max(0, config.inactiveCycle);
	config.activeCycle = std::max(1, config.activeCycle);
	config.csPropagationGenerations = std::max(1,
		std::min(CS_MAX_PROPAGATION_GENERATIONS, config.csPropagationGenerations));
	if (config.runCase == RunCase::Reference2003) {
		config.rodInsertionDepth = 0.0;
	}
	return config;
}

static CsWeightDiagnostics emptyCsWeightDiagnostics() {
	return CsWeightDiagnostics{};
}

static void mergeCsWeightDiagnostics(
	CsWeightDiagnostics& destination,
	const CsWeightDiagnostics& source) {
	destination.observations += source.observations;
	destination.invalidObservations += source.invalidObservations;
	destination.histogramObservations += source.histogramObservations;
	destination.npsAboveReference += source.npsAboveReference;
	destination.sourceAboveReference += source.sourceAboveReference;
	destination.fissionAboveReference += source.fissionAboveReference;
	destination.scoreAboveReference += source.scoreAboveReference;

	destination.maxReferenceWeight = std::max(destination.maxReferenceWeight, source.maxReferenceWeight);
	destination.maxNpsWeight = std::max(destination.maxNpsWeight, source.maxNpsWeight);
	destination.maxSourceWeight = std::max(destination.maxSourceWeight, source.maxSourceWeight);
	destination.maxFissionWeight = std::max(destination.maxFissionWeight, source.maxFissionWeight);
	destination.maxScoreWeight = std::max(destination.maxScoreWeight, source.maxScoreWeight);
	destination.maxNpsRatio = std::max(destination.maxNpsRatio, source.maxNpsRatio);
	destination.maxSourceRatio = std::max(destination.maxSourceRatio, source.maxSourceRatio);
	destination.maxFissionRatio = std::max(destination.maxFissionRatio, source.maxFissionRatio);
	destination.maxScoreRatio = std::max(destination.maxScoreRatio, source.maxScoreRatio);
	destination.maxNpsIncrease = std::max(destination.maxNpsIncrease, source.maxNpsIncrease);
	destination.maxSourceIncrease = std::max(destination.maxSourceIncrease, source.maxSourceIncrease);
	destination.maxFissionIncrease = std::max(destination.maxFissionIncrease, source.maxFissionIncrease);
	destination.maxScoreIncrease = std::max(destination.maxScoreIncrease, source.maxScoreIncrease);
	destination.maxNpsDecrease = std::max(destination.maxNpsDecrease, source.maxNpsDecrease);
	destination.maxSourceDecrease = std::max(destination.maxSourceDecrease, source.maxSourceDecrease);
	destination.maxFissionDecrease = std::max(destination.maxFissionDecrease, source.maxFissionDecrease);
	destination.maxScoreDecrease = std::max(destination.maxScoreDecrease, source.maxScoreDecrease);

	for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) {
		destination.referenceHistogram[bin] += source.referenceHistogram[bin];
		destination.npsHistogram[bin] += source.npsHistogram[bin];
		destination.sourceHistogram[bin] += source.sourceHistogram[bin];
		destination.fissionHistogram[bin] += source.fissionHistogram[bin];
		destination.scoreHistogram[bin] += source.scoreHistogram[bin];
		destination.npsRatioHistogram[bin] += source.npsRatioHistogram[bin];
		destination.sourceRatioHistogram[bin] += source.sourceRatioHistogram[bin];
		destination.scoreRatioHistogram[bin] += source.scoreRatioHistogram[bin];
	}
}

static double csHistogramPercentile(
	const unsigned long long* histogram,
	double probability) {
	unsigned long long total = 0ULL;
	for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) total += histogram[bin];
	if (total == 0ULL) return std::numeric_limits<double>::quiet_NaN();
	const unsigned long long target = static_cast<unsigned long long>(
		std::ceil(probability * static_cast<double>(total)));
	unsigned long long cumulative = 0ULL;
	for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) {
		cumulative += histogram[bin];
		if (cumulative >= std::max(1ULL, target)) {
			if (bin == 0) return std::pow(10.0, CS_WEIGHT_HIST_LOG10_MIN);
			if (bin == CS_WEIGHT_HIST_BINS - 1) return std::numeric_limits<double>::infinity();
			const double log10Midpoint = CS_WEIGHT_HIST_LOG10_MIN +
				(static_cast<double>(bin) - 0.5) * CS_WEIGHT_HIST_LOG10_WIDTH;
			return std::pow(10.0, log10Midpoint);
		}
	}
	return std::numeric_limits<double>::quiet_NaN();
}

static double csHistogramPercentileBounded(
	const unsigned long long* histogram,
	double probability,
	double exactMaximum) {
	return std::min(csHistogramPercentile(histogram, probability), exactMaximum);
}

static double csFraction(unsigned long long numerator, unsigned long long denominator) {
	return denominator > 0ULL
		? static_cast<double>(numerator) / static_cast<double>(denominator)
		: 0.0;
}

static unsigned long long csHistogramTotal(const unsigned long long* histogram) {
	unsigned long long total = 0ULL;
	for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) total += histogram[bin];
	return total;
}

static double csHistogramFractionAbove(
	const unsigned long long* histogram,
	double threshold) {
	const unsigned long long total = csHistogramTotal(histogram);
	if (total == 0ULL || !(threshold > 0.0)) return 0.0;
	const double thresholdLog10 = std::log10(threshold);
	unsigned long long above = 0ULL;
	for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) {
		double midpointLog10 = 0.0;
		if (bin == 0) midpointLog10 = CS_WEIGHT_HIST_LOG10_MIN;
		else if (bin == CS_WEIGHT_HIST_BINS - 1) midpointLog10 = CS_WEIGHT_HIST_LOG10_MAX;
		else midpointLog10 = CS_WEIGHT_HIST_LOG10_MIN +
			(static_cast<double>(bin) - 0.5) * CS_WEIGHT_HIST_LOG10_WIDTH;
		if (midpointLog10 > thresholdLog10) above += histogram[bin];
	}
	return csFraction(above, total);
}

static double actualTopRodDepthCm(const C5G7Geometry& core) {
	double actualDepth = 0.0;

	for (int asmIdx = 0; asmIdx < core.assemblyNo; ++asmIdx) {
		const Assembly& assembly = core.assembly[asmIdx];
		if (assembly.pinCells == nullptr || assembly.xNum <= 0 || assembly.yNum <= 0 || assembly.zNum <= 0 || assembly.length.z <= 0.0) {
			continue;
		}

		const double cellHeight = assembly.length.z / static_cast<double>(assembly.zNum);
		double assemblyRodDepth = 0.0;
		for (int z = assembly.zNum - 1; z >= 0; --z) {
			double layerRodThickness = 0.0;
			for (int y = 0; y < assembly.yNum; ++y) {
				for (int x = 0; x < assembly.xNum; ++x) {
					const int index = pincellStorageIndex(x, y, z, assembly.xNum, assembly.yNum, assembly.xyToMortonRank);
					const Pincell& pinCell = assembly.pinCells[index];
					double pinRodThickness = 0.0;
					if (pinCell.hasAxialMeatSplit && pinCell.axialUpperMeatType == MatType::CR) {
						const double splitFraction = std::max(0.0, std::min(1.0, pinCell.axialMeatSplitFraction));
						pinRodThickness = (1.0 - splitFraction) * cellHeight;
					}
					else if (pinCell.meatType == MatType::CR) {
						pinRodThickness = cellHeight;
					}
					layerRodThickness = std::max(layerRodThickness, pinRodThickness);
				}
			}

			if (layerRodThickness <= 0.0) {
				break;
			}
			assemblyRodDepth += layerRodThickness;
			if (layerRodThickness < cellHeight) {
				break;
			}
		}

		actualDepth = std::max(actualDepth, assemblyRodDepth);
	}

	return actualDepth;
}

int main(int argc, char** argv) {
	cudaDeviceReset();
	RunConfig runConfig = parseRunConfig(argc, argv);
	int initialNum = runConfig.initialNum;
	int bankSize = initialNum * 2;

	
	int numCycle = 0;
	int inactiveCycle = runConfig.inactiveCycle;
	int activeCycle = runConfig.activeCycle;
	numCycle = inactiveCycle + activeCycle;
	int iterLimit = 10000;

	int tallyFetchCycleSpec = 200;

	constexpr int threadPerBlock = C5G7_THREADS_PER_BLOCK;
	int blockPerDim = (bankSize + threadPerBlock - 1) / threadPerBlock;

	double h_multK = 0.0;
	//double unroddedK = 1.14308;
	double unroddedK = 1.14308;
	double aroK = 1.14390;
	double reference2003K = 1.183810;
	double roddedA = 1.12806;
	double roddedB = 1.07777;
	double referenceK = 1.0;
	if (runConfig.runCase == RunCase::Unrodded) {
		referenceK = unroddedK;
	}
	else if (runConfig.runCase == RunCase::Reference2003) {
		referenceK = reference2003K;
	}
	else if (runConfig.runCase == RunCase::RoddedA) {
		referenceK = roddedA;
	}
	else if (runConfig.runCase == RunCase::RoddedB) {
		referenceK = roddedB;
	}
	else if (runConfig.runCase == RunCase::ARO || runConfig.runCase == RunCase::CriticalSearch_Binary) {
		referenceK = aroK;
	}
	else {
		referenceK = 1.0;
	}
	h_multK = referenceK;
	

	double* d_multK = nullptr;
	cudaMalloc(&d_multK, sizeof(double));
	cudaMemcpy(d_multK, &h_multK, sizeof(double), cudaMemcpyHostToDevice);

	const unsigned long long requestedSeedNo = runConfig.seed;
	const unsigned long long seedNo = (requestedSeedNo & 0xFFFFFFFFFFFFULL) | 1ULL;
	if (seedNo != requestedSeedNo) {
		std::cout << "Requested RNG seed: " << requestedSeedNo << "\n";
		std::cout << "RNG seed adjusted to odd 48-bit state: " << seedNo << "\n";
	}
	std::cout << "RNG seed: " << seedNo << "\n";

	GnuAMCM h_RNG(seedNo);
	unsigned long long* h_SeedArr = new unsigned long long[bankSize];
	for (int i = 0; i < bankSize; i++) {
		// Each GPU history starts from an odd 48-bit seed. Subsequent saved
		// states may be even or odd and must not be parity-adjusted again.
		h_SeedArr[i] = ((h_RNG.gen() + i) & 0xFFFFFFFFFFFFULL) | 1ULL;
	}
	unsigned long long* d_SeedArr = nullptr;
	cudaMalloc(&d_SeedArr, bankSize * sizeof(unsigned long long));
	cudaMemcpy(d_SeedArr, h_SeedArr, bankSize * sizeof(unsigned long long), cudaMemcpyHostToDevice);



	std::vector<MatXS> XS;
	
	MatXS h_UO2XS("C5txt/UO2.txt", MatType::UO2);
	MatXS h_MOX4_3("C5txt/Mox4_3.txt", MatType::MOX4_3);
	MatXS h_MOX7_0("C5txt/Mox7_0.txt", MatType::MOX7_0);
	MatXS h_MOX8_7("C5txt/Mox8_7.txt", MatType::MOX8_7);
	MatXS h_FC("C5txt/FC.txt", MatType::FC);
	MatXS h_GT("C5txt/GT.txt", MatType::GT);
	MatXS h_Mod("C5txt/Mod.txt", MatType::MOD);
	MatXS h_CR("C5txt/CR.txt", MatType::CR);
	
	/*
	MatXS h_UO2XS("C5Unity/UO2.txt", MatType::UO2);
	MatXS h_MOX4_3("C5Unity/Mox4_3.txt", MatType::MOX4_3);
	MatXS h_MOX7_0("C5Unity/Mox7_0.txt", MatType::MOX7_0);
	MatXS h_MOX8_7("C5Unity/Mox8_7.txt", MatType::MOX8_7);
	MatXS h_FC("C5Unity/FC.txt", MatType::FC);
	MatXS h_GT("C5Unity/GT.txt", MatType::GT);
	MatXS h_Mod("C5Unity/Mod.txt", MatType::MOD);
	*/
#ifdef XSRESULTDEBUG
	std::cout << h_UO2XS.transXS[0] << " " << h_Mod.transXS[2] << "\n";
#endif

	XS.reserve(8);
	//XS.push_back(h_UO2XS); // uncomfortable
	//XS.push_back(std::move(h_UO2XS)); //h_UO2XS is deprecated
	XS.emplace_back(h_UO2XS);
	XS.emplace_back(h_MOX4_3);
	XS.emplace_back(h_MOX7_0);
	XS.emplace_back(h_MOX8_7);
	XS.emplace_back(h_FC);
	XS.emplace_back(h_GT);
	XS.emplace_back(h_Mod);
	XS.emplace_back(h_CR);

	XSLibrary h_XSLib{};
	MatXSFactory::initialize(h_XSLib, XS);
	XSLibrary* d_XSLib = nullptr;
	GPU_Manager::XSLibDeviceAllocator(&d_XSLib, h_XSLib);

	std::string fileIdentifier;
	C5G7Geometry h_Core{};
	if (runConfig.runCase == RunCase::Unrodded) {
		fileIdentifier = "unrodded";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_unrodded/C5G7CoreGeometry.txt",
			"Geometry_unrodded/UO2Geometry.txt",
			"Geometry_unrodded/MOXGeometry.txt",
			"Geometry_unrodded/CRMODGeometry.txt",
			0.0,
			-1.0);
		h_Core.zMinBoundary = BoundaryCondition::Reflective;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running unrodded case. Rod insertion depth from top: "
			<< runConfig.rodInsertionDepth << " cm\n";
	}
	else if (runConfig.runCase == RunCase::RoddedA) {
		fileIdentifier = "roddedA";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_roddedA/C5G7CoreGeometry.txt",
			"Geometry_roddedA/UO2Geometry.txt",
			"Geometry_roddedA/MOXGeometry.txt",
			"Geometry_roddedA/CRMODGeometry.txt",
			0.0,
			-1.0);
		h_Core.zMinBoundary = BoundaryCondition::Reflective;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running Rodded A case. Rod insertion depth from top: "
			<< runConfig.rodInsertionDepth << " cm\n";
	}
	else if (runConfig.runCase == RunCase::RoddedB) {
		fileIdentifier = "roddedB";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_roddedB/C5G7CoreGeometry.txt",
			"Geometry_roddedB/UO2Geometry.txt",
			"Geometry_roddedB/MOXGeometry.txt",
			"Geometry_roddedB/CRMODGeometry.txt",
			0.0,
			-1.0);
		h_Core.zMinBoundary = BoundaryCondition::Reflective;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running Rodded B case. Rod insertion depth from top: "
			<< runConfig.rodInsertionDepth << " cm\n";
	}
	else if (runConfig.runCase == RunCase::CriticalSearch_Binary) {
		fileIdentifier = "CriticalSearch_Binary";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_ARO/C5G7CoreGeometry.txt",
			"Geometry_ARO/UO2Geometry.txt",
			"Geometry_ARO/MOXGeometry.txt",
			"Geometry_ARO/CRMODGeometry.txt",
			0.0,
			runConfig.rodInsertionDepth);
		h_Core.zMinBoundary = BoundaryCondition::Vacuum;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running extension rod case. Rod insertion depth from top: "
			<< runConfig.rodInsertionDepth << " cm\n";
	}
	else if (runConfig.runCase == RunCase::Reference2003) {
		fileIdentifier = "reference2003";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_reference2003/C5G7CoreGeometry.txt",
			"Geometry_reference2003/UO2Geometry.txt",
			"Geometry_reference2003/MOXGeometry.txt",
			"",
			0.0,
			-1.0);
		h_Core.zMinBoundary = BoundaryCondition::Vacuum;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running reference2003 case. Rod insertion depth from top: 0 cm\n";
	}
	else if (runConfig.runCase == RunCase::ARO) {
		fileIdentifier = "ARO";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_ARO/C5G7CoreGeometry.txt",
			"Geometry_ARO/UO2Geometry.txt",
			"Geometry_ARO/MOXGeometry.txt",
			"Geometry_ARO/CRMODGeometry.txt",
			0.0,
			-1.0);
		h_Core.zMinBoundary = BoundaryCondition::Vacuum;
		h_Core.zMaxBoundary = BoundaryCondition::Vacuum;
		std::cout << "Running ARO case. Rod insertion depth from top: 0 cm\n";
	}
	else if (runConfig.runCase == RunCase::TEST) {
		fileIdentifier = "TEST";
		C5G7GeometryFactory::Initialize(
			h_Core,
			"Geometry_TEST/C5G7CoreGeometry.txt",
			"Geometry_TEST/UO2Geometry.txt",
			"Geometry_TEST/MOXGeometry.txt",
			"",
			0.0,
			-1.0);
		h_Core.xMinBoundary = BoundaryCondition::Reflective;
		h_Core.xMaxBoundary = BoundaryCondition::Reflective;
		h_Core.yMinBoundary = BoundaryCondition::Reflective;
		h_Core.yMaxBoundary = BoundaryCondition::Reflective;
		h_Core.zMinBoundary = BoundaryCondition::Reflective;
		h_Core.zMaxBoundary = BoundaryCondition::Reflective;
		std::cout << "Running TEST case. Rod insertion depth from top: 0 cm\n";
	}
	else {
		throw std::runtime_error("Unhandled run case in geometry initialization.");
	}

	std::cout << h_Core.xMin << " " << h_Core.x << " "
		<< h_Core.yMin << " " << h_Core.y << " "
		<< h_Core.zMin << " " << h_Core.z << "\n";
	double actualRodDepth = actualTopRodDepthCm(h_Core);
	if (runConfig.runCase == RunCase::Unrodded) {
		actualRodDepth = 0.0;
	}
	else if (runConfig.runCase == RunCase::CriticalSearch_Binary) {
		actualRodDepth = std::max(0.0, std::min(runConfig.rodInsertionDepth, h_Core.z - h_Core.zMin));
	}
	std::cout << "Requested rod depth: " << runConfig.rodInsertionDepth
		<< " cm, actual modeled rod depth: " << actualRodDepth << " cm\n";

	CorrelatedSamplingConfig csConfig{};
	if (runConfig.csEnabled) {
		if (runConfig.runCase != RunCase::CriticalSearch_Binary && runConfig.runCase != RunCase::ARO) {
			throw std::invalid_argument("--cs-depth is currently supported only with --case critical or --case aro.");
		}
		const double activeHeight = h_Core.z - h_Core.zMin;
		if (runConfig.csPerturbedRodDepth < actualRodDepth - 1.0e-10) {
			throw std::invalid_argument("--cs-depth must be greater than or equal to the reference --rod-depth. This implementation samples GT -> CR insertion only.");
		}
		if (runConfig.csPerturbedRodDepth > activeHeight + 1.0e-10) {
			throw std::invalid_argument("--cs-depth exceeds the modeled core height.");
		}
		csConfig.enabled = 1;
		csConfig.propagationGenerations = runConfig.csPropagationGenerations;
		csConfig.referenceRodDepth = actualRodDepth;
		csConfig.perturbedRodDepth = runConfig.csPerturbedRodDepth;
		csConfig.coreTopZ = h_Core.z;
		std::cout << "Correlated sampling enabled: reference depth " << csConfig.referenceRodDepth
			<< " cm -> perturbed depth " << csConfig.perturbedRodDepth
			<< " cm, source propagation generations " << csConfig.propagationGenerations << "\n";
	}
	Assembly* d_bufferAssembly = nullptr;
	std::vector<Pincell*> d_bufferPincellVec(h_Core.assemblyNo, nullptr);
	std::vector<int*> d_bufferMortonRankVec(h_Core.assemblyNo, nullptr);
	C5G7Geometry* d_Core = GPU_Manager::CoreDeviceAllocator(h_Core, d_bufferAssembly, d_bufferPincellVec, d_bufferMortonRankVec);

	double h_NeutronNum = static_cast<double>(initialNum);
	double* d_NeutronNum = nullptr;
	double h_NeutronWeightSum = 1.0 * initialNum;
	double initialWeight = 1.0 * initialNum;
	double* d_NeutronWeightSum = nullptr;
	double h_NeutronWeightModifier = 1.0;
	double* d_NeutronWeightModifier = nullptr;
	cudaMalloc(&(d_NeutronNum), sizeof(double));
	cudaMemcpy(d_NeutronNum, &h_NeutronNum, sizeof(double), cudaMemcpyHostToDevice);
	cudaMalloc(&(d_NeutronWeightSum), sizeof(double));
	cudaMemcpy(d_NeutronWeightSum, &h_NeutronWeightSum, sizeof(double), cudaMemcpyHostToDevice);
	cudaMalloc(&(d_NeutronWeightModifier), sizeof(double));
	cudaMemcpy(d_NeutronWeightModifier, &h_NeutronWeightModifier, sizeof(double), cudaMemcpyHostToDevice);


	TallyC5G7Geometry h_CoreTally(h_Core);
	TallyAssembly* d_bufferTallyAssembly = nullptr;
	std::vector<TallyPincell*> d_bufferTallyPincellVec(h_CoreTally.assemblyNo, nullptr);
	std::vector<int*> d_bufferTallyMortonRankVec(h_CoreTally.assemblyNo, nullptr);
	TallyC5G7Geometry* d_CoreTally = GPU_Manager::CoreTallyDeviceAllocator(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec, d_bufferTallyMortonRankVec);

#ifdef TALLYLASTCYCLE
	BenchmarkPinPowerCycleStatistics pinPowerCycleStatistics;
	std::vector<std::vector<double>> cumulativeActivePinTally(h_CoreTally.assemblyNo);
	std::vector<std::vector<double>> cumulativeActiveModTally(h_CoreTally.assemblyNo);
	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const int cellCount = h_CoreTally.tallyAssembly[a].totalPincellNo();
		cumulativeActivePinTally[a].assign(cellCount, 0.0);
		cumulativeActiveModTally[a].assign(cellCount, 0.0);
	}
#endif




	Neutron h_Neutron{ {4.0, 4.0, 200.0}, {0.0, 0.0, 0.0}, 1.0, 1.0 };
	Neutron* d_Neutron = nullptr;
	cudaMalloc(&d_Neutron, sizeof(Neutron));
	cudaMemcpy(d_Neutron, &h_Neutron, sizeof(Neutron), cudaMemcpyHostToDevice);

	NeutronBank h_Bank(bankSize, initialNum, seedNo);


	// fetch K to text file

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << "results/k_history_"
		<< fileIdentifier << "_"
		<< (tm.tm_year + 1900)
		<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
		<< std::setw(2) << std::setfill('0') << tm.tm_mday
		<< "_"
		<< std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec
		<< ".txt";

	std::ostringstream csOss;
	if (csConfig.enabled) {
		csOss << "results/cs_history_" << fileIdentifier << "_"
			<< (tm.tm_year + 1900)
			<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
			<< std::setw(2) << std::setfill('0') << tm.tm_mday << "_"
			<< std::setw(2) << std::setfill('0') << tm.tm_hour
			<< std::setw(2) << std::setfill('0') << tm.tm_min
			<< std::setw(2) << std::setfill('0') << tm.tm_sec << ".txt";
	}
	std::ostringstream csWeightCycleOss;
	std::ostringstream csWeightHistogramOss;
	std::ostringstream csWeightSummaryOss;
	if (csConfig.enabled && runConfig.csTrackWeights) {
		const std::string timestamp =
			std::to_string(tm.tm_year + 1900) +
			(tm.tm_mon + 1 < 10 ? "0" : "") + std::to_string(tm.tm_mon + 1) +
			(tm.tm_mday < 10 ? "0" : "") + std::to_string(tm.tm_mday) + "_" +
			(tm.tm_hour < 10 ? "0" : "") + std::to_string(tm.tm_hour) +
			(tm.tm_min < 10 ? "0" : "") + std::to_string(tm.tm_min) +
			(tm.tm_sec < 10 ? "0" : "") + std::to_string(tm.tm_sec);
		csWeightCycleOss << "results/cs_weight_cycle_" << fileIdentifier << "_" << timestamp << ".csv";
		csWeightHistogramOss << "results/cs_weight_histogram_" << fileIdentifier << "_" << timestamp << ".csv";
		csWeightSummaryOss << "results/cs_weight_summary_" << fileIdentifier << "_" << timestamp << ".txt";
	}

	std::ostringstream oss2;
	oss2 << "results/flux_Tally_"
		<< fileIdentifier << "_"
		<< (tm.tm_year + 1900)
		<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
		<< std::setw(2) << std::setfill('0') << tm.tm_mday
		<< "_"
		<< std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec
		<< ".txt";

#ifdef C5G7_DUMP_FULL_CORE_PINWISE_TALLY
	std::ostringstream oss3;
	oss3 << "results/pinwise_Tally_"
		<< fileIdentifier << "_"
		<< (tm.tm_year + 1900)
		<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
		<< std::setw(2) << std::setfill('0') << tm.tm_mday
		<< "_"
		<< std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec
		<< ".txt";
#endif

#ifdef TALLYLASTCYCLE
	std::ostringstream oss4;
	oss4 << "results/pin_power_cycle_statistics_"
		<< fileIdentifier << "_"
		<< (tm.tm_year + 1900)
		<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
		<< std::setw(2) << std::setfill('0') << tm.tm_mday
		<< "_"
		<< std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec
		<< ".txt";
#endif

	std::ofstream klog(oss.str(), std::ios::out);
	klog << std::fixed << std::setprecision(6);
	std::ofstream csLog;
	std::ofstream csWeightCycleLog;
	if (csConfig.enabled) {
		csLog.open(csOss.str(), std::ios::out);
		csLog << std::setprecision(12);
		csLog << "# cycle k_reference";
		for (int level = 0; level <= csConfig.propagationGenerations; ++level) {
			csLog << " k_cs_L" << level << " delta_k_L" << level;
		}
		csLog << " source_count";
		for (int level = 1; level <= csConfig.propagationGenerations; ++level) {
			csLog << " source_mean_L" << level;
		}
		csLog << "\n";
		if (runConfig.csTrackWeights) {
			csWeightCycleLog.open(csWeightCycleOss.str(), std::ios::out);
			csWeightCycleLog << std::setprecision(12);
			csWeightCycleLog
				<< "cycle,phase,observations,histogram_observations,invalid_observations,"
				<< "max_w,p99_w,p999_w,"
				<< "max_wstar_nps,p99_wstar_nps,p999_wstar_nps,max_wstar_nps_over_w,"
				<< "fraction_wstar_nps_gt_w,max_wstar_nps_minus_w,max_w_minus_wstar_nps,"
				<< "max_wstar_source,p99_wstar_source,p999_wstar_source,max_wstar_source_over_w,"
				<< "fraction_wstar_source_gt_w,max_wstar_source_minus_w,max_w_minus_wstar_source,"
				<< "max_wstar_fission,p99_wstar_fission,p999_wstar_fission,max_wstar_fission_over_w,"
				<< "max_wstar_score,p99_wstar_score,p999_wstar_score,max_wstar_score_over_w,"
				<< "fraction_wstar_score_gt_w,max_wstar_score_minus_w,max_w_minus_wstar_score\n";
		}
	}

	std::ofstream fluxTallyLog(oss2.str(), std::ios::out);
	fluxTallyLog << std::fixed << std::setprecision(6);

#ifdef C5G7_DUMP_FULL_CORE_PINWISE_TALLY
	std::ofstream pinwiseTallyLog(oss3.str(), std::ios::out);
	pinwiseTallyLog << std::fixed << std::setprecision(6);
#endif

#ifdef TALLYLASTCYCLE
	std::ofstream pinPowerCycleStatisticsLog(oss4.str(), std::ios::out);
#endif


	for (int i = 0; i < h_Bank.allocatableNeutronNum; i++) {
		if (i < h_Bank.neutronSize) {
			vec3 randPos = { h_RNG.uniform_open(h_Core.xMin, h_Core.x), h_RNG.uniform_open(h_Core.yMin, h_Core.y), h_RNG.uniform_open(h_Core.zMin, h_Core.z) };
			vec3 centerRandPos = {
				h_RNG.uniform_open(h_Core.xMin, h_Core.xMin + (h_Core.x - h_Core.xMin) / 4.0),
				h_RNG.uniform_open(h_Core.yMin, h_Core.yMin + (h_Core.y - h_Core.yMin) / 4.0),
				h_RNG.uniform_open(h_Core.zMin, h_Core.zMin + (h_Core.z - h_Core.zMin) / 4.0) };
			//h_Bank.neutrons[i] = Neutron(centerRandPos, vec3::randomUnit(h_RNG), static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
			h_Bank.neutrons[i] = Neutron(randPos, vec3::randomUnit(h_RNG), static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
			
		}
		else {
			h_Bank.neutrons[i] = Neutron();
		}
		//h_Bank.neutrons[i] = Neutron({ 0.1, 0.1, 0.1 }, { -1.0, 0.0, 0.0 }, static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		//h_Bank.neutrons[i] = Neutron({ 0.63, 0.63, 0.63 }, { 1.0, 0.0, 0.0 }, static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		//h_Bank.neutrons[i] = Neutron(randPos, vec3::randomUnit(h_RNG), 1, 1.0);
		h_Bank.addedNeutrons[i] = Neutron();
		//h_Bank.addedNeutrons[i].status = false;
		//h_Bank.neutrons[i].printInfo();
		//h_Bank.addedNeutrons[i].printInfo();
	}

	NeutronBank* d_Bank = nullptr;
	cudaMalloc(&d_Bank, sizeof(NeutronBank));
	Neutron* d_bufferNeutrons = nullptr; Neutron* d_bufferAddedNeutrons = nullptr;
	cudaMalloc(&d_bufferNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron));
	cudaMemcpy(d_bufferNeutrons, h_Bank.neutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron), cudaMemcpyHostToDevice);
	cudaMalloc(&d_bufferAddedNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron));
	cudaMemcpy(d_bufferAddedNeutrons, h_Bank.addedNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron), cudaMemcpyHostToDevice);
	NeutronBank tmp_Bank = h_Bank;
	tmp_Bank.neutrons = d_bufferNeutrons;
	tmp_Bank.addedNeutrons = d_bufferAddedNeutrons;
	cudaMemcpy(d_Bank, &tmp_Bank, sizeof(NeutronBank), cudaMemcpyHostToDevice);
	

	Neutron h_testNeutron{ {4.0, 4.0, 200.0}, {0.0, 0.0, 0.0}, 1.0, 1.0 };

	
	//for (int i = 0; i < 10; i++) {	Debug::fuelLayoutDebug(h_Core.assembly[i]);	}
	
	//CPUTest(num, &h_XSLib, &h_Core, &h_Bank, h_SeedArr);
	//GPUTest << <blockPerDim, threadPerBlock >> > (num, d_MatXS, d_Core, d_Bank, d_SeedArr);

	//double tempK = h_multK;
	double previousNumNeutron = h_Bank.getTotalNeutronNum();
	double meanK = 0.0;
	double M2 = 0.0;
	int errorCounter = 0;
	int cycleNum = 0;
	int activeCount = 0;
	std::vector<double> activeKSamples;
	activeKSamples.reserve(activeCycle);


	double h_fissionCount = initialNum;
	double* d_fissionCount = nullptr;
	cudaMalloc(&(d_fissionCount), sizeof(double));
	cudaMemcpy(d_fissionCount, &h_fissionCount, sizeof(double), cudaMemcpyHostToDevice);

	std::array<double, CS_MAX_PROPAGATION_GENERATIONS + 1> h_csFissionScores{};
	std::array<double, CS_MAX_PROPAGATION_GENERATIONS + 1> h_csSourceSums{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csKStatistics{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csDeltaKStatistics{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csDeltaRhoStatistics{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csSourceEffectStatistics{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csPairedNpsStatistics{};
	std::array<RunningStatistics, CS_MAX_PROPAGATION_GENERATIONS + 1> csGenerationIncrementStatistics{};
	std::array<std::vector<double>, CS_MAX_PROPAGATION_GENERATIONS + 1> csDeltaKSamples{};
	std::array<std::vector<double>, CS_MAX_PROPAGATION_GENERATIONS + 1> csDeltaRhoSamples{};
	std::array<std::vector<double>, CS_MAX_PROPAGATION_GENERATIONS + 1> csSourceEffectSamples{};
	std::array<std::vector<double>, CS_MAX_PROPAGATION_GENERATIONS + 1> csGenerationIncrementSamples{};
	for (int level = 0; level <= CS_MAX_PROPAGATION_GENERATIONS; ++level) {
		csDeltaKSamples[level].reserve(activeCycle);
		csDeltaRhoSamples[level].reserve(activeCycle);
		csSourceEffectSamples[level].reserve(activeCycle);
		csGenerationIncrementSamples[level].reserve(activeCycle);
	}
	double h_csSourceCount = 0.0;
	unsigned long long h_csSupportMismatch = 0ULL;
	double* d_csFissionScores = nullptr;
	double* d_csSourceSums = nullptr;
	double* d_csSourceCount = nullptr;
	unsigned long long* d_csSupportMismatch = nullptr;
	CsWeightDiagnostics h_csWeightCycle = emptyCsWeightDiagnostics();
	CsWeightDiagnostics h_csWeightActive = emptyCsWeightDiagnostics();
	CsWeightDiagnostics* d_csWeightDiagnostics = nullptr;
	if (csConfig.enabled) {
		cudaMalloc(&d_csFissionScores, sizeof(double) * h_csFissionScores.size());
		cudaMalloc(&d_csSourceSums, sizeof(double) * h_csSourceSums.size());
		cudaMalloc(&d_csSourceCount, sizeof(double));
		cudaMalloc(&d_csSupportMismatch, sizeof(unsigned long long));
		cudaMemset(d_csSupportMismatch, 0, sizeof(unsigned long long));
		if (runConfig.csTrackWeights) {
			cudaMalloc(&d_csWeightDiagnostics, sizeof(CsWeightDiagnostics));
			cudaMemcpy(d_csWeightDiagnostics, &h_csWeightCycle,
				sizeof(CsWeightDiagnostics), cudaMemcpyHostToDevice);
			std::cout << "CS weight tracking enabled at every fission-scoring collision; "
				<< "histograms sample 1/" << CS_WEIGHT_HIST_SAMPLE_STRIDE
				<< " particle indices.\n";
		}
	}
	double fissionCountBuffer = initialWeight * h_multK;
	double weightBuffer = initialWeight;
	auto t_start = std::chrono::steady_clock::now();


	for (int cycle = 0; cycle < numCycle; cycle++) {
		cycleNum++;

#ifdef TALLYLASTCYCLE
		// Keep one active cycle per in-memory statistical batch. The host
		// accumulates the final full-cell tally separately, so no per-cycle
		// tally files are required.
		ResetCoreTallyOnDevice(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
#endif

		double currentNumNeutron = h_Bank.getTotalNeutronNum();

		cudaMemset(d_NeutronNum, 0, sizeof(double));
		cudaMemset(d_NeutronWeightSum, 0, sizeof(double));
		fetchNeutronNumber << <blockPerDim, threadPerBlock >> > (d_Bank, d_NeutronNum, d_NeutronWeightSum);
		cudaMemcpy(&h_NeutronNum, d_NeutronNum, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_NeutronWeightSum, d_NeutronWeightSum, sizeof(double), cudaMemcpyDeviceToHost);

		//std::cout << " before normalization, NeutronNum: " << h_NeutronNum << ", Neutron Weight Sum: " << h_NeutronWeightSum << ", initialWeight: " << initialWeight << "\n";


		h_NeutronWeightModifier = 1.0 * (initialWeight / h_NeutronWeightSum);
		cudaMemcpy(d_NeutronWeightModifier, &h_NeutronWeightModifier, sizeof(double), cudaMemcpyHostToDevice);
		weightBuffer = h_NeutronWeightSum;

		globalNeutronWeightOffset << <blockPerDim, threadPerBlock >> > (d_Bank, d_NeutronWeightModifier);

		h_NeutronNum = 0.0; cudaMemcpy(d_NeutronNum, &h_NeutronNum, sizeof(double), cudaMemcpyHostToDevice);
		h_NeutronWeightSum = 0.0; cudaMemcpy(d_NeutronWeightSum, &h_NeutronWeightSum, sizeof(double), cudaMemcpyHostToDevice);

		fetchNeutronNumber << <blockPerDim, threadPerBlock >> > (d_Bank, d_NeutronNum, d_NeutronWeightSum);
		cudaMemcpy(&h_NeutronNum, d_NeutronNum, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_NeutronWeightSum, d_NeutronWeightSum, sizeof(double), cudaMemcpyDeviceToHost);

		std::cout << "B4 weight: " << weightBuffer << ", After weight normalization, NeutronNum: " << h_NeutronNum << ", Neutron Weight Sum: " << h_NeutronWeightSum << "\n";
		
		

		//fissionCountBuffer = h_fissionCount;
		h_fissionCount = 0.0; cudaMemcpy(d_fissionCount, &h_fissionCount, sizeof(double), cudaMemcpyHostToDevice); // same as cudaMemset but explicitly set h_fissionCount to 0.0
		if (csConfig.enabled) {
			cudaMemset(d_csFissionScores, 0, sizeof(double) * h_csFissionScores.size());
			cudaMemset(d_csSourceSums, 0, sizeof(double) * h_csSourceSums.size());
			cudaMemset(d_csSourceCount, 0, sizeof(double));
			if (d_csWeightDiagnostics != nullptr) {
				h_csWeightCycle = emptyCsWeightDiagnostics();
				cudaMemcpy(d_csWeightDiagnostics, &h_csWeightCycle,
					sizeof(CsWeightDiagnostics), cudaMemcpyHostToDevice);
			}
		}
		//cudaMemset(&(d_fissionCount), 0.0, sizeof(double));

		cycle_addedNeutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_CoreTally, d_XSLib, d_SeedArr, d_multK, true, iterLimit, d_fissionCount,
			csConfig, d_csFissionScores, d_csSourceSums, d_csSourceCount, d_csSupportMismatch,
			d_csWeightDiagnostics);
#ifdef _DEBUG
		CUDA_KERNEL_CHECK();
#endif
		addedNeutronPassResetter << <blockPerDim, threadPerBlock >> > (d_Bank);
#ifdef _DEBUG
		CUDA_KERNEL_CHECK();
#endif
		cycle_Neutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_CoreTally, d_XSLib, d_SeedArr, d_multK, false, iterLimit, d_fissionCount,
			csConfig, d_csFissionScores, d_csSourceSums, d_csSourceCount, d_csSupportMismatch,
			d_csWeightDiagnostics);
#ifdef _DEBUG
		CUDA_KERNEL_CHECK();
#endif

		if (csConfig.enabled) {
			cudaMemcpy(h_csFissionScores.data(), d_csFissionScores,
				sizeof(double) * h_csFissionScores.size(), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_csSourceSums.data(), d_csSourceSums,
				sizeof(double) * h_csSourceSums.size(), cudaMemcpyDeviceToHost);
			cudaMemcpy(&h_csSourceCount, d_csSourceCount, sizeof(double), cudaMemcpyDeviceToHost);
			if (d_csWeightDiagnostics != nullptr) {
				cudaMemcpy(&h_csWeightCycle, d_csWeightDiagnostics,
					sizeof(CsWeightDiagnostics), cudaMemcpyDeviceToHost);
				if (cycle >= inactiveCycle) {
					mergeCsWeightDiagnostics(h_csWeightActive, h_csWeightCycle);
				}
				csWeightCycleLog
					<< (cycle + 1) << "," << (cycle >= inactiveCycle ? "active" : "inactive") << ","
					<< h_csWeightCycle.observations << ","
					<< h_csWeightCycle.histogramObservations << ","
					<< h_csWeightCycle.invalidObservations << ","
					<< h_csWeightCycle.maxReferenceWeight << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.referenceHistogram, 0.99, h_csWeightCycle.maxReferenceWeight) << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.referenceHistogram, 0.999, h_csWeightCycle.maxReferenceWeight) << ","
					<< h_csWeightCycle.maxNpsWeight << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.npsHistogram, 0.99, h_csWeightCycle.maxNpsWeight) << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.npsHistogram, 0.999, h_csWeightCycle.maxNpsWeight) << ","
					<< h_csWeightCycle.maxNpsRatio << ","
					<< csFraction(h_csWeightCycle.npsAboveReference, h_csWeightCycle.observations) << ","
					<< h_csWeightCycle.maxNpsIncrease << "," << h_csWeightCycle.maxNpsDecrease << ","
					<< h_csWeightCycle.maxSourceWeight << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.sourceHistogram, 0.99, h_csWeightCycle.maxSourceWeight) << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.sourceHistogram, 0.999, h_csWeightCycle.maxSourceWeight) << ","
					<< h_csWeightCycle.maxSourceRatio << ","
					<< csFraction(h_csWeightCycle.sourceAboveReference, h_csWeightCycle.observations) << ","
					<< h_csWeightCycle.maxSourceIncrease << "," << h_csWeightCycle.maxSourceDecrease << ","
					<< h_csWeightCycle.maxFissionWeight << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.fissionHistogram, 0.99, h_csWeightCycle.maxFissionWeight) << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.fissionHistogram, 0.999, h_csWeightCycle.maxFissionWeight) << ","
					<< h_csWeightCycle.maxFissionRatio << ","
					<< h_csWeightCycle.maxScoreWeight << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.scoreHistogram, 0.99, h_csWeightCycle.maxScoreWeight) << ","
					<< csHistogramPercentileBounded(h_csWeightCycle.scoreHistogram, 0.999, h_csWeightCycle.maxScoreWeight) << ","
					<< h_csWeightCycle.maxScoreRatio << ","
					<< csFraction(h_csWeightCycle.scoreAboveReference, h_csWeightCycle.observations) << ","
					<< h_csWeightCycle.maxScoreIncrease << "," << h_csWeightCycle.maxScoreDecrease << "\n";
				csWeightCycleLog.flush();
			}
			normalizeCorrelatedSource << <blockPerDim, threadPerBlock >> > (
				d_Bank, csConfig, d_csSourceSums, d_csSourceCount);
#ifdef _DEBUG
			CUDA_KERNEL_CHECK();
#endif
		}

#ifdef TALLYLASTCYCLE
		if (cycle >= inactiveCycle) {
			GPU_Manager::FetchCoreTallyToHost(
				h_CoreTally,
				d_bufferTallyAssembly,
				d_bufferTallyPincellVec
			);
			pinPowerCycleStatistics.AccumulateActiveCycle(h_Core, h_CoreTally);
			for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
				const int cellCount = h_CoreTally.tallyAssembly[a].totalPincellNo();
				for (int cell = 0; cell < cellCount; ++cell) {
					cumulativeActivePinTally[a][cell] +=
						h_CoreTally.tallyAssembly[a].pinCells[cell].pinTally;
					cumulativeActiveModTally[a][cell] +=
						h_CoreTally.tallyAssembly[a].pinCells[cell].modTally;
				}
			}
		}
#endif


		cudaMemcpy(&h_fissionCount, d_fissionCount, sizeof(double), cudaMemcpyDeviceToHost);

		double rawRatio = h_fissionCount / fissionCountBuffer;
		/*
		if (cycle < 5) {
			// prevents early cycle cycling
			rawRatio = std::max(0.6, std::min(1.4, rawRatio));
		}
		*/
		
		double collisionEstK = rawRatio;


		//std::cout << "Previous fission: " << fissionCountBuffer << ", Current Fission Count: " << h_fissionCount << ",\tCollisionEstimator: " << h_fissionCount / fissionCountBuffer << "\n";
		std::cout << "initial Weight: " << initialWeight	<< ", Current Fission strength: " << h_fissionCount << ", previous Fission strength: " << fissionCountBuffer << ", CollisionEstimator: " << collisionEstK << "\n";
		// this just fetches only the number of neutrons
		cudaMemcpy(&h_Bank, d_Bank, sizeof(NeutronBank), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_multK, d_multK, sizeof(double), cudaMemcpyDeviceToHost);
		
		fissionCountBuffer = h_fissionCount;
		
		currentNumNeutron = h_Bank.getTotalNeutronNum();
		//currentNumNeutron = h_NeutronNum;
		double oldK = h_multK;
		std::cout << "Cycle " << cycle + 1 << ", currentNum: " << currentNumNeutron << " consisting of original + added bank: " << h_Bank.neutronSize << " + " << h_Bank.addedNeutronSize << "\n";
		//h_multK = h_multK * currentNumNeutron / previousNumNeutron;
		h_multK *= collisionEstK;
		previousNumNeutron = currentNumNeutron;
		std::cout << "\tk: " << h_multK << "\n";
		//std::cout << "\t n count : " << h_Bank.neutronSize << " addn count : " << h_Bank.addedNeutronSize << " addN addIndex : " << h_Bank.addedNeutronIndex;

		if (oldK == h_multK) {
			errorCounter++;
			if (errorCounter == 10) {
				std::cout << "Error in GPU Memory - early termination\n";
				break;
			}
		}
		else {
			errorCounter = 0;
		}

		cudaMemcpy(d_multK, &h_multK, sizeof(double), cudaMemcpyHostToDevice);

		//h_NeutronNum = 0.0;

		h_NeutronNum = currentNumNeutron;
		cudaMemcpy(d_NeutronNum, &h_NeutronNum, sizeof(double), cudaMemcpyHostToDevice);
		// i guess you can use memset?
		//cudaMemset(d_NeutronNum, 0.0, sizeof(double));



#ifdef TALLY
		//GPU_Manager::FetchCoreTallyToHost(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
		if (cycle > 180 && cycle < 199) {
			if (cycle % 2 == 0) {
				std::cout << "Fetching info from core structure flux tally, cycle: " << cycle + 1 << "\n";
				GPU_Manager::FetchCoreTallyToHost(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
				DumpCoreTallyZPositionToText(h_CoreTally, fluxTallyLog, cycle + 1, runConfig.tallyZPositionCm, h_multK);
			}
			ResetCoreTallyOnDevice(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
		}

		if (cycle == 1) {
			std::cout << "Fetching info from core structure flux tally, cycle: " << cycle + 1 << "\n";
			GPU_Manager::FetchCoreTallyToHost(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
			// fetch the tally
			DumpCoreTallyZPositionToText(h_CoreTally, fluxTallyLog, cycle + 1, runConfig.tallyZPositionCm, h_multK);
		}
		if (cycle % tallyFetchCycleSpec == 0 && cycle != 0) {
			std::cout << "Fetching info from core structure flux tally, cycle: " << cycle << "\n";
			GPU_Manager::FetchCoreTallyToHost(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
			// fetch the tally
			DumpCoreTallyZPositionToText(h_CoreTally, fluxTallyLog, cycle + 1, runConfig.tallyZPositionCm, h_multK);
		}
		if (cycle % 20 == 0 && cycle < 199) {
			std::cout << "Fetching info from core structure flux tally, cycle: " << cycle + 1 << "\n";
			GPU_Manager::FetchCoreTallyToHost(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
			// fetch the tally
			DumpCoreTallyZPositionToText(h_CoreTally, fluxTallyLog, cycle + 1, runConfig.tallyZPositionCm, h_multK);
		}
		if (cycle % tallyFetchCycleSpec - 1 == 0) {
			ResetCoreTallyOnDevice(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
		}
		if (cycle == 0) {
			ResetCoreTallyOnDevice(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
		}
		if ((cycle + 1) % 20 == 0 && cycle + 1 < 199) {
			ResetCoreTallyOnDevice(h_CoreTally, d_bufferTallyAssembly, d_bufferTallyPincellVec);
		}

#endif



		// Welford Algorithm for calculating the running mean
		klog << (cycle + 1) << " " << h_multK << "\n";
		if (csConfig.enabled) {
			csLog << (cycle + 1) << " " << h_multK;
			const double npsK = h_csFissionScores[0] / initialWeight;
			for (int level = 0; level <= csConfig.propagationGenerations; ++level) {
				const double perturbedK = h_csFissionScores[level] / initialWeight;
				const double deltaK = perturbedK - h_multK;
				csLog << " " << perturbedK << " " << deltaK;
				if (cycle >= inactiveCycle && cycle >= level) {
					csKStatistics[level].add(perturbedK);
					csDeltaKStatistics[level].add(deltaK);
					csDeltaKSamples[level].push_back(deltaK);
					csSourceEffectStatistics[level].add(perturbedK - npsK);
					csSourceEffectSamples[level].push_back(perturbedK - npsK);
					csPairedNpsStatistics[level].add(npsK);
					if (level > 0) {
						const double previousLevelK = h_csFissionScores[level - 1] / initialWeight;
						const double generationIncrement = perturbedK - previousLevelK;
						csGenerationIncrementStatistics[level].add(generationIncrement);
						csGenerationIncrementSamples[level].push_back(generationIncrement);
					}
					if (h_multK > 0.0 && perturbedK > 0.0) {
						const double deltaRho = 1.0 / h_multK - 1.0 / perturbedK;
						csDeltaRhoStatistics[level].add(deltaRho);
						csDeltaRhoSamples[level].push_back(deltaRho);
					}
				}
			}
			csLog << " " << h_csSourceCount;
			for (int level = 1; level <= csConfig.propagationGenerations; ++level) {
				const double sourceMean = h_csSourceCount > 0.0
					? h_csSourceSums[level] / h_csSourceCount : 0.0;
				csLog << " " << sourceMean;
			}
			csLog << "\n";

			const int selectedLevel = csConfig.propagationGenerations;
			std::cout << "\tCS: k_NPS(L0)=" << h_csFissionScores[0] / initialWeight
				<< ", k_total(L" << selectedLevel << ")=" << h_csFissionScores[selectedLevel] / initialWeight
				<< ", delta_k=" << h_csFissionScores[selectedLevel] / initialWeight - h_multK << "\n";
		}
		if (cycle >= inactiveCycle) {
			activeCount += 1;
			double x = h_multK;
			activeKSamples.push_back(x);
			double delta = x - meanK;
			meanK += delta / activeCount;
			double delta2 = x - meanK;
			M2 += delta * delta2;
		}

		
		std::cout << "Merging:\n";
		GPU_Manager::compact_bank_device(d_Bank);
		cudaMemcpy(&h_Bank, d_Bank, sizeof(NeutronBank), cudaMemcpyDeviceToHost);
		
		std::cout << "End Of Cycle\n" << std::flush;
	}
	auto t_end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed = t_end - t_start;
	if (csConfig.enabled) {
		cudaMemcpy(&h_csSupportMismatch, d_csSupportMismatch,
			sizeof(unsigned long long), cudaMemcpyDeviceToHost);
	}
	if (csConfig.enabled && runConfig.csTrackWeights) {
		std::ofstream histogramLog(csWeightHistogramOss.str(), std::ios::out);
		histogramLog << std::setprecision(12)
			<< "bin,log10_lower,log10_upper,w,wstar_nps,wstar_source,wstar_fission,wstar_score,"
			<< "ratio_nps,ratio_source,ratio_score\n";
		for (int bin = 0; bin < CS_WEIGHT_HIST_BINS; ++bin) {
			const double lower = bin == 0 ? -std::numeric_limits<double>::infinity()
				: CS_WEIGHT_HIST_LOG10_MIN + static_cast<double>(bin - 1) * CS_WEIGHT_HIST_LOG10_WIDTH;
			const double upper = bin == CS_WEIGHT_HIST_BINS - 1 ? std::numeric_limits<double>::infinity()
				: CS_WEIGHT_HIST_LOG10_MIN + static_cast<double>(bin) * CS_WEIGHT_HIST_LOG10_WIDTH;
			histogramLog << bin << "," << lower << "," << upper << ","
				<< h_csWeightActive.referenceHistogram[bin] << ","
				<< h_csWeightActive.npsHistogram[bin] << ","
				<< h_csWeightActive.sourceHistogram[bin] << ","
				<< h_csWeightActive.fissionHistogram[bin] << ","
				<< h_csWeightActive.scoreHistogram[bin] << ","
				<< h_csWeightActive.npsRatioHistogram[bin] << ","
				<< h_csWeightActive.sourceRatioHistogram[bin] << ","
				<< h_csWeightActive.scoreRatioHistogram[bin] << "\n";
		}

		std::ofstream summaryLog(csWeightSummaryOss.str(), std::ios::out);
		summaryLog << std::setprecision(12)
			<< "CS weight diagnostic (active cycles only)\n"
			<< "Reference rod depth (cm): " << csConfig.referenceRodDepth << "\n"
			<< "Perturbed rod depth (cm): " << csConfig.perturbedRodDepth << "\n"
			<< "Source propagation level: " << csConfig.propagationGenerations << "\n"
			<< "Fission-score observations (exact maxima/counters): " << h_csWeightActive.observations << "\n"
			<< "Histogram observations (1/" << CS_WEIGHT_HIST_SAMPLE_STRIDE << " particle indices): "
			<< h_csWeightActive.histogramObservations << "\n"
			<< "Invalid observations: " << h_csWeightActive.invalidObservations << "\n\n"
			<< "Definitions:\n"
			<< "  w                  = reference physical particle weight\n"
			<< "  wstar_nps          = w * exp(csLogLikelihood)\n"
			<< "  wstar_source       = wstar_nps * R_L\n"
			<< "  wstar_fission      = wstar_nps * (nuSigmaF*/nuSigmaF)\n"
			<< "  wstar_score        = wstar_fission * R_L\n\n";

		auto writeMetric = [&](const char* name, const unsigned long long* histogram,
			double maximum, double ratioMaximum, unsigned long long aboveReference,
			double increase, double decrease) {
			summaryLog << name << ":\n"
				<< "  p50/p90/p95/p99/p99.9/p99.99 = "
				<< csHistogramPercentileBounded(histogram, 0.50, maximum) << " / "
				<< csHistogramPercentileBounded(histogram, 0.90, maximum) << " / "
				<< csHistogramPercentileBounded(histogram, 0.95, maximum) << " / "
				<< csHistogramPercentileBounded(histogram, 0.99, maximum) << " / "
				<< csHistogramPercentileBounded(histogram, 0.999, maximum) << " / "
				<< csHistogramPercentileBounded(histogram, 0.9999, maximum) << "\n"
				<< "  exact max = " << maximum << "\n";
			if (ratioMaximum >= 0.0) {
				summaryLog << "  exact max(value/w) = " << ratioMaximum << "\n"
					<< "  fraction(value>w) = "
					<< csFraction(aboveReference, h_csWeightActive.observations) << "\n"
					<< "  max(value-w) / max(w-value) = " << increase << " / " << decrease << "\n";
			}
		};
		writeMetric("w", h_csWeightActive.referenceHistogram,
			h_csWeightActive.maxReferenceWeight, -1.0, 0ULL, 0.0, 0.0);
		writeMetric("wstar_nps", h_csWeightActive.npsHistogram,
			h_csWeightActive.maxNpsWeight, h_csWeightActive.maxNpsRatio,
			h_csWeightActive.npsAboveReference, h_csWeightActive.maxNpsIncrease,
			h_csWeightActive.maxNpsDecrease);
		writeMetric("wstar_source", h_csWeightActive.sourceHistogram,
			h_csWeightActive.maxSourceWeight, h_csWeightActive.maxSourceRatio,
			h_csWeightActive.sourceAboveReference, h_csWeightActive.maxSourceIncrease,
			h_csWeightActive.maxSourceDecrease);
		writeMetric("wstar_fission", h_csWeightActive.fissionHistogram,
			h_csWeightActive.maxFissionWeight, h_csWeightActive.maxFissionRatio,
			h_csWeightActive.fissionAboveReference, h_csWeightActive.maxFissionIncrease,
			h_csWeightActive.maxFissionDecrease);
		writeMetric("wstar_score", h_csWeightActive.scoreHistogram,
			h_csWeightActive.maxScoreWeight, h_csWeightActive.maxScoreRatio,
			h_csWeightActive.scoreAboveReference, h_csWeightActive.maxScoreIncrease,
			h_csWeightActive.maxScoreDecrease);

		summaryLog << "\nApproximate sampled exceedance fractions for source-aware particle weight wstar_source:\n";
		for (double threshold : {0.5, 1.0, 2.0, 5.0, 10.0}) {
			summaryLog << "  > " << threshold << ": "
				<< csHistogramFractionAbove(h_csWeightActive.sourceHistogram, threshold) << "\n";
		}
		summaryLog << "\nCandidate upper bounds from sampled wstar_source distribution:\n"
			<< "  p99: " << csHistogramPercentileBounded(h_csWeightActive.sourceHistogram, 0.99, h_csWeightActive.maxSourceWeight) << "\n"
			<< "  p99.9: " << csHistogramPercentileBounded(h_csWeightActive.sourceHistogram, 0.999, h_csWeightActive.maxSourceWeight) << "\n"
			<< "  p99.99: " << csHistogramPercentileBounded(h_csWeightActive.sourceHistogram, 0.9999, h_csWeightActive.maxSourceWeight) << "\n";

		std::cout << "CS weight diagnostics written to:\n"
			<< "  " << csWeightCycleOss.str() << "\n"
			<< "  " << csWeightHistogramOss.str() << "\n"
			<< "  " << csWeightSummaryOss.str() << "\n";
	}

	std::cout << "Total Time: " << elapsed.count() << "seconds. Average: " << elapsed.count() / cycleNum << " seconds per cycle.\n";

	if (activeCount >= 2) {
		double var = M2 / (activeCount - 1);
		double stddev = std::sqrt(var);
		double iid_stderr_mean = stddev / std::sqrt(static_cast<double>(activeCount));
		const StandardErrorEstimate kError = estimateStandardError(activeKSamples, iid_stderr_mean);
		double stderr_mean = kError.conservative;


		if (runConfig.runCase == RunCase::Unrodded) {
			std::cout << "\n\nRunning Unrodded case. Requested rod depth: " << runConfig.rodInsertionDepth
				<< " cm, actual modeled rod depth: " << actualRodDepth << " cm\n";
			klog << "\n\nRunning Unrodded case. Requested rod depth: " << runConfig.rodInsertionDepth
				<< " cm, actual modeled rod depth: " << actualRodDepth << " cm\n";
		}
		else if (runConfig.runCase == RunCase::RoddedA) {
			std::cout << "\n\nRunning Rodded A case\n";
			klog << "\n\nRunning Rodded A case\n";
		}
		else if (runConfig.runCase == RunCase::RoddedB) {
			std::cout << "\n\nRunning Rodded B case\n";
			klog << "\n\nRunning Rodded B case\n";
		}
		else if (runConfig.runCase == RunCase::CriticalSearch_Binary) {
			std::cout << "\n\nRunning Critical Search (Binary) case. Requested rod depth: " << runConfig.rodInsertionDepth
				<< " cm, actual modeled rod depth: " << actualRodDepth << " cm\n";
			klog << "\n\nRunning Critical Search (Binary) case. Requested rod depth: " << runConfig.rodInsertionDepth
				<< " cm, actual modeled rod depth: " << actualRodDepth << " cm\n";
		}
		else if (runConfig.runCase == RunCase::Reference2003) {
			std::cout << "\n\nRunning Reference 2003 case\n";
			klog << "\n\nRunning Reference 2003 case\n";
		}
		else if (runConfig.runCase == RunCase::ARO) {
			std::cout << "\n\nRunning ARO case - only top moderator gets CR\n";
			klog << "\n\nRunning ARO case - only top moderator gets CR\n";
		}

		std::cout << "Active cycles: " << activeCount << "\n";
		std::cout << "Requested rod depth: " << runConfig.rodInsertionDepth << " cm\n";
		std::cout << "Actual modeled rod depth: " << actualRodDepth << " cm\n";
		std::cout << "k_mean: " << std::fixed << std::setprecision(7) << meanK << "\n";
		std::cout << "k_stddev(cycle): " << stddev << "\n";
		std::cout << "Difference with the reference K: " << meanK - referenceK << ", Per cent error: " << (meanK - referenceK) / referenceK * 100 << "\n"; 
		std::cout << "k_stderr(mean): " << stderr_mean << "\n";
		std::cout << "k_stderr(iid-cycle): " << kError.iidCycle << "\n";
		std::cout << "k_stderr(batch-means): " << kError.batchMeans << "\n";
		std::cout << "k_batch_means: " << kError.batchCount << " batches x " << kError.batchSize << " cycles\n";
		std::cout << "Reference K: " << referenceK << "\n";
		std::cout << "Reactivity(pcm): " << (meanK - 1.0) * 100000 / meanK << ", pcm difference: " << ((meanK - 1.0) / meanK) * 100000 - ((referenceK - 1.0) / referenceK) * 100000 << "\n";
		std::cout << "stddev in pcm: " << stddev / (meanK * meanK) * 100000 << ", stderr of mean in pcm: " << stderr_mean / (meanK * meanK) * 100000 << "\n";
		std::cout << "Total Time: " << elapsed.count() << "seconds. Average: " << elapsed.count() / cycleNum << " seconds per cycle.\n";

		

		klog << "\n\nActive cycles: " << activeCount << "\n";
		klog << "Requested rod depth: " << runConfig.rodInsertionDepth << " cm\n";
		klog << "Actual modeled rod depth: " << actualRodDepth << " cm\n";
		klog << "k_mean: " << meanK << "\n";
		klog << "k_stddev(cycle): " << stddev << "\n";
		klog << "k_stderr(mean): " << stderr_mean << "\n";
		klog << "k_stderr(iid-cycle): " << kError.iidCycle << "\n";
		klog << "k_stderr(batch-means): " << kError.batchMeans << "\n";
		klog << "k_batch_means: " << kError.batchCount << " batches x " << kError.batchSize << " cycles\n";
		klog << "Reference K: " << referenceK << "\n";
		klog << "Difference with the reference K: " << meanK - referenceK  << ", Per cent error: " << (meanK - referenceK) / referenceK * 100 << "\n";
		klog << "Reactivity(pcm): " << (meanK - 1.0) * 100000 / meanK << ", pcm difference: " << ((meanK - 1.0) / meanK) * 100000 - ((referenceK - 1.0) / referenceK) * 100000 << "\n";
		klog << "stddev in pcm: " << stddev / (meanK * meanK) * 100000 << ", stderr of mean in pcm: " << stderr_mean / (meanK * meanK) * 100000 << "\n";
		klog << "Total Time: " << elapsed.count() << "seconds. Average: " << elapsed.count() / cycleNum << " seconds per cycle.\n";
	
	}

	if (csConfig.enabled) {
		const int selectedLevel = csConfig.propagationGenerations;
		const RunningStatistics& selectedK = csKStatistics[selectedLevel];
		const RunningStatistics& selectedDeltaK = csDeltaKStatistics[selectedLevel];
		const RunningStatistics& selectedDeltaRho = csDeltaRhoStatistics[selectedLevel];
		const RunningStatistics& npsDeltaK = csDeltaKStatistics[0];
		const RunningStatistics& sourceEffect = csSourceEffectStatistics[selectedLevel];
		const RunningStatistics& lastGenerationIncrement = csGenerationIncrementStatistics[selectedLevel];
		const StandardErrorEstimate npsDeltaKError = estimateStandardError(
			csDeltaKSamples[0], npsDeltaK.standardError());
		const StandardErrorEstimate sourceEffectError = estimateStandardError(
			csSourceEffectSamples[selectedLevel], sourceEffect.standardError());
		const StandardErrorEstimate lastGenerationError = estimateStandardError(
			csGenerationIncrementSamples[selectedLevel], lastGenerationIncrement.standardError());
		const StandardErrorEstimate selectedDeltaKError = estimateStandardError(
			csDeltaKSamples[selectedLevel], selectedDeltaK.standardError());
		const StandardErrorEstimate selectedDeltaRhoError = estimateStandardError(
			csDeltaRhoSamples[selectedLevel], selectedDeltaRho.standardError());

		std::ostringstream csSummary;
		csSummary << std::setprecision(10)
			<< "\nCorrelated-sampling result\n"
			<< "Reference rod depth (cm): " << csConfig.referenceRodDepth << "\n"
			<< "Perturbed rod depth (cm): " << csConfig.perturbedRodDepth << "\n"
			<< "Selected source-propagation level: " << selectedLevel << "\n"
			<< "Usable active cycles: " << selectedK.count << "\n";
		if (selectedK.count > 0) {
			const double pairedReferenceMean = selectedK.mean - selectedDeltaK.mean;
			csSummary << "Paired reference k: " << pairedReferenceMean << "\n"
				<< "Perturbed k, NPS only (paired cycles): " << csPairedNpsStatistics[selectedLevel].mean << "\n"
				<< "Perturbed k, NPS + fission-source propagation (L" << selectedLevel << "): " << selectedK.mean << "\n"
				<< "Delta k, NPS only: " << npsDeltaK.mean
				<< " +/- " << npsDeltaKError.conservative << "\n"
				<< "Additional source effect, L" << selectedLevel << " - L0: " << sourceEffect.mean
				<< " +/- " << sourceEffectError.conservative << "\n"
				<< "Last propagation increment, L" << selectedLevel << " - L" << selectedLevel - 1
				<< ": " << lastGenerationIncrement.mean
				<< " +/- " << lastGenerationError.conservative << "\n"
				<< "Total delta k: " << selectedDeltaK.mean
				<< " +/- " << selectedDeltaKError.conservative << "\n"
				<< "Total delta reactivity (pcm): " << selectedDeltaRho.mean * 1.0e5
				<< " +/- " << selectedDeltaRhoError.conservative * 1.0e5 << "\n"
				<< "CS stderr estimator: max(iid-cycle, non-overlapping batch-means)\n"
				<< "CS batch means: " << selectedDeltaKError.batchCount << " batches x "
				<< selectedDeltaKError.batchSize << " cycles\n"
				<< "Total delta k iid-cycle stderr: " << selectedDeltaKError.iidCycle << "\n"
				<< "Total delta k batch-means stderr: " << selectedDeltaKError.batchMeans << "\n"
				<< "Last propagation increment iid-cycle stderr: " << lastGenerationError.iidCycle << "\n"
				<< "Last propagation increment batch-means stderr: " << lastGenerationError.batchMeans << "\n";
			if (lastGenerationIncrement.count >= 2 &&
				fabs(lastGenerationIncrement.mean) > 2.0 * lastGenerationError.conservative) {
				csSummary << "NOTICE: the last source-propagation increment is statistically resolved; "
					<< "the search driver will assess its practical size and probe-slope consistency.\n";
			}
		}
		csSummary << "Support/non-finite diagnostic count: " << h_csSupportMismatch << "\n";
		if (h_csSupportMismatch != 0ULL) {
			csSummary << "WARNING: nonzero diagnostic count; do not use this CS estimate until the support/weight issue is resolved.\n";
		}

		std::cout << csSummary.str();
		csLog << csSummary.str();
		std::cout << "CS cycle history: " << csOss.str() << "\n";
	}

#ifdef TALLYLASTCYCLE
	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const int cellCount = h_CoreTally.tallyAssembly[a].totalPincellNo();
		for (int cell = 0; cell < cellCount; ++cell) {
			h_CoreTally.tallyAssembly[a].pinCells[cell].pinTally =
				cumulativeActivePinTally[a][cell];
			h_CoreTally.tallyAssembly[a].pinCells[cell].modTally =
				cumulativeActiveModTally[a][cell];
		}
	}
	DumpCoreTallyZPositionToText(
		h_CoreTally,
		fluxTallyLog,
		cycleNum,
		runConfig.tallyZPositionCm,
		meanK
	);
	DumpBenchmarkPinPowerSummaryToText(h_Core, h_CoreTally, fluxTallyLog);
	pinPowerCycleStatistics.DumpToText(pinPowerCycleStatisticsLog);
	std::cout << "Pin-power active-cycle statistics: " << oss4.str() << "\n";

#ifdef C5G7_DUMP_FULL_CORE_PINWISE_TALLY
	DumpFullCorePinwiseTallyToText(
		h_CoreTally,
		pinwiseTallyLog,
		cycleNum,
		activeCount,
		meanK
	);
	std::cout << "Full-core pin-wise tally: " << oss3.str() << "\n";
#endif
	
#endif


	//std::cout << "\n initial Neutron number: " << num << ", for cycle of " << cycleNum - inactiveCycle << ", average k = " << meanK << "\n";

	klog.close();
	if (csLog.is_open()) csLog.close();
	if (csWeightCycleLog.is_open()) csWeightCycleLog.close();
	fluxTallyLog.close();
#ifdef C5G7_DUMP_FULL_CORE_PINWISE_TALLY
	pinwiseTallyLog.close();
#endif
#ifdef TALLYLASTCYCLE
	pinPowerCycleStatisticsLog.close();
#endif

	delete[] h_SeedArr;
	h_SeedArr = nullptr;

	if (d_bufferNeutrons) cudaFree(d_bufferNeutrons);
	if (d_bufferAddedNeutrons) cudaFree(d_bufferAddedNeutrons);
	if (d_Bank) cudaFree(d_Bank);

	if (d_Neutron) cudaFree(d_Neutron);

	for (auto& p : d_bufferPincellVec) {
		if (p) cudaFree(p);
		p = nullptr;
	}

	for (auto& p : d_bufferTallyPincellVec) {
		if (p) cudaFree(p);
		p = nullptr;
	}

	for (auto& p : d_bufferMortonRankVec) {
		if (p) cudaFree(p);
		p = nullptr;
	}

	for (auto& p : d_bufferTallyMortonRankVec) {
		if (p) cudaFree(p);
		p = nullptr;
	}

	if (d_bufferAssembly) cudaFree(d_bufferAssembly);
	if (d_bufferTallyAssembly) cudaFree(d_bufferTallyAssembly);

	if (d_Core) cudaFree(d_Core);
	if (d_CoreTally) cudaFree(d_CoreTally);
	if (d_XSLib) cudaFree(d_XSLib);

	if (d_SeedArr) cudaFree(d_SeedArr);

	if (d_multK) cudaFree(d_multK);
	if (d_NeutronNum) cudaFree(d_NeutronNum);
	if (d_NeutronWeightSum) cudaFree(d_NeutronWeightSum);
	if (d_NeutronWeightModifier) cudaFree(d_NeutronWeightModifier);
	if (d_fissionCount) cudaFree(d_fissionCount);
	if (d_csFissionScores) cudaFree(d_csFissionScores);
	if (d_csSourceSums) cudaFree(d_csSourceSums);
	if (d_csSourceCount) cudaFree(d_csSourceCount);
	if (d_csSupportMismatch) cudaFree(d_csSupportMismatch);
	if (d_csWeightDiagnostics) cudaFree(d_csWeightDiagnostics);

	// Optional: reset device (helps detect leaks in some tools)
	// this would really help
	cudaDeviceReset();
}
