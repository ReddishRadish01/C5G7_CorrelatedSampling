#pragma once

#include "cudaHeader.cuh"
#include <cuda_runtime.h>
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "GpuManager.cuh"
#include "CoreManager.cuh"
#include "XSManager.cuh"
#include "NeutronBankManager.cuh"
#include "Debug.cuh"
#include "TallyConfig.cuh"
#include <iomanip>
#include <array>
#include <limits>
#include <map>
#include <sstream>


class TallyPincell {
public:
	double sideLength = 1.26;	// 1.26 cm. alias with height
	double height = 1.26;
	double radius = 0.54; 
	double pinTally = 0.0;
	double modTally = 0.0;

	H TallyPincell() = default;

	H TallyPincell(double sideLength, double radius = 0.0, double height = 0.0) 
		: sideLength(sideLength), radius(radius), height(height)
	{
		if (2 * radius > sideLength) { this->radius = sideLength / 2.0; }
	}
};

class TallyAssembly {
public:
	vec3 startPos = { 0, 0, 0 };

	vec3 length = { 0, 0, 0 };

	int xNum = 0;
	int yNum = 0;
	int zNum = 0;
	TallyPincell* pinCells = nullptr;
	int* xyToMortonRank = nullptr;
	TallyPincell OOBPincell;

	H TallyAssembly() = default;

	H void Initialize(std::string assemblyTxt, double cellHeight = 0.0, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }) {

		std::ifstream assembly(assemblyTxt);

		if (!assembly.is_open()) {
			throw std::runtime_error("Cannot open: " + assemblyTxt);
		}

		if (!(assembly >> this->xNum >> this->yNum >> this->zNum)) {
			throw std::runtime_error("Failed to read xNum yNum zNum from: " + assemblyTxt);
		}

		//assembly >> this->xNum >> this->yNum >> this->zNum;
		int numPincells = this->xNum * this->yNum * this->zNum;

		double pincellLength = 0.0;
		double pincellHeight = 0.0;
		double radius = 0.0;

		assembly >> pincellLength >> pincellHeight >> radius;

		if (cellHeight != 0.0) {
			pincellHeight = cellHeight;
		}
		int mod = 0;
		assembly >> mod;

		MatType modType = (mod >= 0 && mod <= 7) ? static_cast<MatType>(mod) : MatType::MOD;

		this->pinCells = new TallyPincell[numPincells];
		this->xyToMortonRank = new int[this->xNum * this->yNum];
		buildMortonXYRank(this->xNum, this->yNum, this->xyToMortonRank);

		for (int i = 0; i < xNum * yNum; i++) {
			int x = i % xNum;
			int y = i / xNum;
			for (int j = 0; j < zNum; j++) {
				int index = pincellStorageIndex(x, y, j, xNum, yNum, this->xyToMortonRank);
				this->pinCells[index] = TallyPincell(pincellLength, radius, pincellHeight);
			}
			//std::cout << 
		}

		assembly.close();

		OOBPincell = TallyPincell(pincellLength, radius, pincellHeight);
	}

	H void Initialize_MOD(double cellSize, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }, double cellHeight = 0.0) {
		vec3 length = endPos - startPos;

		double height = 0.0;
		if (cellHeight == 0.0) {
			height = cellSize;
		}
		else {
			height = cellHeight;
		}

		this->startPos = startPos;
		this->length = length;

		if (cellSize <= 0.0 || height <= 0.0 || length.x <= 0.0 || length.y <= 0.0 || length.z <= 0.0) {
			this->xNum = 0;
			this->yNum = 0;
			this->zNum = 0;
			this->pinCells = nullptr;
			this->xyToMortonRank = nullptr;
			OOBPincell = TallyPincell(0.0, 0.0, 0.0);
			return;
		}

		this->xNum = static_cast<int>((length.x + 1.0e-9) / cellSize);
		this->yNum = static_cast<int>((length.y + 1.0e-9) / cellSize);
		this->zNum = static_cast<int>((length.z + 1.0e-9) / height);
		// fuck this floating point shits - the length.z is not enough - 

		int numPincells = this->xNum * this->yNum * this->zNum;
		this->pinCells = new TallyPincell[numPincells];
		this->xyToMortonRank = new int[this->xNum * this->yNum];
		buildMortonXYRank(this->xNum, this->yNum, this->xyToMortonRank);

		for (int i = 0; i < this->xNum * this->yNum; i++) {
			int x = i % this->xNum;
			int y = i / this->xNum;
			for (int j = 0; j < this->zNum; j++) {
				int index = pincellStorageIndex(x, y, j, this->xNum, this->yNum, this->xyToMortonRank);
				this->pinCells[index] = TallyPincell(cellSize, 0, height);
				//std::cout << static_cast<int>(this->pinCells[j].meatType) << " ";
			}
		}

	}

	HD TallyPincell& returnPincellByPos(Neutron& n);
	HD vec3 returnFlooredNeutronPosInPincell(Neutron& n);

	HD int totalPincellNo() {
		return xNum * yNum * zNum;
	}
};

class TallyC5G7Geometry {
public:
	double x; // cm
	double y;
	double z;
	double xMin = 0.0;
	double yMin = 0.0;
	double zMin = 0.0;

	TallyAssembly* tallyAssembly = nullptr;
	TallyAssembly nullAssembly{};
	int assemblyNo = 0;

	TallyC5G7Geometry(C5G7Geometry Core) {
		this->x = Core.x; this->y = Core.y; this->z = Core.z;
		this->xMin = Core.xMin; this->yMin = Core.yMin; this->zMin = Core.zMin;


		this->assemblyNo = Core.assemblyNo;
		tallyAssembly = new TallyAssembly[this->assemblyNo];
		for (int i = 0; i < this->assemblyNo; i++) {
			this->tallyAssembly[i].startPos = Core.assembly[i].startPos;
			this->tallyAssembly[i].length = Core.assembly[i].length;
			this->tallyAssembly[i].xNum = Core.assembly[i].xNum;
			this->tallyAssembly[i].yNum = Core.assembly[i].yNum;
			this->tallyAssembly[i].zNum = Core.assembly[i].zNum;
			int pinCellNo = this->tallyAssembly[i].xNum * this->tallyAssembly[i].yNum * this->tallyAssembly[i].zNum;
			this->tallyAssembly[i].pinCells = new TallyPincell[pinCellNo];
			this->tallyAssembly[i].xyToMortonRank = new int[this->tallyAssembly[i].xNum * this->tallyAssembly[i].yNum];
			buildMortonXYRank(this->tallyAssembly[i].xNum, this->tallyAssembly[i].yNum, this->tallyAssembly[i].xyToMortonRank);
			double sideLength = this->tallyAssembly[i].length.x / this->tallyAssembly[i].xNum;
			double height = this->tallyAssembly[i].length.z / this->tallyAssembly[i].zNum;
			for (int z = 0; z < this->tallyAssembly[i].zNum; z++) {
				for (int y = 0; y < this->tallyAssembly[i].yNum; y++) {
					for (int x = 0; x < this->tallyAssembly[i].xNum; x++) {
						int idx = pincellStorageIndex(x, y, z, this->tallyAssembly[i].xNum, this->tallyAssembly[i].yNum, this->tallyAssembly[i].xyToMortonRank);
						this->tallyAssembly[i].pinCells[idx].sideLength = sideLength;
						this->tallyAssembly[i].pinCells[idx].height = height;
						this->tallyAssembly[i].pinCells[idx].radius = 0.0;
					}
				}
			}
			this->tallyAssembly[i].OOBPincell = TallyPincell(1.26, 0.0, 1.26);
		}
	}
	

	HD TallyAssembly& returnAssemblyByNeutron(Neutron& n) {
		double eps = 1.0e-12;
		double x = n.dirVec.x > 0 ? n.dirVec.x * eps : n.dirVec.x * -eps;
		double y = n.dirVec.y > 0 ? n.dirVec.y * eps : n.dirVec.y * -eps;
		double z = n.dirVec.z > 0 ? n.dirVec.z * eps : n.dirVec.z * -eps;
		vec3 epsVec = { x, y, z };
		//Neutron localN = n;
		//localN.pos = localN.pos + epsVec;
		//n.pos = n.pos + epsVec;
		vec3 localPos = n.pos + epsVec;

		for (int i = 0; i < this->assemblyNo; i++) {
			vec3 endPos = this->tallyAssembly[i].startPos + this->tallyAssembly[i].length;
			if (localPos.x >= this->tallyAssembly[i].startPos.x && localPos.x < endPos.x) {
				if (localPos.y >= this->tallyAssembly[i].startPos.y && localPos.y < endPos.y) {
					if (localPos.z >= this->tallyAssembly[i].startPos.z && localPos.z < endPos.z) {
						return this->tallyAssembly[i];
					}
				}
			}
		}

		// this usually means fucked up - 
		// we are never meant to pass the out-of-bounds neutrons, or nullified neutrons in this function.
		// we never want the code to flow into this far, in this function - idK why it ended up here, maybe put some debugger outputs just in case
		//printf("Neutron Out-Of-Bounds in Position: (%f, %f, %f)\n", n.pos.x, n.pos.y, n.pos.z);
		//n.Nullify();
		// return null assembly 
		return this->nullAssembly;
	}
};

class Tally {
public:


	D inline static void kTally() {
		
	}



};


static inline void DumpCoreTallyZPositionToText(
	const TallyC5G7Geometry& h_CoreTally,
	std::ofstream& out,
	int cycleNo,
	double targetZPositionCm,
	double multiplicationFactor
) {
	out << std::setprecision(17);
	out << "\n\n\nCoreBounds(cm) "
		<< h_CoreTally.xMin << " " << h_CoreTally.x << " "
		<< h_CoreTally.yMin << " " << h_CoreTally.y << " "
		<< h_CoreTally.zMin << " " << h_CoreTally.z << "\n";
	out << "ZSlicePosition(cm) " << targetZPositionCm << "\n";
	out << "AssemblyNo " << h_CoreTally.assemblyNo << "\n";

	out << cycleNo << " th cycle:\n";
	out << multiplicationFactor << " multiplication Factor\n";

	int selectedAssemblyNo = 0;
	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		const double zStart = Asm.startPos.z;
		const double zEnd = zStart + Asm.length.z;
		if (targetZPositionCm >= zStart && targetZPositionCm < zEnd) {
			++selectedAssemblyNo;
		}
	}
	out << "SelectedAssemblyNo " << selectedAssemblyNo << "\n";

	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		const int nx = Asm.xNum;
		const int ny = Asm.yNum;
		const int nz = Asm.zNum;
		if (nx <= 0 || ny <= 0 || nz <= 0 || Asm.length.z <= 0.0) {
			continue;
		}

		const double zStart = Asm.startPos.z;
		const double zEnd = zStart + Asm.length.z;
		if (targetZPositionCm < zStart || targetZPositionCm >= zEnd) {
			continue;
		}

		const double dz = Asm.length.z / static_cast<double>(nz);
		const double localZCellPosition = (targetZPositionCm - zStart) / dz;
		int zCell = static_cast<int>(std::floor(localZCellPosition + 1.0e-10));
		zCell = std::max(0, std::min(nz - 1, zCell));
		const double cellZMin = zStart + zCell * dz;
		const double cellZMax = cellZMin + dz;

		out << "\n";
		out << "Assembly " << a << "\n";
		out << "StartPos " << Asm.startPos.x << " " << Asm.startPos.y << " " << Asm.startPos.z << "\n";
		out << "Length   " << Asm.length.x << " " << Asm.length.y << " " << Asm.length.z << "\n";
		out << "Dims     " << nx << " " << ny << " " << nz << "\n";
		out << "SelectedZCell " << zCell << "\n";
		out << "ZCellBounds(cm) " << cellZMin << " " << cellZMax << "\n";
		//out << "OOB      " << Asm.OOBPincell.pinTally << " " << Asm.OOBPincell.modTally << "\n";
		out << "Format   a k j i pinTally modTally\n";

		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const int idx = pincellStorageIndex(i, j, zCell, nx, ny, Asm.xyToMortonRank);
				const auto& pc = Asm.pinCells[idx];
				out << a << " " << zCell << " " << j << " " << i << " " << pc.pinTally << " " << pc.modTally << "\n";
			}
		}
	}
}

static inline void DumpFullCorePinwiseTallyToText(
	const TallyC5G7Geometry& h_CoreTally,
	std::ofstream& out,
	int completedCycleCount,
	int talliedCycleCount,
	double multiplicationFactor
) {
	out << std::defaultfloat << std::setprecision(17);
#ifdef C5G7_TALLY_PINWISE_FISSION_RATE
	out << "TallyType PinwiseFissionRateCollisionEstimator\n";
#elif defined(C5G7_TALLY_FISSION_SOURCE_DIST)
	out << "TallyType FissionSourceDistributionCollisionEstimator\n";
#endif
	out << "Scope full_core_all_assemblies_all_z_cells\n";
	out << "CompletedCycleCount " << completedCycleCount << "\n";
	out << "TalliedCycleCount " << talliedCycleCount << "\n";
	out << "MultiplicationFactor " << multiplicationFactor << "\n";
	out << "CoreBounds(cm) "
		<< h_CoreTally.xMin << " " << h_CoreTally.x << " "
		<< h_CoreTally.yMin << " " << h_CoreTally.y << " "
		<< h_CoreTally.zMin << " " << h_CoreTally.z << "\n";
	out << "AssemblyNo " << h_CoreTally.assemblyNo << "\n";

	long long totalCellCount = 0;
	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		totalCellCount += static_cast<long long>(Asm.xNum)
			* static_cast<long long>(Asm.yNum)
			* static_cast<long long>(Asm.zNum);
	}
	out << "TotalTallyCellCount " << totalCellCount << "\n";

	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		const int nx = Asm.xNum;
		const int ny = Asm.yNum;
		const int nz = Asm.zNum;
		if (nx <= 0 || ny <= 0 || nz <= 0 || Asm.pinCells == nullptr) {
			continue;
		}

		out << "\nAssembly " << a << "\n";
		out << "StartPos(cm) " << Asm.startPos.x << " " << Asm.startPos.y << " " << Asm.startPos.z << "\n";
		out << "Length(cm) " << Asm.length.x << " " << Asm.length.y << " " << Asm.length.z << "\n";
		out << "Dims " << nx << " " << ny << " " << nz << "\n";
		out << "Format a z_index j i pinTally modTally\n";

		for (int zCell = 0; zCell < nz; ++zCell) {
			for (int j = 0; j < ny; ++j) {
				for (int i = 0; i < nx; ++i) {
					const int idx = pincellStorageIndex(i, j, zCell, nx, ny, Asm.xyToMortonRank);
					const auto& pc = Asm.pinCells[idx];
					out << a << " " << zCell << " " << j << " " << i
						<< " " << pc.pinTally << " " << pc.modTally << "\n";
				}
			}
		}
	}
}

struct BenchmarkPinPowerAccumulator {
	double xCenter = 0.0;
	double yCenter = 0.0;
	std::array<double, 3> sliceRaw{ 0.0, 0.0, 0.0 };
	double coreRaw = 0.0;
};

struct BenchmarkPinPowerExtrema {
	double minRaw = std::numeric_limits<double>::infinity();
	double maxRaw = -std::numeric_limits<double>::infinity();
	double minX = 0.0;
	double minY = 0.0;
	double maxX = 0.0;
	double maxY = 0.0;
};

static inline void UpdateBenchmarkPinPowerExtrema(
	BenchmarkPinPowerExtrema& extrema,
	double raw,
	double xCenter,
	double yCenter
) {
	if (raw < extrema.minRaw) {
		extrema.minRaw = raw;
		extrema.minX = xCenter;
		extrema.minY = yCenter;
	}
	if (raw > extrema.maxRaw) {
		extrema.maxRaw = raw;
		extrema.maxX = xCenter;
		extrema.maxY = yCenter;
	}
}

static inline bool IsBenchmarkFuelMaterial(MatType material) {
	return material == MatType::UO2
		|| material == MatType::MOX4_3
		|| material == MatType::MOX7_0
		|| material == MatType::MOX8_7;
}

static inline MatType BenchmarkMaterialAtPinCenter(const Pincell& pinCell) {
	if (pinCell.hasAxialMeatSplit && 0.5 >= pinCell.axialMeatSplitFraction) {
		return pinCell.axialUpperMeatType;
	}
	return pinCell.meatType;
}

using BenchmarkPinKey = std::pair<long long, long long>;

struct BenchmarkPinCycleMoment {
	double xCenter = 0.0;
	double yCenter = 0.0;
	long long sampleCount = 0;
	double mean = 0.0;
	double M2 = 0.0;
};

class BenchmarkPinPowerCycleStatistics {
public:
	long long activeCycleCount = 0;
	std::map<BenchmarkPinKey, BenchmarkPinCycleMoment> pins;

	void AccumulateActiveCycle(
		const C5G7Geometry& h_Core,
		const TallyC5G7Geometry& h_CoreTally
	) {
		std::map<BenchmarkPinKey, BenchmarkPinPowerAccumulator> cyclePins;

		for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
			const auto& Asm = h_CoreTally.tallyAssembly[a];
			const auto& GeometryAsm = h_Core.assembly[a];
			if (Asm.xNum <= 0 || Asm.yNum <= 0 || Asm.zNum <= 0
				|| GeometryAsm.xNum != Asm.xNum
				|| GeometryAsm.yNum != Asm.yNum
				|| GeometryAsm.zNum != Asm.zNum) {
				continue;
			}

			const double dx = Asm.length.x / static_cast<double>(Asm.xNum);
			const double dy = Asm.length.y / static_cast<double>(Asm.yNum);
			for (int zCell = 0; zCell < Asm.zNum; ++zCell) {
				for (int j = 0; j < Asm.yNum; ++j) {
					for (int i = 0; i < Asm.xNum; ++i) {
						const int geometryIdx = pincellStorageIndex(
							i, j, zCell,
							GeometryAsm.xNum, GeometryAsm.yNum,
							GeometryAsm.xyToMortonRank
						);
						if (!IsBenchmarkFuelMaterial(
							BenchmarkMaterialAtPinCenter(GeometryAsm.pinCells[geometryIdx])
						)) {
							continue;
						}

						const int tallyIdx = pincellStorageIndex(
							i, j, zCell,
							Asm.xNum, Asm.yNum,
							Asm.xyToMortonRank
						);
						const double raw = Asm.pinCells[tallyIdx].pinTally;
						if (raw < 0.0 || !std::isfinite(raw)) {
							continue;
						}

						const double xCenter = Asm.startPos.x + (i + 0.5) * dx;
						const double yCenter = Asm.startPos.y + (j + 0.5) * dy;
						const BenchmarkPinKey key{
							static_cast<long long>(std::llround(xCenter * 1.0e9)),
							static_cast<long long>(std::llround(yCenter * 1.0e9))
						};
						auto& pin = cyclePins[key];
						pin.xCenter = xCenter;
						pin.yCenter = yCenter;
						pin.coreRaw += raw;
					}
				}
			}
		}

		double cycleTotalRaw = 0.0;
		for (const auto& entry : cyclePins) {
			cycleTotalRaw += entry.second.coreRaw;
		}
		if (cyclePins.empty() || !(cycleTotalRaw > 0.0)) {
			return;
		}

		const double normalization =
			static_cast<double>(cyclePins.size()) / cycleTotalRaw;
		++activeCycleCount;
		for (const auto& entry : cyclePins) {
			const double normalizedPinPower = entry.second.coreRaw * normalization;
			auto& moment = pins[entry.first];
			moment.xCenter = entry.second.xCenter;
			moment.yCenter = entry.second.yCenter;
			++moment.sampleCount;
			const double delta = normalizedPinPower - moment.mean;
			moment.mean += delta / static_cast<double>(moment.sampleCount);
			const double delta2 = normalizedPinPower - moment.mean;
			moment.M2 += delta * delta2;
		}
	}

	void DumpToText(std::ofstream& out) const {
		out << std::setprecision(17);
		out << "PinPowerActiveCycleStatistics\n";
		if (pins.empty() || activeCycleCount <= 0) {
			out << "Status no_active_cycle_pin_samples\n";
			return;
		}

		out << "Status ok\n";
		out << "TallyQuantity full_length_pinwise_fission_rate\n";
		out << "SampleUnit active_cycle\n";
		out << "ActiveCycleCount " << activeCycleCount << "\n";
		out << "FuelPinCount " << pins.size() << "\n";
		out << "FuelSelection geometry_material_mask\n";
		out << "Normalization each_active_cycle_total_power_equals_fuel_pin_count\n";
		out << "StandardErrorDefinition cycle_sample_stddev/sqrt(active_cycle_count)\n";
		out << "Format x_cm y_cm mean_normalized_pin_power cycle_stddev standard_error relative_standard_error_percent\n";

		const BenchmarkPinCycleMoment* minimum = nullptr;
		const BenchmarkPinCycleMoment* maximum = nullptr;
		for (const auto& entry : pins) {
			const auto& pin = entry.second;
			const double variance = pin.sampleCount > 1
				? pin.M2 / static_cast<double>(pin.sampleCount - 1)
				: 0.0;
			const double stddev = std::sqrt(std::max(0.0, variance));
			const double standardError = pin.sampleCount > 0
				? stddev / std::sqrt(static_cast<double>(pin.sampleCount))
				: 0.0;
			const double relativePercent = pin.mean != 0.0
				? 100.0 * standardError / std::abs(pin.mean)
				: 0.0;
			out << pin.xCenter << " " << pin.yCenter << " "
				<< pin.mean << " " << stddev << " "
				<< standardError << " " << relativePercent << "\n";

			if (minimum == nullptr || pin.mean < minimum->mean) {
				minimum = &pin;
			}
			if (maximum == nullptr || pin.mean > maximum->mean) {
				maximum = &pin;
			}
		}

		auto dumpExtremum = [&out](const char* label, const BenchmarkPinCycleMoment& pin) {
			const double variance = pin.sampleCount > 1
				? pin.M2 / static_cast<double>(pin.sampleCount - 1)
				: 0.0;
			const double stddev = std::sqrt(std::max(0.0, variance));
			const double standardError = pin.sampleCount > 0
				? stddev / std::sqrt(static_cast<double>(pin.sampleCount))
				: 0.0;
			const double relativePercent = pin.mean != 0.0
				? 100.0 * standardError / std::abs(pin.mean)
				: 0.0;
			out << label
				<< " MeanNormalizedPinPower " << pin.mean
				<< " CycleStdDev " << stddev
				<< " StandardError " << standardError
				<< " RelativeStandardErrorPercent " << relativePercent
				<< " XY(cm) " << pin.xCenter << " " << pin.yCenter
				<< "\n";
		};
		dumpExtremum("CoreMinimumPinPowerStatistics", *minimum);
		dumpExtremum("CoreMaximumPinPowerStatistics", *maximum);
	}
};

static inline void DumpBenchmarkPinPowerSummaryToText(
	const C5G7Geometry& h_Core,
	const TallyC5G7Geometry& h_CoreTally,
	std::ofstream& out
) {
	if (h_Core.assemblyNo != h_CoreTally.assemblyNo) {
		out << "\nPinPowerBenchmarkSummary\n";
		out << "Status geometry_tally_assembly_count_mismatch\n";
		return;
	}

	double fuelZMin = std::numeric_limits<double>::infinity();
	double fuelZMax = -std::numeric_limits<double>::infinity();

	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		const auto& GeometryAsm = h_Core.assembly[a];
		if (Asm.xNum <= 0 || Asm.yNum <= 0 || Asm.zNum <= 0 || Asm.length.z <= 0.0
			|| GeometryAsm.xNum != Asm.xNum
			|| GeometryAsm.yNum != Asm.yNum
			|| GeometryAsm.zNum != Asm.zNum) {
			continue;
		}
		const double dz = Asm.length.z / static_cast<double>(Asm.zNum);
		for (int zCell = 0; zCell < Asm.zNum; ++zCell) {
			bool layerHasFuel = false;
			for (int j = 0; j < Asm.yNum && !layerHasFuel; ++j) {
				for (int i = 0; i < Asm.xNum; ++i) {
					const int geometryIdx = pincellStorageIndex(
						i,
						j,
						zCell,
						GeometryAsm.xNum,
						GeometryAsm.yNum,
						GeometryAsm.xyToMortonRank
					);
					if (IsBenchmarkFuelMaterial(BenchmarkMaterialAtPinCenter(GeometryAsm.pinCells[geometryIdx]))) {
						layerHasFuel = true;
						break;
					}
				}
			}
			if (layerHasFuel) {
				const double zLow = Asm.startPos.z + zCell * dz;
				fuelZMin = std::min(fuelZMin, zLow);
				fuelZMax = std::max(fuelZMax, zLow + dz);
			}
		}
	}

	if (!std::isfinite(fuelZMin) || !std::isfinite(fuelZMax) || fuelZMax <= fuelZMin) {
		out << "\nPinPowerBenchmarkSummary\n";
		out << "Status no_geometry_fuel_cells\n";
		return;
	}

	const double sliceHeight = (fuelZMax - fuelZMin) / 3.0;
	std::map<BenchmarkPinKey, BenchmarkPinPowerAccumulator> pins;

	for (int a = 0; a < h_CoreTally.assemblyNo; ++a) {
		const auto& Asm = h_CoreTally.tallyAssembly[a];
		const auto& GeometryAsm = h_Core.assembly[a];
		if (Asm.xNum <= 0 || Asm.yNum <= 0 || Asm.zNum <= 0
			|| GeometryAsm.xNum != Asm.xNum
			|| GeometryAsm.yNum != Asm.yNum
			|| GeometryAsm.zNum != Asm.zNum) {
			continue;
		}
		const double dx = Asm.length.x / static_cast<double>(Asm.xNum);
		const double dy = Asm.length.y / static_cast<double>(Asm.yNum);
		const double dz = Asm.length.z / static_cast<double>(Asm.zNum);

		for (int zCell = 0; zCell < Asm.zNum; ++zCell) {
			const double zCenter = Asm.startPos.z + (zCell + 0.5) * dz;
			if (zCenter < fuelZMin || zCenter >= fuelZMax) {
				continue;
			}
			int sliceIndex = static_cast<int>((zCenter - fuelZMin) / sliceHeight);
			sliceIndex = std::max(0, std::min(2, sliceIndex));

			for (int j = 0; j < Asm.yNum; ++j) {
				for (int i = 0; i < Asm.xNum; ++i) {
					const int geometryIdx = pincellStorageIndex(
						i,
						j,
						zCell,
						GeometryAsm.xNum,
						GeometryAsm.yNum,
						GeometryAsm.xyToMortonRank
					);
					if (!IsBenchmarkFuelMaterial(BenchmarkMaterialAtPinCenter(GeometryAsm.pinCells[geometryIdx]))) {
						continue;
					}

					const int idx = pincellStorageIndex(i, j, zCell, Asm.xNum, Asm.yNum, Asm.xyToMortonRank);
					const double raw = Asm.pinCells[idx].pinTally;
					if (raw < 0.0 || !std::isfinite(raw)) {
						continue;
					}

					const double xCenter = Asm.startPos.x + (i + 0.5) * dx;
					const double yCenter = Asm.startPos.y + (j + 0.5) * dy;
					const BenchmarkPinKey key{
						static_cast<long long>(std::llround(xCenter * 1.0e9)),
						static_cast<long long>(std::llround(yCenter * 1.0e9))
					};
					auto& pin = pins[key];
					pin.xCenter = xCenter;
					pin.yCenter = yCenter;
					pin.sliceRaw[sliceIndex] += raw;
					pin.coreRaw += raw;
				}
			}
		}
	}

	double totalRaw = 0.0;
	for (const auto& entry : pins) {
		totalRaw += entry.second.coreRaw;
	}
	if (pins.empty() || !(totalRaw > 0.0)) {
		out << "\nPinPowerBenchmarkSummary\n";
		out << "Status no_aggregated_fuel_pin_tallies\n";
		return;
	}

	const double normalization = static_cast<double>(pins.size()) / totalRaw;
	std::array<BenchmarkPinPowerExtrema, 3> sliceExtrema{};
	BenchmarkPinPowerExtrema coreExtrema{};
	for (const auto& entry : pins) {
		const auto& pin = entry.second;
		for (int sliceIndex = 0; sliceIndex < 3; ++sliceIndex) {
			UpdateBenchmarkPinPowerExtrema(
				sliceExtrema[sliceIndex],
				pin.sliceRaw[sliceIndex],
				pin.xCenter,
				pin.yCenter
			);
		}
		UpdateBenchmarkPinPowerExtrema(coreExtrema, pin.coreRaw, pin.xCenter, pin.yCenter);
	}

	std::ostringstream summary;
	summary << std::setprecision(17);
	summary << "\nPinPowerBenchmarkSummary\n";
	summary << "Status ok\n";
	summary << "TallyQuantity pinwise_fission_rate_collision_estimator\n";
	summary << "FuelSelection geometry_material_mask\n";
	summary << "Normalization total_power_equals_fuel_pin_count\n";
	summary << "FuelPinCount " << pins.size() << "\n";
	summary << "FuelZBounds(cm) " << fuelZMin << " " << fuelZMax << "\n";
	summary << "RawTotalFissionTally " << totalRaw << "\n";
	summary << "NormalizationFactor " << normalization << "\n";
	for (int sliceIndex = 0; sliceIndex < 3; ++sliceIndex) {
		const double zLow = fuelZMin + sliceIndex * sliceHeight;
		const double zHigh = zLow + sliceHeight;
		const auto& extrema = sliceExtrema[sliceIndex];
		summary << "PinPowerSlice " << (sliceIndex + 1)
			<< " ZBounds(cm) " << zLow << " " << zHigh
			<< " MinimumRaw " << extrema.minRaw
			<< " MinimumNormalized " << extrema.minRaw * normalization
			<< " MinimumXY(cm) " << extrema.minX << " " << extrema.minY
			<< " MaximumRaw " << extrema.maxRaw
			<< " MaximumNormalized " << extrema.maxRaw * normalization
			<< " MaximumXY(cm) " << extrema.maxX << " " << extrema.maxY
			<< "\n";
	}
	summary << "PinPowerCore"
		<< " MinimumRaw " << coreExtrema.minRaw
		<< " MinimumNormalized " << coreExtrema.minRaw * normalization
		<< " MinimumXY(cm) " << coreExtrema.minX << " " << coreExtrema.minY
		<< " MaximumRaw " << coreExtrema.maxRaw
		<< " MaximumNormalized " << coreExtrema.maxRaw * normalization
		<< " MaximumXY(cm) " << coreExtrema.maxX << " " << coreExtrema.maxY
		<< "\n";
	summary << "CoreMinimumPinPower " << coreExtrema.minRaw * normalization
		<< " XY(cm) " << coreExtrema.minX << " " << coreExtrema.minY << "\n";
	summary << "CoreMaximumPinPower " << coreExtrema.maxRaw * normalization
		<< " XY(cm) " << coreExtrema.maxX << " " << coreExtrema.maxY << "\n";

	out << summary.str();
	std::cout << summary.str();
}

__global__ void ZeroTallyPincellsKernel(TallyPincell* pcs, int n);

void ResetCoreTallyOnDevice(const TallyC5G7Geometry& h_CoreTally, TallyAssembly* d_bufferTallyAssembly, const std::vector<TallyPincell*>& d_bufferTallyPincellVec);
