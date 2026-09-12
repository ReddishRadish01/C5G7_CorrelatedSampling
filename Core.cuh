#pragma once
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"

constexpr double meatOrModEps = 1.0e-4;

HD inline void updateBestDistance(double t, double& best, double epsT) {
	if (isfinite(t) && t > epsT && t < best) {
		best = t;
	}
}

HD inline unsigned int inflate2D(unsigned int n) {
	n = (n | (n << 8)) & 0x00FF00FFu;
	n = (n | (n << 4)) & 0x0F0F0F0Fu;
	n = (n | (n << 2)) & 0x33333333u;
	n = (n | (n << 1)) & 0x55555555u;
	return n;
}

HD inline unsigned int morton2D(unsigned int x, unsigned int y) {
	return inflate2D(x) | (inflate2D(y) << 1);
}

H inline void buildMortonXYRank(int xNum, int yNum, int* xyToMortonRank) {
	if (xNum <= 0 || yNum <= 0 || xyToMortonRank == nullptr) {
		return;
	}

	const int xyNum = xNum * yNum;
	std::vector<int> logicalXY(xyNum);
	for (int i = 0; i < xyNum; ++i) {
		logicalXY[i] = i;
	}

	std::sort(logicalXY.begin(), logicalXY.end(), [xNum](int a, int b) {
		const int ax = a % xNum;
		const int ay = a / xNum;
		const int bx = b % xNum;
		const int by = b / xNum;
		const unsigned int am = morton2D(static_cast<unsigned int>(ax), static_cast<unsigned int>(ay));
		const unsigned int bm = morton2D(static_cast<unsigned int>(bx), static_cast<unsigned int>(by));
		if (am != bm) {
			return am < bm;
		}
		return a < b;
	});

	for (int rank = 0; rank < xyNum; ++rank) {
		xyToMortonRank[logicalXY[rank]] = rank;
	}
}

HD inline int pincellStorageIndex(int x, int y, int z, int xNum, int yNum, const int* xyToMortonRank) {
	const int logicalXY = y * xNum + x;
	const int xyNum = xNum * yNum;
	const int rank = (xyToMortonRank != nullptr) ? xyToMortonRank[logicalXY] : logicalXY;
	return z * xyNum + rank;
}

enum class BoundaryCondition {
	Vacuum,
	Reflective
};

class Pincell {
public:
	double sideLength = 1.26;	// 1.26 cm. alias with height
	double height = 1.26;
	double radius = 0.54;	// 0.54 cm radius
	MatType meatType = MatType::Unknown;
	MatType modType = MatType::MOD;
	vec2 centerOffset = { 0.0, 0.0 };
	double axialMeatSplitFraction = 0.0;
	MatType axialUpperMeatType = MatType::Unknown;
	bool hasAxialMeatSplit = false;
	// This remains true after a GT cell is changed to CR, allowing a second
	// rod-tip position to be evaluated without constructing another geometry.
	bool isControlRodChannel = false;

	H Pincell() = default;

	H Pincell(double sideLength, double radius, double height = 0.0, MatType meat = MatType::Unknown, MatType mod = MatType::MOD, vec2 centerOffset = { 0.0, 0.0 })
		: sideLength(sideLength), radius(radius), height(height), meatType(meat), modType(mod), centerOffset(centerOffset),
		isControlRodChannel(meat == MatType::GT || meat == MatType::CR)
	{
		// exceptions: if radius is bigger than sideLength / 2 ? remove it
		if (2 * radius > sideLength) { this->radius = sideLength / 2.0; }
	}

	HD MatType meatOrMod(vec3 pincellLocalPos) const {
		// this is for moderator block
		if (radius <= 1.0e-8) { return modType; }

		// rest is for regular pincell
		vec2 center = { this->sideLength / 2.0 + centerOffset.x, this->sideLength / 2.0 + centerOffset.y };
		if ((pincellLocalPos.x - center.x) * (pincellLocalPos.x - center.x) + (pincellLocalPos.y - center.y) * (pincellLocalPos.y - center.y) <= (radius * radius)) {
			if (hasAxialMeatSplit && pincellLocalPos.z >= axialMeatSplitFraction * height) {
				return axialUpperMeatType;
			}
			return meatType;
		}
		else 
			return modType;
	}


	HD double DTC(vec3 localPos, Neutron& n, XSLibrary* XSLib, GnuAMCM& RNG);
	HD double DTS(vec3 localPos, Neutron& n);

	HD bool isInsidePin(vec3 pincellLocalPos) const {
		if (radius == 0.0) { return false; }

		// rest is for regular pincell
		vec2 center = { this->sideLength / 2.0 + centerOffset.x, this->sideLength / 2.0 + centerOffset.y };
		if ((pincellLocalPos.x - center.x) * (pincellLocalPos.x - center.x) + (pincellLocalPos.y - center.y) * (pincellLocalPos.y - center.y) >= (radius * radius))
			return false;
		else
			return true;
	}

	HD MatType materialAtRodDepth(vec3 pincellLocalPos, double globalZ,
		double coreTopZ, double rodDepth) const;
	HD double csDeltaOpticalThickness(vec3 localStart, double globalStartZ,
		const Neutron& n, double length, const CorrelatedSamplingConfig& config,
		XSLibrary* XSLib) const;
};

class Assembly {
public:
	vec3 startPos = { 0, 0, 0 };

	vec3 length = { 0, 0, 0 };

	int xNum = 0;
	int yNum = 0;
	int zNum = 0;
	Pincell* pinCells = nullptr;
	int* xyToMortonRank = nullptr;
	Pincell OOBPincell;

	H Assembly() = default;

	H Assembly(vec3 startPos, vec3 length, int x, int y, int z)
		: startPos(startPos), length(length), xNum(x), yNum(y), zNum(z)
	{
	}

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

		// 시발 material 하드코딩 씨발인데
		MatType modType = (mod >= 0 && mod <= 7) ? static_cast<MatType>(mod) : MatType::MOD;
		
		this->pinCells = new Pincell[numPincells];
		this->xyToMortonRank = new int[this->xNum * this->yNum];
		buildMortonXYRank(this->xNum, this->yNum, this->xyToMortonRank);

		for (int i = 0; i < xNum * yNum; i++) {
			int meat = 0;
			assembly >> meat;
			MatType meatType = (meat >= 0 && meat <= 7) ? static_cast<MatType>(meat) : MatType::MOD;
			int x = i % xNum;
			int y = i / xNum;
			for (int j = 0; j < zNum; j++) {
				int index = pincellStorageIndex(x, y, j, xNum, yNum, this->xyToMortonRank);
				this->pinCells[index] = Pincell(pincellLength, radius, pincellHeight, meatType, modType);
				//std::cout << static_cast<int>(this->pinCells[j].meatType) << " ";
			}
			//std::cout << 
		}

		assembly.close();

		OOBPincell = Pincell(0.0, 0.0, 0.0);
	}
	
	H void Initialize_MOD(double cellSize, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }, double cellHeight = 0.0, double pincellHeight = 0.0) {
		vec3 length = endPos - startPos;

		double height = pincellHeight;
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
			OOBPincell = Pincell(0.0, 0.0, 0.0);
			return;
		}

		this->xNum = static_cast<int>((length.x + 1.0e-9) / cellSize);
		this->yNum = static_cast<int>((length.y + 1.0e-9) / cellSize);
		this->zNum = static_cast<int>((length.z + 1.0e-9) / height);
		// fuck this floating point shits - the length.z is not enough - 

		int numPincells = this->xNum * this->yNum * this->zNum;
		this->pinCells = new Pincell[numPincells];
		this->xyToMortonRank = new int[this->xNum * this->yNum];
		buildMortonXYRank(this->xNum, this->yNum, this->xyToMortonRank);

		for (int i = 0; i < this->xNum * this->yNum; i++) {
			int x = i % this->xNum;
			int y = i / this->xNum;
			for (int j = 0; j < this->zNum; j++) {
				int index = pincellStorageIndex(x, y, j, this->xNum, this->yNum, this->xyToMortonRank);
				this->pinCells[index] = Pincell(cellSize, 0.0, height, MatType::MOD, MatType::MOD);
				//std::cout << static_cast<int>(this->pinCells[j].meatType) << " ";
			}
		}

	}

	H void ApplyTopControlRodInsertion(double insertedDepth, MatType guideType = MatType::GT, MatType rodType = MatType::CR) {
		if (this->pinCells == nullptr || this->xNum <= 0 || this->yNum <= 0 || this->zNum <= 0) {
			return;
		}

		const double cellHeight = this->pinCells[0].height;
		if (cellHeight <= 0.0) {
			return;
		}

		const double activeHeight = cellHeight * this->zNum;
		const double clampedDepth = std::max(0.0, std::min(insertedDepth, activeHeight));
		const double rodBottomZ = activeHeight - clampedDepth;
		const double boundaryEps = std::max(1.0, activeHeight) * 1.0e-12;

		for (int z = 0; z < this->zNum; ++z) {
			const double cellBottomZ = static_cast<double>(z) * cellHeight;
			const double cellTopZ = cellBottomZ + cellHeight;
			const bool fullyGuide = rodBottomZ >= cellTopZ - boundaryEps;
			const bool fullyRodded = rodBottomZ <= cellBottomZ + boundaryEps;
			const double splitFraction = std::max(0.0, std::min(1.0, (rodBottomZ - cellBottomZ) / cellHeight));
			for (int y = 0; y < this->yNum; ++y) {
				for (int x = 0; x < this->xNum; ++x) {
					const int index = pincellStorageIndex(x, y, z, this->xNum, this->yNum, this->xyToMortonRank);
					Pincell& pinCell = this->pinCells[index];
					if (pinCell.meatType != guideType && pinCell.meatType != rodType) {
						continue;
					}
					pinCell.isControlRodChannel = true;

					pinCell.hasAxialMeatSplit = false;
					pinCell.axialMeatSplitFraction = 0.0;
					pinCell.axialUpperMeatType = MatType::Unknown;

					if (fullyGuide) {
						pinCell.meatType = guideType;
					}
					else if (fullyRodded) {
						pinCell.meatType = rodType;
					}
					else {
						pinCell.meatType = guideType;
						pinCell.hasAxialMeatSplit = true;
						pinCell.axialMeatSplitFraction = splitFraction;
						pinCell.axialUpperMeatType = rodType;
					}
				}
			}
		}
	}

	H void MakePincellStorageIndependent() {
		if (this->pinCells == nullptr) {
			return;
		}

		const int count = this->totalPincellNo();
		Pincell* independentPinCells = new Pincell[count];
		std::copy(this->pinCells, this->pinCells + count, independentPinCells);
		this->pinCells = independentPinCells;
	}

	H void ApplyTopControlRodInsertionFromGlobalTop(double globalTopZ, double insertedDepth,
		MatType guideType = MatType::GT, MatType rodType = MatType::CR) {
		if (insertedDepth < 0.0 || this->length.z <= 0.0) {
			return;
		}

		const double globalRodBottomZ = globalTopZ - insertedDepth;
		const double assemblyTopZ = this->startPos.z + this->length.z;
		const double localInsertedDepth = std::max(0.0, std::min(this->length.z, assemblyTopZ - globalRodBottomZ));
		this->ApplyTopControlRodInsertion(localInsertedDepth, guideType, rodType);
	}
	

	HD int totalPincellNo();
	HD Pincell& returnPincellByIndex(int x, int y, int z);
	HD Pincell& returnPincellByPos(const Neutron& n);
	HD vec3 returnFlooredNeutronPosInPincell(const Neutron& n);
	HD double DTC(Neutron& n, XSLibrary* XSLib, GnuAMCM& RNG);
	HD double DTS(Neutron& n);

	H void SetBounds(vec3 newStartPos, vec3 newLength) {
		this->startPos = newStartPos;
		this->length = newLength;
	}
};

class C5G7Geometry {
public:
	double x; // cm, max coordinate
	double y;
	double z;
	double xMin = 0.0;
	double yMin = 0.0;
	double zMin = 0.0;

	//PinCell* pinCells = nullptr;
	Assembly* assembly = nullptr;
	Assembly nullAssembly{};

	int assemblyNo = 0;

	BoundaryCondition xMinBoundary = BoundaryCondition::Reflective;
	BoundaryCondition yMinBoundary = BoundaryCondition::Reflective;
	BoundaryCondition zMinBoundary = BoundaryCondition::Reflective;
	BoundaryCondition xMaxBoundary = BoundaryCondition::Vacuum;
	BoundaryCondition yMaxBoundary = BoundaryCondition::Vacuum;
	BoundaryCondition zMaxBoundary = BoundaryCondition::Vacuum;

	// I have to put something here --- ahhh

	C5G7Geometry()
		: x(0.0), y(0.0), z(0.0), xMin(0.0), yMin(0.0), zMin(0.0)
	{
	}

	// lazy initialization - we are not going to fully construct here. Instead, we will use the C5G7GeometryFactory to initialize.
	C5G7Geometry(std::string CoreProfile)
	{
		std::ifstream coreProfile(CoreProfile);
		x = 0.6426;
		y = 0.6426;
		z = 2.1420;
	}

	// force bump into the neutron's directing assembly - this is to avoid the neutron being on the boundary of two assemblies, 
	// which will cause problems in the returnAssemblyByNeutron function.
	HD Assembly& returnAssemblyByNeutron(Neutron& n) {
		double eps = 1.0e-13;
		
		
		//Neutron localN = n;
		//localN.pos = localN.pos + epsVec;
		//you push the neutron inside.
		//n.pos = n.pos + epsVec;

		// really need boundary management 

		for (int i = 0; i < this->assemblyNo; i++) {
			vec3 endPos = this->assembly[i].startPos + this->assembly[i].length;
			if (n.pos.x >= this->assembly[i].startPos.x && n.pos.x < endPos.x) {
				if (n.pos.y >= this->assembly[i].startPos.y && n.pos.y < endPos.y) {
					if (n.pos.z >= this->assembly[i].startPos.z && n.pos.z < endPos.z) {
						return this->assembly[i];
					}
				}
			}
		}



		// this usually means fucked up - 
		// we are never meant to pass the out-of-bounds neutrons, or nullified neutrons in this function.
		// we never want the code to flow into this far, in this function - idK why it ended up here, maybe put some debugger outputs just in case

		
		printf("Neutron stuck\t");
		double x = n.dirVec.x * eps;
		double y = n.dirVec.y * eps;
		double z = n.dirVec.z * eps;
		vec3 epsVec = { x, y, z };

		n.pos = n.pos + epsVec;

		for (int i = 0; i < this->assemblyNo; i++) {
			vec3 endPos = this->assembly[i].startPos + this->assembly[i].length;
			if (n.pos.x >= this->assembly[i].startPos.x && n.pos.x < endPos.x) {
				if (n.pos.y >= this->assembly[i].startPos.y && n.pos.y < endPos.y) {
					if (n.pos.z >= this->assembly[i].startPos.z && n.pos.z < endPos.z) {
						return this->assembly[i];
					}
				}
			}
		}

		printf("Neutron Out-Of-Bounds in Position: (%f, %f, %f)\n", n.pos.x, n.pos.y, n.pos.z);
		//n.Nullify();
		// return null assembly 
		return this->nullAssembly;
		
		
	}

	HD bool isVacuumBoundaryHit(vec3 pos, double eps) const {
		return (pos.x <= this->xMin + eps && xMinBoundary == BoundaryCondition::Vacuum) ||
			(pos.y <= this->yMin + eps && yMinBoundary == BoundaryCondition::Vacuum) ||
			(pos.z <= this->zMin + eps && zMinBoundary == BoundaryCondition::Vacuum) ||
			(pos.x >= this->x - eps && xMaxBoundary == BoundaryCondition::Vacuum) ||
			(pos.y >= this->y - eps && yMaxBoundary == BoundaryCondition::Vacuum) ||
			(pos.z >= this->z - eps && zMaxBoundary == BoundaryCondition::Vacuum);
	}

	HD bool isReflectiveBoundaryHit(vec3 pos, double eps) const {
		return (pos.x <= this->xMin + eps && xMinBoundary == BoundaryCondition::Reflective) ||
			(pos.y <= this->yMin + eps && yMinBoundary == BoundaryCondition::Reflective) ||
			(pos.z <= this->zMin + eps && zMinBoundary == BoundaryCondition::Reflective) ||
			(pos.x >= this->x - eps && xMaxBoundary == BoundaryCondition::Reflective) ||
			(pos.y >= this->y - eps && yMaxBoundary == BoundaryCondition::Reflective) ||
			(pos.z >= this->z - eps && zMaxBoundary == BoundaryCondition::Reflective);
	}

	HD Assembly& unstuckNeutron(Neutron& n) {
		double eps = 1.0e-12;
		//printf("Neutron stuck, passed to unstuckNeutron - modifying it\t");
		double xBump = 0.0;
		double yBump = 0.0;
		double zBump = 0.0;

		double bigEps = 1.0E-8;
		double localCoreX = n.pos.x - this->xMin;
		double localCoreY = n.pos.y - this->yMin;
		double localCoreZ = n.pos.z - this->zMin;
		if (fmod(localCoreX, 21.42) < eps || fmod(localCoreX, 21.42) > 21.42 - eps) {
			xBump = n.dirVec.x > 0 ? bigEps : -bigEps;
		}
		if (fmod(localCoreY, 21.42) < eps || fmod(localCoreY, 21.42) > 21.42 - eps) {
			yBump = n.dirVec.y > 0 ? bigEps : -bigEps;
		}

		if (fmod(localCoreZ, 1.26) < eps || fmod(localCoreZ, 1.26) > 1.26 - eps) {
			zBump = n.dirVec.z > 0 ? bigEps : -bigEps;
		}

		/*
		double x = n.dirVec.x > 0 ? n.dirVec.x * eps : n.dirVec.x * -eps;
		double y = n.dirVec.y > 0 ? n.dirVec.y * eps : n.dirVec.y * -eps;
		double z = n.dirVec.z > 0 ? n.dirVec.z * eps : n.dirVec.z * -eps;
		*/
		double x = n.dirVec.x * eps;
		double y = n.dirVec.y * eps;
		double z = n.dirVec.z * eps;

		vec3 epsVec = { x + xBump, y + yBump, z + zBump };

		n.pos = n.pos + epsVec;
		printf("Neutron stuck, passed to unstuckNeutron - stuck neutron bumped to (%f, %f, %f)\n", n.pos.x, n.pos.y, n.pos.z);

		for (int i = 0; i < this->assemblyNo; i++) {
			vec3 endPos = this->assembly[i].startPos + this->assembly[i].length;
			if (n.pos.x >= this->assembly[i].startPos.x && n.pos.x < endPos.x) {
				if (n.pos.y >= this->assembly[i].startPos.y && n.pos.y < endPos.y) {
					if (n.pos.z >= this->assembly[i].startPos.z && n.pos.z < endPos.z) {
						return this->assembly[i];
					}
				}
			}
		}

		printf("UnstuckNeutron Failed, neutron nullified at (%f,%f,%f)\n", n.pos.x, n.pos.y, n.pos.z);
		n.Nullify();
		return this->nullAssembly;
	}

};


H static void fuelLayoutDebug(Assembly& assembly) {
	std::cout << "\nAssembly layout info: startPos: (" << assembly.startPos.x << ", " << assembly.startPos.y << ", " << assembly.startPos.z << "), length: (" << assembly.length.x << ", " << assembly.length.y << ", " << assembly.length.z << "), xNum: " << assembly.xNum << ", yNum: " << assembly.yNum << ", zNum: " << assembly.zNum << "\n";
	const double trackingSideLength = (assembly.xNum > 0) ? assembly.length.x / assembly.xNum : assembly.pinCells[0].sideLength;
	const double trackingHeight = (assembly.zNum > 0) ? assembly.length.z / assembly.zNum : assembly.pinCells[0].height;
	std::cout << "Pincell info: sideLength: " << trackingSideLength << ", radius: " << assembly.pinCells[0].radius << ", height: " << trackingHeight << ", meatType: " << static_cast<int>(assembly.pinCells[0].meatType) << ", modType: " << static_cast<int>(assembly.pinCells[0].modType) << "\n";
	std::cout << "Assembly layout info: total number of pincells in this assembly: " << assembly.totalPincellNo() << "\n";
	std::cout << "for the z=0 layer:\n";
	for (int k = 0; k < 1; k++) {
		for (int ii = 0; ii < assembly.yNum; ii++) {
			for (int ij = 0; ij < assembly.xNum; ij++) {
				int index = pincellStorageIndex(ij, ii, k, assembly.xNum, assembly.yNum, assembly.xyToMortonRank);
				std::cout << static_cast<int>(assembly.pinCells[index].meatType) << " ";
			}
			std::cout << "\n";
		}
	}
}


class C5G7GeometryFactory {
public:
	H static void Initialize(C5G7Geometry& Core, std::string totalCoreProfile, std::string UO2Geometry, std::string MOXGeometry, std::string CRMODGeometry = "", double cellHeight = 0.0, double controlRodInsertionDepth = -1.0) {

		double pincellSize = 0.0;
		double pincellHeight = 0.0;
		
		std::ifstream totalCore(totalCoreProfile); 
		std::ifstream lineFetch(totalCoreProfile);
		std::string line;
		int lineNo = 0;
		while (std::getline(lineFetch, line)) {
			bool isBlank = std::all_of(line.begin(), line.end(),
				[](unsigned char c) { return std::isspace(c); });
			if (!isBlank) lineNo++;
		}
		lineFetch.close();

		Assembly asmUO2{};
		asmUO2.Initialize(UO2Geometry);
		Assembly asmMOX{};
		asmMOX.Initialize(MOXGeometry);
		Assembly asmCRUO2{};
		asmCRUO2.Initialize(UO2Geometry);
		asmCRUO2.ApplyTopControlRodInsertion(1.0e300);
		Assembly asmCRMOX{};
		asmCRMOX.Initialize(MOXGeometry);
		asmCRMOX.ApplyTopControlRodInsertion(1.0e300);
		Assembly asmCRMOD{};
		if (!CRMODGeometry.empty()) {
			asmCRMOD.Initialize(CRMODGeometry);
		}

		Assembly asmMOD{};
		std::cout << asmUO2.length.z / asmUO2.zNum << "\n";
		asmMOD.Initialize_MOD(pincellSize, { 0,0,0 }, { 0,0,0 }, asmUO2.length.z / asmUO2.zNum);
		

		totalCore >> Core.x >> Core.y >> Core.z;
		totalCore >> pincellSize >> pincellHeight;

		std::vector<Assembly> assemblyVec;
		assemblyVec.reserve(lineNo - 2);
		//std::string line;

		Core.assembly = new Assembly[lineNo - 2];

		for (int i = 0; i < lineNo - 2; i++) {
			Assembly assembly{};
			std::string assemblyType;
			bool supportsDynamicRodInsertion = false;

			totalCore >> assemblyType;

			double startPosX = 0.0; double endPosX = 0.0;
			double startPosY = 0.0; double endPosY = 0.0;
			double startPosZ = 0.0; double endPosZ = 0.0;
			totalCore >> startPosX >> endPosX >> startPosY >> endPosY >> startPosZ >> endPosZ;

			double lengthX = endPosX - startPosX;
			double lengthY = endPosY - startPosY;
			double lengthZ = endPosZ - startPosZ;
			vec3 startPos{startPosX, startPosY, startPosZ};
			vec3 endPos(endPosX, endPosY, endPosZ);
			vec3 length{lengthX, lengthY, lengthZ};

			if (assemblyType == "UO2") {
				assemblyVec.emplace_back(asmUO2);
				assemblyVec[i].SetBounds(startPos, length);
				supportsDynamicRodInsertion = true;
			}
			else if (assemblyType == "CRUO2") {
				assemblyVec.emplace_back(asmCRUO2);
				assemblyVec[i].SetBounds(startPos, length);
			}
			else if (assemblyType == "MOX") {
				assemblyVec.emplace_back(asmMOX);
				assemblyVec[i].SetBounds(startPos, length);
				supportsDynamicRodInsertion = true;
			}
			else if (assemblyType == "CRMOX") {
				assemblyVec.emplace_back(asmCRMOX);
				assemblyVec[i].SetBounds(startPos, length);
			}
			else if (assemblyType == "CRMOD") {
				if (asmCRMOD.pinCells == nullptr) {
					throw std::runtime_error("Core profile uses CRMOD, but no CRMOD geometry file was provided.");
				}
				assemblyVec.emplace_back(asmCRMOD);
				assemblyVec[i].SetBounds(startPos, length);
				supportsDynamicRodInsertion = true;
			}
			else if (assemblyType == "MOD") {
				assemblyVec.emplace_back(asmMOD);
				assemblyVec[i].Initialize_MOD(pincellSize, startPos, endPos, pincellHeight);
			}
			else {
				throw std::runtime_error("Unknown assembly type in core profile: " + assemblyType);
			}

			if (controlRodInsertionDepth >= 0.0 && supportsDynamicRodInsertion) {
				assemblyVec[i].MakePincellStorageIndependent();
				assemblyVec[i].ApplyTopControlRodInsertionFromGlobalTop(Core.z, controlRodInsertionDepth);
			}
			fuelLayoutDebug(assemblyVec[i]);

		} 
		//std::cout << "line size: " << lineNo - 2 << ", vector size: " << assemblyVec.size() << "\n";

		
		if (lineNo - 2 != assemblyVec.size()) {
			delete[] Core.assembly;
			Core.assembly = new Assembly[assemblyVec.size()];
		}

		std::copy(assemblyVec.begin(), assemblyVec.end(), Core.assembly);

		Core.assemblyNo = static_cast<int>(assemblyVec.size());
		if (!assemblyVec.empty()) {
			Core.xMin = assemblyVec[0].startPos.x;
			Core.yMin = assemblyVec[0].startPos.y;
			Core.zMin = assemblyVec[0].startPos.z;
			Core.x = assemblyVec[0].startPos.x + assemblyVec[0].length.x;
			Core.y = assemblyVec[0].startPos.y + assemblyVec[0].length.y;
			Core.z = assemblyVec[0].startPos.z + assemblyVec[0].length.z;

			for (const Assembly& assembly : assemblyVec) {
				vec3 endPos = assembly.startPos + assembly.length;
				Core.xMin = std::min(Core.xMin, assembly.startPos.x);
				Core.yMin = std::min(Core.yMin, assembly.startPos.y);
				Core.zMin = std::min(Core.zMin, assembly.startPos.z);
				Core.x = std::max(Core.x, endPos.x);
				Core.y = std::max(Core.y, endPos.y);
				Core.z = std::max(Core.z, endPos.z);
			}
		}

		totalCore.close();
	}
	
};
