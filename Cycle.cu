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

//#define NEUTRONIDXDEBUG
//#define REFLECTIONDEBUG
#define OUTBOUNDDEBUG

//#define TALLYINCYCLE

#if defined(TALLYINCYCLE) && \
	(defined(C5G7_TALLY_PINWISE_FISSION_RATE) || defined(C5G7_TALLY_FISSION_SOURCE_DIST))
#error TALLYINCYCLE cannot be combined with a fission tally mode.
#endif


constexpr double globaleps = 1.0e-10;

HD inline double sign_nonzero(double v) {
    return (v >= 0.0) ? 1.0 : -1.0;
}

HD inline bool near_multiple(double x, double pitch, double tol) {
    if (pitch <= 0.0) return false;

    double r = fmod(x, pitch);
    if (r < 0.0) r += pitch;

    return (r < tol) || ((pitch - r) < tol);
}

HD inline bool bumpAcrossFlatCellPlane(Neutron& n, const Assembly& a, double tol) {
    bool bumped = false;

    vec3 local = n.pos - a.startPos;

    double dx = a.length.x / a.xNum;
    double dy = a.length.y / a.yNum;
    double dz = a.length.z / a.zNum;

    if (fabs(n.dirVec.x) > 1.0e-14 &&
        local.x > tol && local.x < a.length.x - tol &&
        near_multiple(local.x, dx, tol)) {
        n.pos.x += sign_nonzero(n.dirVec.x) * tol;
        bumped = true;
    }

    if (fabs(n.dirVec.y) > 1.0e-14 &&
        local.y > tol && local.y < a.length.y - tol &&
        near_multiple(local.y, dy, tol)) {
        n.pos.y += sign_nonzero(n.dirVec.y) * tol;
        bumped = true;
    }

    if (fabs(n.dirVec.z) > 1.0e-14 &&
        local.z > tol && local.z < a.length.z - tol &&
        near_multiple(local.z, dz, tol)) {
        n.pos.z += sign_nonzero(n.dirVec.z) * tol;
        bumped = true;
    }

    return bumped;
}


G void globalNeutronWeightOffset(NeutronBank* Bank, double* weightModifier) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx >= Bank->allocatableNeutronNum) {
		return;
	}
	if (!Bank->neutrons[idx].isNullified()) {
		Bank->neutrons[idx].weight *= *weightModifier;
	}
	if (!Bank->addedNeutrons[idx].isNullified()) {
		Bank->addedNeutrons[idx].weight *= *weightModifier;
	}
}


__inline__ __device__
double warpReduceSum(double val) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

__inline__ __device__
double blockReduceSum(double val) {
    static __shared__ double shared[32]; // max 1024 threads / 32 warps

    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;

    val = warpReduceSum(val);

    if (lane == 0) {
        shared[wid] = val;
    }

    __syncthreads();

    val = 0.0;

    if (threadIdx.x < blockDim.x / warpSize) {
        val = shared[lane];
    }

    if (wid == 0) {
        val = warpReduceSum(val);
    }

    return val;
}

G void fetchNeutronNumber(
    NeutronBank* Bank,
    double* d_NeutronNum,
    double* d_NeutronWeightSum
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    double localNum = 0.0;
    double localWeight = 0.0;

    if (idx < Bank->allocatableNeutronNum) {
        if (!Bank->neutrons[idx].isNullified()) {
            localNum += 1.0;
            localWeight += Bank->neutrons[idx].weight;
        }

        if (!Bank->addedNeutrons[idx].isNullified()) {
            localNum += 1.0;
            localWeight += Bank->addedNeutrons[idx].weight;
        }
    }

    double blockNum = blockReduceSum(localNum);
    double blockWeight = blockReduceSum(localWeight);

    if (threadIdx.x == 0) {
        atomicAdd(d_NeutronNum, blockNum);
        atomicAdd(d_NeutronWeightSum, blockWeight);
    }
}

/*
// Naive atomicAdd approach - not efficient
G void fetchNeutronNumber(NeutronBank* Bank, double* d_NeutronNum, double* d_NeutronWeightSum) {
    // here maybe we can use parallel reduction?

    //int tid = threadIdx.x;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }

    if (!Bank->neutrons[idx].isNullified()) {
        atomicAdd(d_NeutronNum, 1.0);
        atomicAdd(d_NeutronWeightSum, Bank->neutrons[idx].weight);
    }

    if (!Bank->addedNeutrons[idx].isNullified()) {
        atomicAdd(d_NeutronNum, 1.0);
        atomicAdd(d_NeutronWeightSum, Bank->addedNeutrons[idx].weight);
    }
}
*/

G void weightParallelReductionSum(NeutronBank* Bank, double*) {
    
}


G void cycle_Neutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes, double* fissionWeight,
    CorrelatedSamplingConfig csConfig, double* csFissionScores, double* csSourceSums,
    double* csSourceCount, unsigned long long* csSupportMismatch,
    CsWeightDiagnostics* csWeightDiagnostics) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }
    GnuAMCM RNG(seedNo[idx]);

    double eps = 1.0e-13;
    int reflectionCounter = 0;

    if (!Bank->neutrons[idx].isNullified()) {
        for (int i = 0; i < iterTimes; i++) {


            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->neutrons[idx]);
            // currently, currentAssembly pushes the neutron inside the assembly by its dirvec direction
            double fissionSourceScore = 0.0;
            double fissionRateScore = 0.0;
#ifdef TALLYINCYCLE
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->neutrons[idx]);
#endif
#if defined(C5G7_TALLY_PINWISE_FISSION_RATE) || defined(C5G7_TALLY_FISSION_SOURCE_DIST)
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->neutrons[idx]);
#endif


            if (Bank->neutrons[idx].isNullified()) {
                atomicAdd(&(Bank->neutronSize), -1);
#ifdef OUTBOUNDDEBUG
                printf("Neutron idx %d nullfied because Assembly positioning got fucked\n", idx);
#endif
                seedNo[idx] = RNG.gen();
                return;
            }

            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->neutrons[idx]);
            vec3 flooredNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);
#ifdef TALLYINCYCLE
            TallyPincell currentTallyPincell = currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]);
#endif
#if defined(C5G7_TALLY_PINWISE_FISSION_RATE) || defined(C5G7_TALLY_FISSION_SOURCE_DIST)
            TallyPincell currentTallyPincell = currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]);
#endif

            // handling OOB pincells returned from Assembly indexing error - OOB pincell returned from function returnPincellByPos()
            if (currentPincell.sideLength < 1e-4 || currentPincell.height < 1e-4) {
                // on-surface exceptions - bump it into other cell, and continue
                // this is error exception control - almost no neutron ends up in this branch.
                /*
                if ((Bank->neutrons[idx].pos.x >= currentAssembly.startPos.x + currentAssembly.length.x - eps * 100) && (Bank->neutrons[idx].pos.x <= currentAssembly.startPos.x + currentAssembly.length.x + eps * 100) ||
                    (Bank->neutrons[idx].pos.y >= currentAssembly.startPos.y + currentAssembly.length.y - eps * 100) && (Bank->neutrons[idx].pos.y <= currentAssembly.startPos.y + currentAssembly.length.y + eps * 100) ||
                    (Bank->neutrons[idx].pos.z >= currentAssembly.startPos.z + currentAssembly.length.z - eps * 100) && (Bank->neutrons[idx].pos.z <= currentAssembly.startPos.z + currentAssembly.length.z + eps * 100)) {
                    Bank->neutrons[idx].updateWithLength(eps * 10000);
                    //printf("lmao you leak here");
                    continue;
                }
                
                printf("error in returning pincell of neutron idx %d, - nullifying this neutron\n", idx);
                Bank->neutrons[idx].Nullify();
                atomicAdd(&(Bank->neutronSize), -1);
                */

				currentAssembly = Core->unstuckNeutron(Bank->neutrons[idx]);
				currentPincell = currentAssembly.returnPincellByPos(Bank->neutrons[idx]);
                flooredNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);

                /*
                seedNo[idx] = RNG.gen();
                return;
                */
            }


            double DTC = currentAssembly.DTC(Bank->neutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->neutrons[idx]);
			const double csSegmentLength = DTC < DTS ? DTC : DTS;
			if (csConfig.enabled) {
				Bank->neutrons[idx].csLogLikelihood -= currentPincell.csDeltaOpticalThickness(
					flooredNeutronPos, Bank->neutrons[idx].pos.z, Bank->neutrons[idx],
					csSegmentLength, csConfig, XSLib);
			}

            if (DTC < DTS) {    // reaction
                Bank->neutrons[idx].updateWithLength(DTC);
				const vec3 collisionLocalPos = flooredNeutronPos + Bank->neutrons[idx].dirVec * DTC;
				const MatType perturbedMaterial = csConfig.enabled
					? currentPincell.materialAtRodDepth(collisionLocalPos, Bank->neutrons[idx].pos.z,
						csConfig.coreTopZ, csConfig.perturbedRodDepth)
					: currentPincell.meatOrMod(flooredNeutronPos);
                //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);

                // implicit reaction - in order to accomodate neutron weight and some RR shits
                Interaction::reaction_implicit(
                    Bank->neutrons[idx],
                    Bank,
                    XSLib,
                    currentPincell,
                    flooredNeutronPos,
                    RNG,
                    k_mult,
                    passFlag,
                    false,
                    fissionWeight,
                    fissionSourceScore,
                    fissionRateScore,
					perturbedMaterial,
					csConfig,
					csFissionScores,
					csSourceSums,
					csSourceCount,
					csSupportMismatch,
					csWeightDiagnostics,
					(idx % CS_WEIGHT_HIST_SAMPLE_STRIDE) == 0
                );

#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), DTC); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), DTC); }
#endif

#ifdef C5G7_TALLY_PINWISE_FISSION_RATE
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), fissionRateScore); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), fissionRateScore); }
#elif defined(C5G7_TALLY_FISSION_SOURCE_DIST)
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), fissionSourceScore); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), fissionSourceScore); }
#endif

                Interaction::russianRoulette(Bank->neutrons[idx], Bank, RNG, 0.25, 0.5, false);
                if (Bank->neutrons[idx].isNullified()) {
                    seedNo[idx] = RNG.gen();
                    return;
                }

            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                //vec3 updatedPos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTC;
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), DTS); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), DTS); }
#endif
                vec3 afterDTSPos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTS;
                vec3 afterDTCPos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTC;

                // ** Neutrons leaked through vaccum boundary
                
                if (Core->isVacuumBoundaryHit(afterDTSPos, eps * 100)) {
                    Bank->neutrons[idx].Nullify();
                    atomicAdd(&(Bank->neutronSize), -1);
                    seedNo[idx] = RNG.gen();
                    return;
                }

                
                // ** Neutrons reflected at boundaries - push neutron to the surface, bump it (with reflection function) and continue the loop
                //if (afterDTSPos.x < eps || afterDTSPos.y < eps || afterDTSPos.z < eps || afterDTSPos.x >= Core->x - eps || afterDTSPos.y >= Core->y - eps || afterDTSPos.z >= Core->z - eps) {
                if (Core->isReflectiveBoundaryHit(afterDTSPos, eps * 100)) {
                    if (true) {
                        //Bank->neutrons[idx].printInfo_Kernel(idx);
                        Interaction::reflection(Bank->neutrons[idx], Core, DTS, DTC, afterDTSPos, eps * 100);
#ifdef REFLECTIONDEBUG
                        reflectionCounter++;
                        if (reflectionCounter > 1) {
                            printf("Neutron reflected twice on ");
                            Bank->neutrons[idx].printInfo_Kernel(idx);
                        }
#endif
                        // this reflection already bumps the neutron to the surface location, with some room for error
                        if (Bank->neutrons[idx].isNullified()) {
                            // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                            atomicAdd(&(Bank->neutronSize), -1);
                            seedNo[idx] = RNG.gen();
                            return;
                        }
                        //continue;
                    }
                }
                else {
                    // Pin-moderator trans boundary situation - seems like neutron gets stuck with DTS of zero - thus even with eps push at the line above, it still gets stuck.
                    if (DTS < eps * 10) {
                        const double boundaryBump = 1.0e-8;

                        bool bumped = bumpAcrossFlatCellPlane(Bank->neutrons[idx], currentAssembly, boundaryBump);

                        if (!bumped) {
                            // fallback for pin/mod cylindrical boundary or any unresolved near-surface case
                            Bank->neutrons[idx].updateWithLength(boundaryBump);
                        }

                        continue;
                                              
                        /*
                        double localX = fmod(afterDTSPos.x, 1.26);
                        double localY = fmod(afterDTSPos.y, 1.26);
                        double distSq = (localX - 0.63) * (localX - 0.63) + (localY - 0.63) * (localY - 0.63);
                        double radiusSq = 0.54 * 0.54;
                        if (distSq < radiusSq + eps * 100 && distSq > radiusSq - eps * 100) {
                            //printf("Pushed neutron at the rod-moderator boundary, because DTS was so small\n");
                            //Bank->neutrons[idx].printInfo_Kernel(idx);
                            Bank->neutrons[idx].updateWithLength(eps * 1.0e+4);
                            //Bank->neutrons[idx].printInfo_Kernel(idx);
                        }
                        */
                    }
                    else {
                        // neutron is inside boundary, no reflection, and entered different cell(Trans-surface)
                        Bank->neutrons[idx].updateWithLength(DTS + eps * 100);
                    }

                }
                
                
            }
            //printf("idx %d neutron on %d loop\n", idx, i);

            if (i == iterTimes * 7 / 10 && DTS < eps * 10) {
                printf("idx %d neutron on %d loop struggling: (%f,%f,%f), (%f,%f,%f), DTC:%e, DTS:%e\n", idx, i, Bank->neutrons[idx].pos.x, Bank->neutrons[idx].pos.y, Bank->neutrons[idx].pos.z, Bank->neutrons[idx].dirVec.x, Bank->neutrons[idx].dirVec.y, Bank->neutrons[idx].dirVec.z, DTC, DTS);
                Bank->neutrons[idx].updateWithLength(eps * 1.0e+3);
            }

            if (i == iterTimes - 1) {
                //Bank->neutrons[idx].printInfo_Kernel(idx);
                //printf("idx %d neutron (%f,%f,%f) didn't reacted after %d loops..? weight: %f\n", idx, Bank->neutrons[idx].pos.x, Bank->neutrons[idx].pos.y, Bank->neutrons[idx].pos.z, iterTimes, Bank->neutrons[idx].weight);
                //Bank->neutrons[idx].updateWithLength(1.0e-7);

                printf("STUCK neutron idx %d: pos=(%.17g,%.17g,%.17g), dir=(%.17g,%.17g,%.17g), weight=%.17g, DTS:%.17g \n",
                    idx,
                    Bank->neutrons[idx].pos.x, Bank->neutrons[idx].pos.y, Bank->neutrons[idx].pos.z,
                    Bank->neutrons[idx].dirVec.x, Bank->neutrons[idx].dirVec.y, Bank->neutrons[idx].dirVec.z,
                    Bank->neutrons[idx].weight, DTS);

                Bank->neutrons[idx].updateWithLength(1.0e-7);
                seedNo[idx] = RNG.gen();
                return;
                
                //vec3 bump = { 1.0e-3, 1.0e-3, 1.0e-3 };
                //vec3 randPos = { RNG.uniform(0, Core->x), RNG.uniform(0, Core->y), RNG.uniform(0, Core->z) };
                //Bank->neutrons[idx].pos = Bank->neutrons[idx].pos - bump;
                //Bank->neutrons[idx].pos = randPos;
                //Bank->neutrons[idx].dirVec = vec3::randomUnit(RNG);
                //seedNo[idx] = RNG.gen();
                //Bank->neutrons[idx].Nullify();
                //atomicAdd(&(Bank->neutronSize), -1);
                //return;
            }
        }


        //printf("idx %d neutron didn't reacted after %d loops..?\n", idx, iterTimes);
        seedNo[idx] = RNG.gen();
        return;

    }
    seedNo[idx] = RNG.gen();
    return;

}

G void addedNeutronPassResetter(NeutronBank* Neutrons) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= Neutrons->allocatableNeutronNum) { return; }

    if (!Neutrons->addedNeutrons[idx].isNullified()) {
        Neutrons->addedNeutrons[idx].passFlag = false;
    }
}


G void cycle_addedNeutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes, double* fissionWeight,
    CorrelatedSamplingConfig csConfig, double* csFissionScores, double* csSourceSums,
    double* csSourceCount, unsigned long long* csSupportMismatch,
    CsWeightDiagnostics* csWeightDiagnostics) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }

    GnuAMCM RNG(seedNo[idx]);
    double eps = 1.0e-13;
    int reflectionCounter = 0;

    if (!Bank->addedNeutrons[idx].isNullified()) {
        if (Bank->addedNeutrons[idx].passFlag) {
            Bank->addedNeutrons[idx].passFlag = false;
            return;
        }
    }

    if (!Bank->addedNeutrons[idx].isNullified()) {
        for (int i = 0; i < iterTimes; i++) {


            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
            double fissionSourceScore = 0.0;
            double fissionRateScore = 0.0;
#ifdef TALLYINCYCLE
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
#endif
#if defined(C5G7_TALLY_PINWISE_FISSION_RATE) || defined(C5G7_TALLY_FISSION_SOURCE_DIST)
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
#endif

            if (Bank->addedNeutrons[idx].isNullified()) {
                atomicAdd(&(Bank->addedNeutronSize), -1);
#ifdef OUTBOUNDDEBUG
                printf("Neutron idx %d nullfied because Assembly positioning got fucked\n", idx);
#endif
                seedNo[idx] = RNG.gen();
                return;
            }

            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
            vec3 flooredAddedNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);
                
#ifdef TALLYINCYCLE
            TallyPincell currentTallyPincell = currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
#endif
#if defined(C5G7_TALLY_PINWISE_FISSION_RATE) || defined(C5G7_TALLY_FISSION_SOURCE_DIST)
            TallyPincell currentTallyPincell = currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
#endif

            // OOB Pincell handling
            if (currentPincell.sideLength < 1e-4 || currentPincell.height < 1e-4) {
                /*
                if ((Bank->addedNeutrons[idx].pos.x >= currentAssembly.startPos.x + currentAssembly.length.x - eps * 100) && (Bank->addedNeutrons[idx].pos.x <= currentAssembly.startPos.x + currentAssembly.length.x + eps * 100) ||
                    (Bank->addedNeutrons[idx].pos.y >= currentAssembly.startPos.y + currentAssembly.length.y - eps * 100) && (Bank->addedNeutrons[idx].pos.y <= currentAssembly.startPos.y + currentAssembly.length.y + eps * 100) ||
                    (Bank->addedNeutrons[idx].pos.z >= currentAssembly.startPos.z + currentAssembly.length.z - eps * 100) && (Bank->addedNeutrons[idx].pos.z <= currentAssembly.startPos.z + currentAssembly.length.z + eps * 100)) {
                    Bank->addedNeutrons[idx].updateWithLength(eps * 10000);
                    continue;
                    printf("Lmao you leak here");
                }



                printf("error in returning pincell of adddedNeutron idx %d, - nullifying this neutron\n", idx);
                Bank->addedNeutrons[idx].Nullify();
                atomicAdd(&(Bank->addedNeutronSize), -1);
                seedNo[idx] = RNG.gen();
                return;
                */
                currentAssembly = Core->unstuckNeutron(Bank->addedNeutrons[idx]);
                currentPincell = currentAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
                flooredAddedNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);
            }

            double DTC = currentAssembly.DTC(Bank->addedNeutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->addedNeutrons[idx]);
			const double csSegmentLength = DTC < DTS ? DTC : DTS;
			if (csConfig.enabled) {
				Bank->addedNeutrons[idx].csLogLikelihood -= currentPincell.csDeltaOpticalThickness(
					flooredAddedNeutronPos, Bank->addedNeutrons[idx].pos.z, Bank->addedNeutrons[idx],
					csSegmentLength, csConfig, XSLib);
			}

            if (DTC < DTS) {    // reaction
                Bank->addedNeutrons[idx].updateWithLength(DTC);
				const vec3 collisionLocalPos = flooredAddedNeutronPos + Bank->addedNeutrons[idx].dirVec * DTC;
				const MatType perturbedMaterial = csConfig.enabled
					? currentPincell.materialAtRodDepth(collisionLocalPos, Bank->addedNeutrons[idx].pos.z,
						csConfig.coreTopZ, csConfig.perturbedRodDepth)
					: currentPincell.meatOrMod(flooredAddedNeutronPos);
                //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);

                Interaction::reaction_implicit(
                    Bank->addedNeutrons[idx],
                    Bank,
                    XSLib,
                    currentPincell,
                    flooredAddedNeutronPos,
                    RNG,
                    k_mult,
                    passFlag,
                    true,
                    fissionWeight,
                    fissionSourceScore,
                    fissionRateScore,
					perturbedMaterial,
					csConfig,
					csFissionScores,
					csSourceSums,
					csSourceCount,
					csSupportMismatch,
					csWeightDiagnostics,
					(idx % CS_WEIGHT_HIST_SAMPLE_STRIDE) == 0
                );

#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), DTC); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), DTC); }
#endif
#ifdef C5G7_TALLY_PINWISE_FISSION_RATE
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), fissionRateScore); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), fissionRateScore); }
#elif defined(C5G7_TALLY_FISSION_SOURCE_DIST)
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), fissionSourceScore); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), fissionSourceScore); }
#endif

                Interaction::russianRoulette(Bank->addedNeutrons[idx], Bank, RNG, 0.25, 0.5, true);
                if (Bank->addedNeutrons[idx].isNullified()) {
                    seedNo[idx] = RNG.gen();
                    return;
                }
            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                //vec3 updatedPos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTC;
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), DTS); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), DTS); }
#endif
                vec3 afterDTSPos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTS;
                vec3 afterDTCPos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTC;

                // handle vaccum boundary neutrons;
                //if (afterDTSPos.x >= Core->x || afterDTSPos.y >= Core->y || afterDTSPos.z >= Core->z) {
                
                if (Core->isVacuumBoundaryHit(afterDTSPos, eps * 100)) {
                    Bank->addedNeutrons[idx].Nullify();
                    atomicAdd(&(Bank->addedNeutronSize), -1);
                    seedNo[idx] = RNG.gen();
                    return;
                }
                

                // Reflection
                //if (afterDTSPos.x < eps || afterDTSPos.y < eps || afterDTSPos.z < eps || afterDTSPos.x >= Core->x - eps || afterDTSPos.y >= Core->y - eps || afterDTSPos.z >= Core->z - eps) {
                if (Core->isReflectiveBoundaryHit(afterDTSPos, eps * 100)) {
                    if (true) {
						//Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                        Interaction::reflection(Bank->addedNeutrons[idx], Core, DTS, DTC, afterDTSPos, eps * 100);
#ifdef REFLECTIONDEBUG
                        if (reflectionCounter > 1) {
                            printf("Neutron reflected twice on ");
                            Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                        }
#endif
                        if (Bank->addedNeutrons[idx].isNullified()) {
                            // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                            atomicAdd(&(Bank->addedNeutronSize), -1);
                            seedNo[idx] = RNG.gen();
                            return;
                        }
                    }
                }
                else {
                    // Pin-moderator trans boundary situation - seems like neutron gets stuck with DTS of zero - thus even with eps push at the line above, it still gets stuck.
                    if (DTS < eps * 10) {
                        /*
						double localX = fmod(afterDTSPos.x, 1.26);
                        double localY = fmod(afterDTSPos.y, 1.26);
						double distSq = (localX - 0.63) * (localX - 0.63) + (localY - 0.63) * (localY - 0.63);
                        double radiusSq = 0.54 * 0.54;
                        if (distSq < radiusSq + eps * 100 && distSq > radiusSq - eps * 100) {
                            printf("Pushed neutron at the rod-moderator boundary, because DTS was so small\n");
                            //Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                            Bank->addedNeutrons[idx].updateWithLength(eps * 1.0e+4);
                            //Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                        }
                        */
                        const double boundaryBump = 1.0e-8;

                        bool bumped = bumpAcrossFlatCellPlane(Bank->addedNeutrons[idx], currentAssembly, boundaryBump);

                        if (!bumped) {
                            // fallback for pin/mod cylindrical boundary or any unresolved near-surface case
                            Bank->addedNeutrons[idx].updateWithLength(boundaryBump);
                        }

                        continue;
                    }
                    else {
                        // neutron is inside boundary, no reflection, and entered different cell(Trans-surface)
                        Bank->addedNeutrons[idx].updateWithLength(DTS + eps * 100);
                    }
                }
                
                
            }
            //printf("idx %d neutron on %d loop\n", idx, i);


            if (i == iterTimes * 6 / 10 && DTS < eps * 10) {
                printf("idx %d neutron on %d loop struggling: (%f,%f,%f), (%f,%f,%f), DTC:%e, DTS:%e\n", idx, i, Bank->addedNeutrons[idx].pos.x, Bank->addedNeutrons[idx].pos.y, Bank->addedNeutrons[idx].pos.z, Bank->addedNeutrons[idx].dirVec.x, Bank->addedNeutrons[idx].dirVec.y, Bank->addedNeutrons[idx].dirVec.z, DTC, DTS);
                Bank->addedNeutrons[idx].updateWithLength(eps * 1.0e+4);
            }

            if (i == iterTimes - 1) {
                //Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                //printf("idx %d addedNeutron (%f,%f,%f) didn't reacted after %d loops..? weight: %f\n", idx, Bank->addedNeutrons[idx].pos.x, Bank->addedNeutrons[idx].pos.y, Bank->addedNeutrons[idx].pos.z, iterTimes, Bank->addedNeutrons[idx].weight);
                //Bank->addedNeutrons[idx].updateWithLength(1.0e-4);

                printf("STUCK addedNeutron idx %d: pos=(%.17g,%.17g,%.17g), dir=(%.17g,%.17g,%.17g), weight=%.17g\n",
                    idx,
                    Bank->addedNeutrons[idx].pos.x, Bank->addedNeutrons[idx].pos.y, Bank->addedNeutrons[idx].pos.z,
                    Bank->addedNeutrons[idx].dirVec.x, Bank->addedNeutrons[idx].dirVec.y, Bank->addedNeutrons[idx].dirVec.z,
                    Bank->addedNeutrons[idx].weight);

                Bank->addedNeutrons[idx].updateWithLength(1.0e-8);
                seedNo[idx] = RNG.gen();
                return;
                
                //vec3 bump = { 1.0e-3, 1.0e-3, 1.0e-3 };
                //vec3 randPos = { RNG.uniform(0, Core->x), RNG.uniform(0, Core->y), RNG.uniform(0, Core->z) };
                //Bank->addedNeutrons[idx].pos = Bank->addedNeutrons[idx].pos - bump;
                //Bank->addedNeutrons[idx].pos = randPos;
                //Bank->addedNeutrons[idx].dirVec = vec3::randomUnit(RNG);
                //seedNo[idx] = RNG.gen();
                //Bank->addedNeutrons[idx].Nullify();
                //atomicAdd(&(Bank->addedNeutronSize), -1);
                //return;
            }
        }
        //printf("idx %d neutron didn't reacted after %d loops..?\n", idx, iterTimes);
        seedNo[idx] = RNG.gen();
        return;
    }
    seedNo[idx] = RNG.gen();
    return;


}

G void normalizeCorrelatedSource(NeutronBank* Bank, CorrelatedSamplingConfig csConfig,
    double* csSourceSums, double* csSourceCount) {
    const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (!csConfig.enabled || idx >= Bank->allocatableNeutronNum) return;

    Neutron& child = Bank->addedNeutrons[idx];
    if (child.isNullified() || !child.csNeedsSourceNormalization) return;

    child.csSourceRatio[0] = 1.0;
    const double sourceCount = *csSourceCount;
    for (int level = 1; level <= csConfig.propagationGenerations; ++level) {
        const double sum = csSourceSums[level];
        child.csSourceRatio[level] = (sourceCount > 0.0 && sum > 0.0)
            ? child.csSourceRatio[level] * sourceCount / sum
            : 0.0;
    }
    for (int level = csConfig.propagationGenerations + 1;
        level <= CS_MAX_PROPAGATION_GENERATIONS; ++level) {
        child.csSourceRatio[level] = 1.0;
    }
    child.csNeedsSourceNormalization = false;
}
