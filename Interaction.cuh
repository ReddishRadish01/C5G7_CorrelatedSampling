#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "XSManager.cuh"
#include "CoreManager.cuh"

//#define INTERACTIONDEBUG 
//#define FISSIONPRINT



class Interaction {
public:
	HD static void scatter(Neutron& neutron, double outEnergy, GnuAMCM& RNG) {
		//neutron.updateWithLength(1.0e-10);
		neutron.energy = outEnergy;
		neutron.dirVec = vec3::randomUnit(RNG);
		//printf(" I SCATTER\n");
	}

	HD static void absorption(Neutron& n) {
		n.Nullify();
	}

	D static void fission(Neutron& n, NeutronBank* Bank, MatXS& matXS, GnuAMCM& RNG, double* k_mult, bool passFlag) {
		double nu = matXS.nu[static_cast<int>(n.energy) - 1];
		int fissionNum = static_cast<int>(nu / *k_mult + RNG.uniform(0.0, 1.0));

		n.dirVec = vec3::randomUnit(RNG);
		int fissionE = XSManager::returnFissionNeutronEnergy(matXS, RNG);
		int addIndex = atomicAdd(&(Bank->addedNeutronIndex), fissionNum - 1);
		int hardCap = addIndex + (fissionNum - 1);
		if (hardCap >= Bank->allocatableNeutronNum) {
			return;
		}

		atomicAdd(&(Bank->addedNeutronSize), fissionNum - 1);
		n.energy = fissionE;

		for (int i = 0; i < fissionNum - 1; i++) {
			//Bank->addedNeutrons[addIndex + i].status = true;
			int fissionE = XSManager::returnFissionNeutronEnergy(matXS, RNG);
			Bank->addedNeutrons[addIndex + i].reInitialize(n.pos, vec3::randomUnit(RNG), fissionE, 1.0, passFlag);
		}

	}

	D static void fission_forced(
		Neutron& n,
		NeutronBank* Bank,
		MatXS& matXS,
		GnuAMCM& RNG,
		double* k_mult,
		bool passFlag,
		double* fissionStrength,
		double& fissionSourceScore,
		double& fissionRateScore,
		MatXS& perturbedMatXS,
		const CorrelatedSamplingConfig& csConfig,
		double* csFissionScores,
		double* csSourceSums,
		double* csSourceCount,
		unsigned long long* csSupportMismatch
	) {
		
		int g = static_cast<int>(n.energy);
		double nu = matXS.nu[g - 1];
		if (nu == 0.0) { return; }

		//double motherWeight = n.weight;
		
		double sigC = collisionXS(matXS, g);	// this is transportXS: capture + fission + scattrer

		//
		double daughterWeight = n.weight * matXS.fisXS[g-1] / sigC;
		fissionRateScore = daughterWeight;
		fissionSourceScore = nu * daughterWeight;

		double csFissionRatio = 1.0;
		bool csValid = true;
		if (csConfig.enabled) {
			const double referenceNuSigmaF = nu * matXS.fisXS[g - 1];
			const double perturbedNuSigmaF = perturbedMatXS.nu[g - 1] * perturbedMatXS.fisXS[g - 1];
			if (n.csLogLikelihood < -700.0) {
				csFissionRatio = 0.0;
			}
			else if (!isfinite(n.csLogLikelihood) || n.csLogLikelihood > 700.0) {
				atomicAdd(csSupportMismatch, 1ULL);
				csValid = false;
			}
			else if (referenceNuSigmaF <= 0.0) {
				if (perturbedNuSigmaF > 0.0) atomicAdd(csSupportMismatch, 1ULL);
				csValid = false;
			}
			else {
				csFissionRatio = exp(n.csLogLikelihood) * perturbedNuSigmaF / referenceNuSigmaF;
				if (!isfinite(csFissionRatio)) {
					atomicAdd(csSupportMismatch, 1ULL);
					csValid = false;
				}
			}

			if (csValid) {
				for (int level = 0; level <= csConfig.propagationGenerations; ++level) {
					atomicAdd(&csFissionScores[level], fissionSourceScore * csFissionRatio * n.csSourceRatio[level]);
				}
			}
		}

		int fissionNum = static_cast<int>(nu / *k_mult * daughterWeight + RNG.uniform(0.0, 1.0));
		// for Fission source collision estimator, you have to sum the source weight, not the Bernoulli result of source weight.
		atomicAdd(fissionStrength, fissionSourceScore);
		//printf("fissionNum: %d\n", fissionNum);
		//n.dirVec = vec3::randomUnit(RNG);
		int addIndex = 0;

		if (fissionNum <= 0) {
			//printf("No fission neutron added: \n");
			return;
		}
		else {
			//printf("New neutron added from forced fision: %d. Parent Neutron weight: %f \n", fissionNum, motherWeight);
			//if (fissionNum >= 2) { printf("fissionNum greater than 2: %d\n", fissionNum); }
			addIndex = atomicAdd(&(Bank->addedNeutronIndex), fissionNum);
			if (addIndex + fissionNum > Bank->allocatableNeutronNum) {
				//printf("addIndex overflow: %d\n", addIndex);
				return;
			}
			atomicAdd(&(Bank->addedNeutronSize), fissionNum);
			if (csConfig.enabled) {
				atomicAdd(csSourceCount, static_cast<double>(fissionNum));
				for (int level = 1; level <= csConfig.propagationGenerations; ++level) {
					const double sourceRatio = csValid ? csFissionRatio * n.csSourceRatio[level - 1] : 0.0;
					atomicAdd(&csSourceSums[level], static_cast<double>(fissionNum) * sourceRatio);
				}
			}
			for (int i = 0; i < fissionNum; i++) {
				int fissionE = XSManager::returnFissionNeutronEnergy(matXS, RNG);
				Neutron& child = Bank->addedNeutrons[addIndex + i];
				child.reInitialize(n.pos, vec3::randomUnit(RNG), fissionE, 1.0, passFlag);
				if (csConfig.enabled) {
					for (int level = 1; level <= csConfig.propagationGenerations; ++level) {
						child.csSourceRatio[level] = csValid ? csFissionRatio * n.csSourceRatio[level - 1] : 0.0;
					}
					child.csNeedsSourceNormalization = true;
				}
			}
		}

	}

	D static void russianRoulette(Neutron& n, NeutronBank* Bank, GnuAMCM& RNG, double weightCutoff, double newWeight, bool addTrue) {
		if (n.weight < weightCutoff) {
			//double oldWeight = n.weight;
			double rouletteRNG = RNG.uniform(0.0, 1.0);
			double survivalProbability = n.weight / newWeight;
			if (rouletteRNG < survivalProbability) {
				//printf("Neutron survived RR with weight of %f\n", oldWeight);
				n.setWeight(newWeight);
			}
			else {
				n.Nullify();
				//printf("Neutron Nullified with weight of %f\n", oldWeight);
				if (addTrue == true) {
					atomicAdd(&(Bank->addedNeutronSize), -1);
				}
				else {
					atomicAdd(&(Bank->neutronSize), -1);
				}

			}
		}
	}

	D static void reaction_implicit(
		Neutron& n,
		NeutronBank* Bank,
		XSLibrary* XSLib,
		Pincell currentPincell,
		vec3 localPos,
		GnuAMCM& RNG,
		double* k_mult,
		bool passFlag,
		bool add,
		double* fissionWeight,
		double& fissionSourceScore,
		double& fissionRateScore,
		MatType perturbedMaterial,
		const CorrelatedSamplingConfig& csConfig,
		double* csFissionScores,
		double* csSourceSums,
		double* csSourceCount,
		unsigned long long* csSupportMismatch
	) {
		MatType currentMat = currentPincell.meatOrMod(localPos);
		double outEnergy = XSManager::scatteringEnergy(XSLib, currentMat, RNG, n.energy);
		MatXS& referenceXS = XSLib->returnMatXSByType(currentMat);
		MatXS& perturbedXS = XSLib->returnMatXSByType(perturbedMaterial);
		Interaction::fission_forced(
			n,
			Bank,
			referenceXS,
			RNG,
			k_mult,
			passFlag,
			fissionWeight,
			fissionSourceScore,
			fissionRateScore,
			perturbedXS,
			csConfig,
			csFissionScores,
			csSourceSums,
			csSourceCount,
			csSupportMismatch
		);

		MatXS& xs = referenceXS;
		int g = static_cast<int>(n.energy);
		double sigS = scatterSum(xs, g);
		//double sigT = xs.capXS[g - 1] + xs.fisXS[g - 1] + sigS; // isn't this just xs.transXS[g-1]?  yes, but we want to be explicit here.
		double sigT = xs.transXS[g - 1];
		n.setWeight(n.weight * sigS / sigT);

		if (csConfig.enabled) {
			const int outGroup = static_cast<int>(outEnergy);
			const double referenceScatter = referenceXS.elsXS[g - 1][outGroup - 1];
			const double perturbedScatter = perturbedXS.elsXS[g - 1][outGroup - 1];
			if (referenceScatter <= 0.0) {
				if (perturbedScatter > 0.0) atomicAdd(csSupportMismatch, 1ULL);
				n.csLogLikelihood = -1.0e300;
			}
			else if (perturbedScatter <= 0.0) {
				n.csLogLikelihood = -1.0e300;
			}
			else if (isfinite(n.csLogLikelihood)) {
				n.csLogLikelihood += log(perturbedScatter / referenceScatter);
			}
		}

		// implicit capture: get only the elastic XS out of all the XS sum
		//double elasXS = XSLib->returnMatXSByType(currentMat).transXS[static_cast<int>(n.energy)-1] - XSLib->returnMatXSByType(currentMat).capXS[static_cast<int>(n.energy) - 1] - XSLib->returnMatXSByType(currentMat).fisXS[static_cast<int>(n.energy) - 1];
		//double elasXS = XSLib->returnMatXSByType(currentMat).transXS[static_cast<int>(n.energy) - 1] - XSLib->returnMatXSByType(currentMat).absXS[static_cast<int>(n.energy) - 1];
		//n.setWeight(n.weight * elasXS / XSLib->returnMatXSByType(currentMat).transXS[static_cast<int>(n.energy) - 1]);
		//InteractionType interactionT = XSManager::returnInteracitonType(XSLib, currentMat, RNG, n.energy, outEnergy);
		//printf("currentEnergy: %f, outEnergy: %f\n", n.energy, outEnergy);
		Interaction::scatter(n, outEnergy, RNG);
	}

	D static void reaction(Neutron& n, NeutronBank* Bank, XSLibrary* XSLib, Pincell currentPincell, vec3 localPos, GnuAMCM& RNG, double* k_mult, bool passFlag, bool add) {
		MatType currentMat = currentPincell.meatOrMod(localPos);
		double outEnergy = 0.0;
		Interaction::fission(n, Bank, XSLib->returnMatXSByType(currentMat), RNG, k_mult, passFlag);

		InteractionType interactionT = XSManager::returnInteracitonType(XSLib, currentMat, RNG, n.energy, outEnergy);
		if (interactionT == InteractionType::nel) {
			
			Interaction::scatter(n, outEnergy, RNG);
		}
		else if (interactionT == InteractionType::ng) {

			Interaction::absorption(n);
			if (add == true) { // absorption in addedneutron
				atomicAdd(&(Bank->addedNeutronSize), -1);
			}
			else {
				atomicAdd(&(Bank->neutronSize), -1);
			}
		}
		else if (interactionT == InteractionType::nf) {
			Interaction::absorption(n);
			if (add == true) { // absorption in addedneutron
				atomicAdd(&(Bank->addedNeutronSize), -1);
			}
			else {
				atomicAdd(&(Bank->neutronSize), -1);
			}
		}

	}

	D static void reflection(Neutron& n, C5G7Geometry* Core, double DTS, double DTC, vec3 updatedSurfacePos, double eps) {

		//vec3 afterDTCPos = n.pos + n.dirVec * DTC;
		eps *= 10000;
		vec3 reflectVec = n.dirVec;

		if (updatedSurfacePos.x <= Core->xMin + eps &&
			Core->xMinBoundary == BoundaryCondition::Reflective &&
			n.dirVec.x < 0.0) {
			reflectVec.x = -reflectVec.x;
		}
		if (updatedSurfacePos.y <= Core->yMin + eps &&
			Core->yMinBoundary == BoundaryCondition::Reflective &&
			n.dirVec.y < 0.0) {
			reflectVec.y = -reflectVec.y;
		}
		if (updatedSurfacePos.z <= Core->zMin + eps &&
			Core->zMinBoundary == BoundaryCondition::Reflective &&
			n.dirVec.z < 0.0) {
			reflectVec.z = -reflectVec.z;
		}
		if (updatedSurfacePos.x >= Core->x - eps &&
			Core->xMaxBoundary == BoundaryCondition::Reflective &&
			n.dirVec.x > 0.0) {
			reflectVec.x = -reflectVec.x;
		}
		if (updatedSurfacePos.y >= Core->y - eps &&
			Core->yMaxBoundary == BoundaryCondition::Reflective &&
			n.dirVec.y > 0.0) {
			reflectVec.y = -reflectVec.y;
		}
		if (updatedSurfacePos.z >= Core->z - eps &&
			Core->zMaxBoundary == BoundaryCondition::Reflective &&
			n.dirVec.z > 0.0) {
			reflectVec.z = -reflectVec.z;
		}

		vec3 collisionPos = n.pos + n.dirVec * DTS * (1.0 - eps * 10);
		//vec3 collisionPos = n.pos + reflectVec * eps * 10000;
		//n.printInfo_Kernel(0);

		n.pos = collisionPos + reflectVec * eps;
		n.dirVec = reflectVec;
		//n.printInfo_Kernel(0);
	}
};
