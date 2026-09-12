#pragma once
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


G void fetchNeutronNumber(NeutronBank* Bank, double* d_NeutronNum, double* d_NeutronWeightSum);
G void globalNeutronWeightOffset(NeutronBank* Bank, double* weightOffset);

G void cycle_Neutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes, double* fissionCount,
	CorrelatedSamplingConfig csConfig, double* csFissionScores, double* csSourceSums,
	double* csSourceCount, unsigned long long* csSupportMismatch);
G void addedNeutronPassResetter(NeutronBank* Neutrons);
G void cycle_addedNeutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes, double* fissionCount,
	CorrelatedSamplingConfig csConfig, double* csFissionScores, double* csSourceSums,
	double* csSourceCount, unsigned long long* csSupportMismatch);
G void normalizeCorrelatedSource(NeutronBank* Bank, CorrelatedSamplingConfig csConfig,
	double* csSourceSums, double* csSourceCount);
