#include "main.h"
#include <thread>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>
#include <inttypes.h>
#include <iterator>
#include <utility>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <set>
#include <vector>
#include <map>
#include <unordered_set>
#include <mutex>
using namespace std;

__constant__ ulong cols;		  // 8 bytes
__constant__ ulong rows;		  // 8 bytes
__device__ unsigned long long int numResultKendalls = 0;
__device__ unsigned long long int numResultSpearmans = 0;
__device__ unsigned long long int numResultNMI = 0;

__global__ void kendallTwoGenes(double *aResultKendalls, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor, float correctionThreshold)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
	ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
	if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{		
		long r1 = 0;
		long r2 = -1;
		long auxPat = pattern - rows + 1;
		if (auxPat < 0)
		{
			r2 = auxPat + rows;
		}
		for (ulong j = rows - 2; r2 == -1; j--)
		{
			auxPat = auxPat - j;
			r1++;
			if (auxPat < 0)
			{
				r2 = (j + auxPat) + (r1 + 1);
			}
		}

		if (r1 < rows && r2 < rows)
		{			
			int iConcordant = 0, iDiscordant = 0;
			double dKendall = -1;

			for (int iCol1 = 0; iCol1 < cols; iCol1++)
			{
				float fValueR1 = *(mDataGPU + r1 * cols + iCol1);
				float fValueR2 = *(mDataGPU + r2 * cols + iCol1);
				
				for (int iCol2 = 0; iCol2 < cols; iCol2++)
				{
					if (*(mDataGPU + r1 * cols + iCol2) > fValueR1)
					{
						if (*(mDataGPU + r2 * cols + iCol2) > fValueR2)
						{
							iConcordant += 1;
						}
						else
						{
							iDiscordant += 1;
						}
					}
				}
			}

			if (iConcordant + iDiscordant != 0)
			{ // Control division by zero
				dKendall = (double)(iConcordant - iDiscordant) / (iConcordant + iDiscordant);
				if (dKendall < 0)
				{ // Absolute value
					dKendall = dKendall * -1;
				}
				if (dKendall <= correctionThreshold)
				{
					dKendall = -1;
				}
			}

			*(aResultKendalls + idTh) = dKendall;
		}
	}
}

__global__ void spearmanCalcfDiGenesOne(double *fDiSpearman, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
		long r1 = 0;
		long r2 = -1;
		long auxPat = pattern - rows + 1;
		if (auxPat < 0)
		{
			r2 = auxPat + rows;
		}
		for (ulong j = rows - 2; r2 == -1; j--)
		{
			auxPat = auxPat - j;
			r1++;
			if (auxPat < 0)
			{
				r2 = (j + auxPat) + (r1 + 1);
			}
		}

		if (r1 < rows && r2 < rows)
		{
			for (int iConditions = 0; iConditions < cols; iConditions++)
			{
				int numEqualOrdX = 0;
				for (int iCont = 0; iCont < cols; iCont++)
				{
					if (iCont != iConditions)
					{
						if (*(mDataGPU + r1 * cols + iCont) < *(mDataGPU + r1 * cols + iConditions))
						{
							*(fDiSpearman + idTh * cols + iConditions) += 1;
						}
						else if (*(mDataGPU + r1 * cols + iCont) == *(mDataGPU + r1 * cols + iConditions))
						{
							numEqualOrdX += 1;
						}
					}
				}

				if (numEqualOrdX == 0)
				{
					*(fDiSpearman + idTh * cols + iConditions) = *(fDiSpearman + idTh * cols + iConditions) + 1;
				}
				else
				{
					*(fDiSpearman + idTh * cols + iConditions) = *(fDiSpearman + idTh * cols + iConditions) + 1 + (numEqualOrdX / 2.0);
				}
			}
		}
	}
}

__global__ void spearmanCalcfDiGenesTwo(double *fDiSpearman, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
		long r1 = 0;
		long r2 = -1;
		long auxPat = pattern - rows + 1;
		if (auxPat < 0)
		{
			r2 = auxPat + rows;
		}
		for (ulong j = rows - 2; r2 == -1; j--)
		{
			auxPat = auxPat - j;
			r1++;
			if (auxPat < 0)
			{
				r2 = (j + auxPat) + (r1 + 1);
			}
		}

		if (r1 < rows && r2 < rows)
		{
			for (int iConditions = 0; iConditions < cols; iConditions++)
			{
				int numEqualOrdY = 0;
				float fValueR2 = *(mDataGPU + r2 * cols + iConditions);
				double dFi = *(fDiSpearman + idTh * cols + iConditions);
				for (int iCont = 0; iCont < cols; iCont++)
				{
					if (iCont != iConditions)
					{
						float fValueCont = *(mDataGPU + r2 * cols + iCont);
						if (fValueCont < fValueR2)
						{
							dFi -= 1;
						}
						else if (fValueCont == fValueR2)
						{
							numEqualOrdY += 1;
						}
					}
				}
				
				if (numEqualOrdY == 0)
				{
					dFi = dFi - 1;
				}
				else
				{
					dFi = dFi - 1 - (numEqualOrdY / 2.0);
				}

				*(fDiSpearman + idTh * cols + iConditions) = dFi;
			}
		}
	}
}

__global__ void spearmanCalc(double *aResultSpearmans, double *fDiSpearman, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor, float correctionThreshold)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
			double diSquare = 0;
			for (int iConditions = 0; iConditions < cols; iConditions++)
			{
				diSquare = diSquare + (*(fDiSpearman + idTh * cols + iConditions) * *(fDiSpearman + idTh * cols + iConditions));
			}
			double dSpearman = 1 - ((6 * diSquare) / (cols * ((cols * cols) - 1)));
			if(dSpearman < 0){
				dSpearman = dSpearman * -1;
			}
			if(dSpearman <= correctionThreshold){
				dSpearman = -1;
			}
			*(aResultSpearmans + idTh) = dSpearman;
	}
}

__global__ void nmiCalcMutualInformation(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, int *mDataNormalizedGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
		long r1 = 0;
		long r2 = -1;
		long auxPat = pattern - rows + 1;
		if (auxPat < 0)
		{
			r2 = auxPat + rows;
		}
		for (ulong j = rows - 2; r2 == -1; j--)
		{
			auxPat = auxPat - j;
			r1++;
			if (auxPat < 0)
			{
				r2 = (j + auxPat) + (r1 + 1);
			}
		}

		if (r1 < rows && r2 < rows)
		{
			int maxVal = *(mDataNormalizedGPU + r1 * (cols + 1) + cols);

			// Clean dNMIResults by GPU device
			for (int iColumn = 0; iColumn < 11; ++iColumn) {
				*(dNMIResults + idTh * 11 + iColumn) = 0;
			}	
			
			// dNMIResults
			// [0] --> Mutual information and NMI // [1] -->  entropyGen1 // [2] --> entropyGen2 // [3 - 10] --> probMaps (calcMI)
			for (int iColumn = 0; iColumn < cols; ++iColumn) {
				int valGen1Column = *(mDataNormalizedGPU + r1 * (cols + 1) + iColumn);
				int valGen2Column = *(mDataNormalizedGPU + r2 * (cols + 1) + iColumn);	

				*(dNMIResults + idTh * 11 + (valGen1Column + 3)) = *(dNMIResults + idTh * 11 + (valGen1Column + 3)) + 1;
				*(dNMIResults + idTh * 11 + (valGen2Column + 5)) = *(dNMIResults + idTh * 11 + (valGen2Column + 5)) + 1;
				*(dNMIResults + idTh * 11 + ((valGen1Column + maxVal * valGen2Column) + 7)) = *(dNMIResults + idTh * 11 + ((valGen1Column + maxVal * valGen2Column) + 7)) + 1;
			}

			for (int iCont = 0; iCont < 8; iCont++) {
				*(dNMIResults + idTh * 11 + (3 + iCont)) = *(dNMIResults + idTh * 11 + (3 + iCont)) / cols;
			}			

			double mi = 0;
			for (int iCont = 0; iCont < 4; iCont++) {
				if (*(dNMIResults + idTh * 11 + (7 + iCont)) > 0 && *(dNMIResults + idTh * 11 + ((iCont%maxVal)+3)) > 0 && *(dNMIResults + idTh * 11 + (iCont/maxVal)+5) > 0) {
					mi += *(dNMIResults + idTh * 11 + (7 + iCont)) * logf(*(dNMIResults + idTh * 11 + (7 + iCont)) / *(dNMIResults + idTh * 11 + ((iCont%maxVal)+3)) / *(dNMIResults + idTh * 11 + (iCont/maxVal)+5));
				}
			}		

			mi = mi / logf(2);			
			*(dNMIResults + idTh * 11 + 0) = mi;
		}
	}
}

__global__ void nmiCalcEntropy(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, int *mDataNormalizedGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
		long r1 = 0;
		long r2 = -1;
		long auxPat = pattern - rows + 1;
		if (auxPat < 0)
		{
			r2 = auxPat + rows;
		}
		for (ulong j = rows - 2; r2 == -1; j--)
		{
			auxPat = auxPat - j;
			r1++;
			if (auxPat < 0)
			{
				r2 = (j + auxPat) + (r1 + 1);
			}
		}

		if (r1 < rows && r2 < rows)
		{
			// Clean auxiliar data
			*(dNMIResults + idTh * 11 + 3) = 0;
			*(dNMIResults + idTh * 11 + 4) = 0;
			*(dNMIResults + idTh * 11 + 5) = 0;
			*(dNMIResults + idTh * 11 + 6) = 0;
			*(dNMIResults + idTh * 11 + 7) = 0;
			*(dNMIResults + idTh * 11 + 8) = 0;
			*(dNMIResults + idTh * 11 + 9) = 0;
			*(dNMIResults + idTh * 11 + 10) = 0;

			// dNMIResults
			// [0] --> Mutual information and NMI // [1] -->  entropyGen1 // [2] --> entropyGen2 // [3 - 10] --> probMaps (calcMI)
			for (int iColumn = 0; iColumn < cols; ++iColumn) {
				int valGen1Column = *(mDataNormalizedGPU + r1 * (cols + 1) + iColumn);
				int valGen2Column = *(mDataNormalizedGPU + r2 * (cols + 1) + iColumn);

				*(dNMIResults + idTh * 11 + (valGen1Column + 3)) = *(dNMIResults + idTh * 11 + (valGen1Column + 3)) + 1;
				*(dNMIResults + idTh * 11 + (valGen2Column + 5)) = *(dNMIResults + idTh * 11 + (valGen2Column + 5)) + 1;
			}

			for (int iCont = 0; iCont < 4; iCont++) {
				*(dNMIResults + idTh * 11 + (3 + iCont)) = *(dNMIResults + idTh * 11 + (3 + iCont)) / cols;
			}			

			float dEntropyG1, dEntropyG2 = 0;
			for (int iCont = 0; iCont < 2; iCont++) {
				float varAuxG1 = *(dNMIResults + idTh * 11 + (3 + iCont));
				float varAuxG2 = *(dNMIResults + idTh * 11 + (5 + iCont));
				if(varAuxG1 > 0){
					dEntropyG1 = dEntropyG1 - varAuxG1 * logf(varAuxG1);
				}

				if(varAuxG2 > 0){
					dEntropyG2 = dEntropyG2 - varAuxG2 * logf(varAuxG2);
				}
			}		

			dEntropyG1 = dEntropyG1 / logf(2);		
			dEntropyG2 = dEntropyG2 / logf(2);	
			*(dNMIResults + idTh * 11 + 1) = dEntropyG1;
			*(dNMIResults + idTh * 11 + 2) = dEntropyG2;
		}
	}
}

__global__ void nmiCalcValue(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor, float correctionThreshold){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs && pattern < (pairsPerGpuPrevious+totalFor))
	{
		float dNMI = -1;
		float denom = (*(dNMIResults + idTh * 11 + 1) + *(dNMIResults + idTh * 11 + 2));
		if(denom != 0){
			dNMI = 2.0 * *(dNMIResults + idTh * 11 + 0) / (*(dNMIResults + idTh * 11 + 1) + *(dNMIResults + idTh * 11 + 2));
		}

		if(dNMI > 300){
			dNMI = -1;
		}
		
		if(dNMI < correctionThreshold){
			dNMI = -1;
		}

		*(dNMIResults + idTh * 11 + 0) = dNMI;
	}
}

// #############
// # GPU FUNCS #
// #############
void getNumPairs()
{
	maxPairs = 0;
	for (int i = 0; i < ulRowsData; i++)
	{
		for (int j = i + 1; j < ulRowsData; j++)
		{
			maxPairs++;
		}
	}
}

void prepareGpu1D(ulong lNumber)
{
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, device);
	lastBlocksGrid = 1;
	maxIteratorGPU = 0;
	maxThreadsPerBlock = lNumber; // Case 1: 0 < lNumber <= prop.maxThreadsPerBlock
	if (lNumber > prop.maxThreadsPerBlock)
	{ // Case 2: lNumber > prop.maxThreadsPerBlock && Supported GPU in a for
		maxThreadsPerBlock = prop.maxThreadsPerBlock;
		maxBlocksPerGrid = lNumber / prop.maxThreadsPerBlock;
		lastBlocksGrid = lNumber / prop.maxThreadsPerBlock;
		if (lNumber % prop.maxThreadsPerBlock != 0)
		{
			maxBlocksPerGrid++;
			lastBlocksGrid++;
		}
		
		if (maxBlocksPerGrid > prop.maxGridSize[1])
		{ // Case 3: Not supported GPU with a for --> Split patterns in multiple for
			maxIteratorGPU = maxBlocksPerGrid / prop.maxGridSize[1];
			lastBlocksGrid = maxBlocksPerGrid - (maxIteratorGPU * prop.maxGridSize[1]);
			maxBlocksPerGrid = prop.maxGridSize[1];
		}
	}
}

long long *getPairsFiltered(int id, ulong numPatterns)
{
	long long *patFiltered;
	cudaMalloc((void **)&patFiltered, numPatterns * sizeof(long long));
	cudaMemset(patFiltered, -1, numPatterns * sizeof(long long));
	return patFiltered;
}

void threadsPerDevice(int id, cudaStream_t s, ulong chunks,
					  ulong pairsPerGpuPrevious, ulong pairsPerRun, float *mDataGPU, int *mDataNormalizedGPU,
					  mutex *m)
{

	cudaSetDevice(id);
	ulong totalPairs = 0;
	ofstream fileResults;
	
	if (bOutput)
	{
		fileResults.open("results/results_GPU_" + to_string(id) + ".csv");
	}
	
	for (ulong largeScale = 0; largeScale < chunks; largeScale++)
	{
		cudaSetDevice(id);
		// 1) Validation
		double *aResultKendalls, *aResultSpearmans, *fDiSpearman;
		float *dNMIResults;
		unsigned long long int numResultKendallsCpu = 0;
		unsigned long long int numResultSpearmanCpu = 0;
		unsigned long long int numResultNMICpu = 0;
		cudaMalloc((void **)&aResultKendalls, pairsPerRun * sizeof(double));
		cudaMemset(aResultKendalls, 0, pairsPerRun * sizeof(double));
		cudaMalloc((void **)&aResultSpearmans, pairsPerRun * sizeof(double));
		cudaMemset(aResultSpearmans, 0, pairsPerRun * sizeof(double));
		cudaMalloc((void **)&fDiSpearman, pairsPerRun * ulColsData * sizeof(double));
		cudaMemset(fDiSpearman, 0, pairsPerRun * ulColsData * sizeof(double));
		cudaMalloc((void **)&dNMIResults, pairsPerRun * 11 * sizeof(float));
		cudaMemset(dNMIResults, 0, pairsPerRun * 11 * sizeof(float));

		// Kendall
		prepareGpu1D(pairsPerRun);
		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			kendallTwoGenes<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(aResultKendalls, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdKendall);
		}
		kendallTwoGenes<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(aResultKendalls, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdKendall);

		// Spearman
		prepareGpu1D(pairsPerRun);
		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCalcfDiGenesOne<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		spearmanCalcfDiGenesOne<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCalcfDiGenesTwo<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		spearmanCalcfDiGenesTwo<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCalc<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(aResultSpearmans, fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdSpearman);
		}
		spearmanCalc<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(aResultSpearmans, fDiSpearman, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdSpearman);
		
		// NMI
		prepareGpu1D(pairsPerRun);
		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			nmiCalcMutualInformation<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, mDataNormalizedGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		nmiCalcMutualInformation<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, mDataNormalizedGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			nmiCalcEntropy<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, mDataNormalizedGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		nmiCalcEntropy<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, mDataNormalizedGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			nmiCalcValue<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdNMI);
		}
		nmiCalcValue<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid, correctionThresholdNMI);

		// 2) Save results
		if (bOutput)
		{
			// Transfer global memory to RAM
			double *aResultKendallsCpu = (double *)malloc(pairsPerRun * sizeof(double));
			cudaMemcpy(aResultKendallsCpu, aResultKendalls, pairsPerRun * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpyFromSymbol(&numResultKendallsCpu, numResultKendalls,
								 sizeof(unsigned long long int), 0, cudaMemcpyDeviceToHost);

			double *aResultSpearmansCpu = (double *)malloc(pairsPerRun * sizeof(double));
			cudaMemcpy(aResultSpearmansCpu, aResultSpearmans, pairsPerRun * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpyFromSymbol(&numResultKendallsCpu, numResultKendalls,
								 sizeof(unsigned long long int), 0, cudaMemcpyDeviceToHost);

			float *dNMIResultsCpu = (float *)malloc(pairsPerRun * 11 * sizeof(float));
			cudaMemcpy(dNMIResultsCpu, dNMIResults, pairsPerRun * 11 * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpyFromSymbol(&numResultNMICpu, numResultNMI,
								 sizeof(unsigned long long int), 0, cudaMemcpyDeviceToHost);

			fileResults << "g1;g2;kendall_tau;spearman;nmi"
						<< "\n";
			for (ulong iRow = 0; iRow < pairsPerRun; iRow++)
			{
				ulong pattern = iRow + pairsPerGpuPrevious;
				if (pattern < maxPairs && pattern < (pairsPerGpuPrevious + (maxThreadsPerBlock*maxBlocksPerGrid)))
				{
					long r1 = 0;
					long r2 = -1;
					long auxPat = pattern - ulRowsData + 1;
					if (auxPat < 0)
					{
						r2 = auxPat + ulRowsData;
					}
					for (ulong j = ulRowsData - 2; r2 == -1; j--)
					{
						auxPat = auxPat - j;
						r1++;
						if (auxPat < 0)
						{
							r2 = (j + auxPat) + (r1 + 1);
						}
					}

					int iMajorVoting = 0;
					double dKendall = *(aResultKendallsCpu + iRow);
					double dSpearman = *(aResultSpearmansCpu + iRow);
					float dNMI = *(dNMIResultsCpu + iRow * 11 + 0);
					if(dKendall != -1){
						iMajorVoting++;
					}
	
					if(dSpearman != -1){
						iMajorVoting++;
					}

					if(dNMI != -1){
						iMajorVoting++;
					}
					
					if (iMajorVoting >= 2)
					{
						fileResults << r1 << ";" << r2 << ";" << dKendall << ";" << dSpearman << ";" << dNMI << "\n";
					}
				}
			}

			free(aResultKendallsCpu);
			free(aResultSpearmansCpu);
			free(dNMIResultsCpu);
		}

		cudaFree(aResultKendalls);
		cudaFree(aResultSpearmans);
		cudaFree(fDiSpearman);
		cudaFree(dNMIResults);
		numResultKendallsCpu = 0;
		numResultSpearmanCpu = 0;
		numResultNMICpu = 0;
		cudaMemcpyToSymbol(numResultKendalls, &numResultKendallsCpu, sizeof(unsigned long long int), 0, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(numResultSpearmans, &numResultSpearmanCpu, sizeof(unsigned long long int), 0, cudaMemcpyHostToDevice);	
		cudaMemcpyToSymbol(numResultNMI, &numResultNMICpu, sizeof(unsigned long long int), 0, cudaMemcpyHostToDevice);			
		totalPairs += pairsPerRun;
	}

	if (bOutput)
	{
		fileResults.close();
	}
}

double runAlgorithm()
{
	for (int i = 0; i < iDeviceCount; ++i)
	{
		cudaSetDevice(i);
		cudaMemcpyToSymbol(*(&cols), &ulColsData, sizeof(ulong), 0,
						   cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(*(&rows), &ulRowsData, sizeof(ulong), 0,
						   cudaMemcpyHostToDevice);
	}

	getNumPairs();

	// 2) PREPARING LARGE-SCALE DATA: CHUNKS
	cudaStream_t s[iDeviceCount];
	thread threads[iDeviceCount];
	ulong chunks[iDeviceCount], pairsPerRun[iDeviceCount];
	ulong pairsPerGpu = maxPairs / iDeviceCount;
	ulong restPairsPerGpu = maxPairs % iDeviceCount;

	for (int i = 0; i < iDeviceCount; ++i)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&s[i]);
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		double availableMemory = ((3 * prop.totalGlobalMem) / 4 - (ulRowsData * ulColsData * sizeof(float)) - (ulRowsData * (ulColsData+1) * sizeof(int))); // mData + mDataNormalized
		double sizeResultKendall = 0, sizeResultsSpearman = 0, fDiSpearman = 0;
		float dNMIResults = 0;
		sizeResultKendall = (pairsPerGpu * sizeof(double));
		sizeResultsSpearman = (pairsPerGpu * sizeof(double));
		fDiSpearman = (pairsPerGpu * ulColsData * sizeof(double));
		dNMIResults = (pairsPerGpu * 11 * sizeof(float)); // [0] --> Mutual information and NMI [1] --> entropyGen1 [2] --> entropyGen2 [3 - 10] --> probMaps (calcMI) [3 - 4] --> calcEnropy
		chunks[i] = ((sizeResultKendall + sizeResultsSpearman + fDiSpearman + dNMIResults) / availableMemory) + 1;
		pairsPerRun[i] = pairsPerGpu / chunks[i];

		if (pairsPerGpu % chunks[i] != 0)
		{
			pairsPerRun[i]++;
		}
		if (iDeviceCount > 1 && maxPairs % iDeviceCount != 0 && i == iDeviceCount - 1)
		{
			pairsPerRun[i] += restPairsPerGpu;
		}
	}

	ulong pairsPerGpuPrevious = 0;

	mutex m;
	struct timeval stop, start;
	float *mDataGPU;
	int *mDataNormalizedGPU;
	for (int i = 0; i < iDeviceCount; ++i)
	{
		cudaSetDevice(i);
		cudaMallocHost((void **)&mDataGPU, ulRowsData * ulColsData * sizeof(float));
		cudaMemcpy(mDataGPU, mData, ulRowsData * ulColsData * sizeof(float),
				   cudaMemcpyHostToDevice);
		cudaMallocHost((void **)&mDataNormalizedGPU, ulRowsData * (ulColsData+1) * sizeof(int));
		cudaMemcpy(mDataNormalizedGPU, mDataNormalized, ulRowsData * (ulColsData+1) * sizeof(int),
							  cudaMemcpyHostToDevice);
	}
	for (int i=0; i < iDeviceCount; ++i)
	{
		if (i > 0)
		{
			pairsPerGpuPrevious += chunks[i - 1] * pairsPerRun[i - 1];
		}

		threads[i] = thread(threadsPerDevice, i, s[i], chunks[i], pairsPerGpuPrevious, pairsPerRun[i], mDataGPU, mDataNormalizedGPU, &m);
	}

	gettimeofday(&start, NULL);
	for (auto &th : threads)
	{
		th.join();
	}
	gettimeofday(&stop, NULL);


	free(mData);
	for (int i = 0; i < iDeviceCount; ++i)
    {
        cudaSetDevice(i);
        cudaDeviceReset();
    }
	return (((stop.tv_sec - start.tv_sec) * 1000.0) + ((stop.tv_usec - start.tv_usec) / 1000.0)) / 1000.0;
}
