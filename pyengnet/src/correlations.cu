#include "main.h"
using namespace std;

__constant__ ulong cols;		  // 8 bytes
__constant__ ulong rows;		  // 8 bytes
__constant__ int maxValueDatasetGPU;		  // 8 bytes
__device__ unsigned long long int numResultKendalls = 0;
__device__ unsigned long long int numResultSpearmans = 0;
__device__ unsigned long long int numResultNMI = 0;

// Kendall
__global__ void kendallTwoGenes(ulong lCombination, double *aResultKendalls, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
	ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
	if (pattern < maxPairs)
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
			int iConcordant = 0, iDiscordant = 0, tiersGene1 = 0, tiersGene2 = 0;
			double dKendall = -1;

			// 1) Calc maxValue index of Col1
			int iPosMaxValue = 0;
			float fMaxValue = *(mDataGPU + r1 * cols + 0);
			for (int iCol1 = 1; iCol1 < cols; iCol1++)
			{
				if(fMaxValue < *(mDataGPU + r1 * cols + iCol1)){
					fMaxValue = *(mDataGPU + r1 * cols + iCol1);
					iPosMaxValue = iCol1;
				}
			}

			// 2) Calc concordant and discordant pairs	
			for (int iCol1 = 0; iCol1 < cols; iCol1++)
			{
				if(iCol1 != iPosMaxValue)
				{
					for (int iCol2 = 0; iCol2 < cols; iCol2++)
					{		
																				
						if(iCol1 != iCol2 && *(mDataGPU + r1 * cols + iCol1) < *(mDataGPU + r1 * cols + iCol2))
						{								
							if (*(mDataGPU + r2 * cols + iCol2) > *(mDataGPU + r2 * cols + iCol1))
							{
								iConcordant += 1;
							}
							
							if(*(mDataGPU + r2 * cols + iCol2) < *(mDataGPU + r2 * cols + iCol1))
							{
								iDiscordant += 1;
							}
						}											
					}
				}				
			}

			// 3) Calc tiers
			for (int iCol1 = 0; iCol1 < cols; iCol1++)
			{
				for (int iCol2 = iCol1 + 1; iCol2 < cols; iCol2++)
				{
					// Gene 1
					if(*(mDataGPU + r1 * cols + iCol1) == *(mDataGPU + r1 * cols + iCol2)){
						tiersGene1 += 1;
					}

					if(*(mDataGPU + r2 * cols + iCol1) == *(mDataGPU + r2 * cols + iCol2)){
						tiersGene2 += 1;
					}
				}
			}

			dKendall = (double)(iConcordant - iDiscordant) / (double)(sqrtf((lCombination - tiersGene1) * (lCombination - tiersGene2)));

			*(aResultKendalls + idTh) = dKendall;
		}		
	}
}

// Spearman
/*__global__ void spearmanCalcRanks(double *fSpearmanStats, double *fSpearmanRankG1, double *fSpearmanRankG2, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs)
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
			double dMeanG1 = 0, dMeanG2 = 0;
			*(fSpearmanStats + idTh * 4 + 0) = 0;
			*(fSpearmanStats + idTh * 4 + 1) = 0;
			for (int iConditions = 0; iConditions < cols; iConditions++)
			{
				double numEqualOrdG1 = 1.0, readG1 = 1.0, numEqualOrdG2 = 1.0, readG2 = 1.0;
				for (int iCont = 0; iCont < cols; iCont++)
				{
					if (iCont != iConditions)
					{
						if (*(mDataGPU + r1 * cols + iCont) < *(mDataGPU + r1 * cols + iConditions))
						{
							readG1 += 1.0;
						}
						else if (*(mDataGPU + r1 * cols + iCont) == *(mDataGPU + r1 * cols + iConditions))
						{
							numEqualOrdG1 += 1.0;
						}

						if (*(mDataGPU + r2 * cols + iCont) < *(mDataGPU + r2 * cols + iConditions))
						{
							readG2 += 1.0;
						}
						else if (*(mDataGPU + r2 * cols + iCont) == *(mDataGPU + r2 * cols + iConditions))
						{
							numEqualOrdG2 += 1.0;
						}
					}
				}
				*(fSpearmanRankG1 + idTh * cols + iConditions) = readG1 + ((numEqualOrdG1-1.0) * 0.5);
				*(fSpearmanRankG2 + idTh * cols + iConditions) = readG2 + ((numEqualOrdG2-1.0) * 0.5);
				*(fSpearmanStats + idTh * 4 + 0) += *(fSpearmanRankG1 + idTh * cols + iConditions);
				*(fSpearmanStats + idTh * 4 + 1) += *(fSpearmanRankG2 + idTh * cols + iConditions);				
			}

			*(fSpearmanStats + idTh * 4 + 0) = *(fSpearmanStats + idTh * 4 + 0) / cols;
			*(fSpearmanStats + idTh * 4 + 1) = *(fSpearmanStats + idTh * 4 + 1) / cols;
		}
	}
}*/

__global__ void spearmanCalcRanks(double *fSpearmanStats, double *fSpearmanRankG1, double *fSpearmanRankG2, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs)
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
			*(fSpearmanStats + idTh * 4 + 0) = 0;
			*(fSpearmanStats + idTh * 4 + 1) = 0;
			for (int iConditions = 0; iConditions < cols; iConditions++)
			{
				double numEqualOrdG1 = 1.0, readG1 = 1.0, numEqualOrdG2 = 1.0, readG2 = 1.0;
				for (int iCont = 0; iCont < cols; iCont++)
				{
					if (iCont != iConditions)
					{
						if (*(mDataGPU + r1 * cols + iCont) < *(mDataGPU + r1 * cols + iConditions))
						{
							readG1 += 1.0;
						}
						else if (*(mDataGPU + r1 * cols + iCont) == *(mDataGPU + r1 * cols + iConditions))
						{
							numEqualOrdG1 += 1.0;
						}

						if (*(mDataGPU + r2 * cols + iCont) < *(mDataGPU + r2 * cols + iConditions))
						{
							readG2 += 1.0;
						}
						else if (*(mDataGPU + r2 * cols + iCont) == *(mDataGPU + r2 * cols + iConditions))
						{
							numEqualOrdG2 += 1.0;
						}
					}
				}
				*(fSpearmanRankG1 + idTh * cols + iConditions) = readG1 + ((numEqualOrdG1-1.0) * 0.5);
				*(fSpearmanRankG2 + idTh * cols + iConditions) = readG2 + ((numEqualOrdG2-1.0) * 0.5);
				*(fSpearmanStats + idTh * 4 + 0) += *(fSpearmanRankG1 + idTh * cols + iConditions);
				*(fSpearmanStats + idTh * 4 + 1) += *(fSpearmanRankG2 + idTh * cols + iConditions);				
			}

			*(fSpearmanStats + idTh * 4 + 0) = *(fSpearmanStats + idTh * 4 + 0) / cols;
			*(fSpearmanStats + idTh * 4 + 1) = *(fSpearmanStats + idTh * 4 + 1) / cols;
		}
	}
}

__global__ void spearmanCovariance(double *fSpearmanStats, double *aResultSpearmans, double *fSpearmanRankG1, double *fSpearmanRankG2, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
    ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
    if (pattern < maxPairs)
	{
		//double dCovariance = 0;
		*(fSpearmanStats + idTh * 4 + 2) = 0;
		for (int iConditions = 0; iConditions < cols; iConditions++)
		{
			*(fSpearmanStats + idTh * 4 + 2) += (*(fSpearmanRankG1 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 0)) * (*(fSpearmanRankG2 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 1));
		}

		*(aResultSpearmans + idTh) = *(fSpearmanStats + idTh * 4 + 2) / (cols - 1);
	}
}

__global__ void spearmanCalcValue(double *fSpearmanStats, double *aResultSpearmans, double *fSpearmanRankG1, double *fSpearmanRankG2, ulong maxPairs, int id, ulong pairsPerGpuPrevious, float *mDataGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor)
{
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
    ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
    if (pattern < maxPairs)
	{
		//double dStdG1 = 0, dStdG2 = 0;
		*(fSpearmanStats + idTh * 4 + 2) = 0;
		*(fSpearmanStats + idTh * 4 + 3) = 0;
		for (int iConditions = 0; iConditions < cols; iConditions++)
		{
			*(fSpearmanStats + idTh * 4 + 2) += (*(fSpearmanRankG1 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 0)) * (*(fSpearmanRankG1 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 0));
			*(fSpearmanStats + idTh * 4 + 3) += (*(fSpearmanRankG2 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 1)) * (*(fSpearmanRankG2 + idTh * cols + iConditions) - *(fSpearmanStats + idTh * 4 + 1));
		}

		*(fSpearmanStats + idTh * 4 + 2) = (double)(sqrtf(*(fSpearmanStats + idTh * 4 + 2) / (cols - 1)));
		*(fSpearmanStats + idTh * 4 + 3) = (double)(sqrtf(*(fSpearmanStats + idTh * 4 + 3) / (cols - 1)));

		*(aResultSpearmans + idTh) = *(aResultSpearmans + idTh) / (*(fSpearmanStats + idTh * 4 + 2) * *(fSpearmanStats + idTh * 4 + 3));
	}
}

// NMI
__global__ void nmiCalcMutualInformation(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, int *mDataNormalizedGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs)
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

			// Clean dNMIResults by GPU device
			for (int i = 3; i < ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3); i++){
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) = -1;
			}	
			
			// dNMIResults
			// [0] --> Mutual information and NMI // [1] -->  entropyGen1 // [2] --> entropyGen2 // [3] - ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3)] --> probMaps (calcMI)
			for (int iColumn = 0; iColumn < cols; ++iColumn) {
				int valGen1Column = *(mDataNormalizedGPU + r1 * (cols + 1) + iColumn);				
				int valGen2Column = *(mDataNormalizedGPU + r2 * (cols + 1) + iColumn);

				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (valGen1Column + 3)) = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (valGen1Column + 3)) + 1;
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (maxValueDatasetGPU + 3 + valGen2Column)) = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (maxValueDatasetGPU + 3 + valGen2Column)) + 1;
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((valGen1Column + maxValueDatasetGPU * valGen2Column) + (maxValueDatasetGPU + maxValueDatasetGPU) + 3)) = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((valGen1Column + maxValueDatasetGPU * valGen2Column) + (maxValueDatasetGPU + maxValueDatasetGPU) + 3)) + 1;
			}

			for (int i = 3; i < ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3); i++) {
				if(*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) != -1){
					*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) = (*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) + 1) / cols;
				}				
			}	

			double mi = 0;
			for (int iCont = 0; iCont < (maxValueDatasetGPU * maxValueDatasetGPU); iCont++) {
				if(*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((maxValueDatasetGPU + maxValueDatasetGPU) + 3 + iCont)) > 0){
					float doubleValue = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((maxValueDatasetGPU + maxValueDatasetGPU) + 3 + iCont));
					float doubleValue2 = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((iCont%maxValueDatasetGPU)+3));
					float doubleValue3 = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + ((iCont/maxValueDatasetGPU)+3+maxValueDatasetGPU));
					if(doubleValue > 0 && doubleValue2 > 0 && doubleValue3 > 0){
						mi += doubleValue * logf(doubleValue / doubleValue2 / doubleValue3);
					}					
				}
			}	
			
			mi = mi / logf(2);	
			*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 0) = mi;
		}
	}
}

__global__ void nmiCalcEntropy(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, int *mDataNormalizedGPU, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
    ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
    if (pattern < maxPairs)
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
			int maxValGene1 = *(mDataNormalizedGPU + r1 * (cols + 1) + cols);
			int maxValGene2 = *(mDataNormalizedGPU + r2 * (cols + 1) + cols);
			
			// Clean auxiliar data
			for (int i = 3; i < ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3); i++){
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) = -1;
			}

			// dNMIResults
			// [0] --> Mutual information and NMI // [1] -->  entropyGen1 // [2] --> entropyGen2 // [3 - ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3)] --> probMaps (calcMI)
			for (int iColumn = 0; iColumn < cols; ++iColumn) {
				int valGen1Column = *(mDataNormalizedGPU + r1 * (cols + 1) + iColumn);				
				int valGen2Column = *(mDataNormalizedGPU + r2 * (cols + 1) + iColumn);
				
				// valGen1Column
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (valGen1Column + 3)) = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (valGen1Column + 3)) + 1;
				
				// valGen2Column
				*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (maxValueDatasetGPU + 3 + valGen2Column)) = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (maxValueDatasetGPU + 3 + valGen2Column)) + 1;
			}

			for (int i = 3; i < ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3); i++) {
				if(*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) != -1){
					*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) = (*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + i) + 1) / cols;
				}				
			}	

			float dEntropyG1 = 0;
			for(int i = 0; i < maxValGene1; i++){
				float varAuxG1 = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (3 + i));
				if(varAuxG1 > 0){
					dEntropyG1 = dEntropyG1 - varAuxG1 * logf(varAuxG1);
				}
			}

			float dEntropyG2 = 0;
			for(int i = 0; i < maxValGene2; i++){				
				float varAuxG2 = *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + (3 + i + maxValueDatasetGPU));
				if(varAuxG2 > 0){
					dEntropyG2 = dEntropyG2 - varAuxG2 * logf(varAuxG2);
				}
			}

			dEntropyG1 = dEntropyG1 / logf(2);		
			dEntropyG2 = dEntropyG2 / logf(2);	
			*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 1) = dEntropyG1;
			*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 2) = dEntropyG2;
		}
	}
}

__global__ void nmiCalcValue(float *dNMIResults, ulong maxPairs, int id, ulong pairsPerGpuPrevious, ulong totalPairs, ulong pairsPerRun, int iter, ulong totalFor){
	ulong idTh = blockIdx.x * blockDim.x + threadIdx.x;
        ulong pattern = idTh + (totalFor * (iter - 1)) + pairsPerGpuPrevious + totalPairs;
        if (pattern < maxPairs)
		{
			float dNMI = 0;
			float denom = (*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 1) + *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 2));
			if(denom > 0){
				dNMI = 2.0 * *(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 0) / denom;
			}

			*(dNMIResults + idTh * ((maxValueDatasetGPU + maxValueDatasetGPU) + (maxValueDatasetGPU * maxValueDatasetGPU) + 3) + 0) = dNMI;
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
		double *aResultKendalls, *aResultSpearmans, *fSpearmanRankG1, *fSpearmanRankG2, *fSpearmanStats;
		float *dNMIResults;
		unsigned long long int numResultKendallsCpu = 0;
		unsigned long long int numResultSpearmanCpu = 0;
		unsigned long long int numResultNMICpu = 0;
		cudaMalloc((void **)&aResultKendalls, pairsPerRun * sizeof(double));
		cudaMemset(aResultKendalls, 0, pairsPerRun * sizeof(double));
		cudaMalloc((void **)&aResultSpearmans, pairsPerRun * sizeof(double));
		cudaMemset(aResultSpearmans, 0, pairsPerRun * sizeof(double));
		cudaMalloc((void **)&fSpearmanRankG1, pairsPerRun * ulColsData * sizeof(double));
		cudaMemset(fSpearmanRankG1, 0, pairsPerRun * ulColsData * sizeof(double));
		cudaMalloc((void **)&fSpearmanRankG2, pairsPerRun * ulColsData * sizeof(double));
		cudaMemset(fSpearmanRankG2, 0, pairsPerRun * ulColsData * sizeof(double));
		cudaMalloc((void **)&fSpearmanStats, pairsPerRun * 4 * sizeof(double));
		cudaMemset(fSpearmanStats, 0, pairsPerRun * 4 * sizeof(double));
		cudaMalloc((void **)&dNMIResults, pairsPerRun * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) * sizeof(float));
		cudaMemset(dNMIResults, 0, pairsPerRun * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) * sizeof(float));

		// Kendall
		ulong lCombination = 0;
		for(ulong i=1; i < ulColsData; i++){
			lCombination += i;
		}
		prepareGpu1D(pairsPerRun);
		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			kendallTwoGenes<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(lCombination, aResultKendalls, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		kendallTwoGenes<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(lCombination, aResultKendalls, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		// Spearman
		prepareGpu1D(pairsPerRun);
		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCalcRanks<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		spearmanCalcRanks<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCovariance<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, aResultSpearmans, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		spearmanCovariance<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, aResultSpearmans, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

		for (int i = 1; i <= maxIteratorGPU; i++)
		{
			spearmanCalcValue<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, aResultSpearmans, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		spearmanCalcValue<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(fSpearmanStats, aResultSpearmans, fSpearmanRankG1, fSpearmanRankG2, maxPairs, id, pairsPerGpuPrevious, mDataGPU, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);
		
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
			nmiCalcValue<<<maxBlocksPerGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, totalPairs, pairsPerRun, i, maxThreadsPerBlock * maxBlocksPerGrid);
		}
		nmiCalcValue<<<lastBlocksGrid, maxThreadsPerBlock, 0, s>>>(dNMIResults, maxPairs, id, pairsPerGpuPrevious, totalPairs, pairsPerRun, maxIteratorGPU + 1, maxThreadsPerBlock * maxBlocksPerGrid);

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
			
			float *dNMIResultsCpu = (float *)malloc(pairsPerRun * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) * sizeof(float));
			cudaMemcpy(dNMIResultsCpu, dNMIResults, pairsPerRun * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpyFromSymbol(&numResultNMICpu, numResultNMI,
								 sizeof(unsigned long long int), 0, cudaMemcpyDeviceToHost);

			if(largeScale == 0){
				fileResults << "g1;g2;nmi;kendall;spearman"
						<< "\n";
			}
			for (ulong iRow = 0; iRow < pairsPerRun; iRow++)
			{
				ulong pattern = iRow + pairsPerGpuPrevious;

				if (pattern < maxPairs)
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
					float dNMI = *(dNMIResultsCpu + iRow * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) + 0);

					if(dKendall < 0){
						dKendall = dKendall * -1;
					}

					if(dSpearman < 0){
						dSpearman = dSpearman * -1;
					}

					if(dNMI < 0){
						dNMI = dNMI * -1;
					}

					if(dKendall >= correctionThresholdKendall){
						iMajorVoting++;
					}
	
					if(dSpearman >= correctionThresholdSpearman){
						iMajorVoting++;
					}

					if(dNMI >= correctionThresholdNMI){
						iMajorVoting++;
					}
					
					if (iMajorVoting >= 2)
					{
						fileResults << geneNames[r1] << ";" << geneNames[r2] << ";" << dNMI << ";" << dKendall << ";" << dSpearman << "\n";
					}
				}
			}

			free(aResultKendallsCpu);
			free(aResultSpearmansCpu);
			free(dNMIResultsCpu);
		}

		cudaFree(aResultKendalls);
		cudaFree(aResultSpearmans);
		cudaFree(fSpearmanRankG1);
		cudaFree(fSpearmanRankG2);
		cudaFree(fSpearmanStats);
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
		cudaMemcpyToSymbol(*(&maxValueDatasetGPU), &maxValueDataset, sizeof(int), 0,
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
		double sizeResultKendall = 0, sizeResultsSpearman = 0, fSpearmanRankG1 = 0, fSpearmanRankG2 = 0, fSpearmanStats = 0;
		float dNMIResults = 0;
		sizeResultKendall = (pairsPerGpu * sizeof(double));
		sizeResultsSpearman = (pairsPerGpu * sizeof(double));
		fSpearmanRankG1 = (pairsPerGpu * ulColsData * sizeof(double));
		fSpearmanRankG2 = (pairsPerGpu * ulColsData * sizeof(double));
		fSpearmanStats = (pairsPerGpu * 4 * sizeof(double)); // [1] --> Mean of rank G1; [2] --> Mean of rank G2; ; [3] --> aux ; [4] --> aux
		dNMIResults = (pairsPerGpu * ((maxValueDataset + maxValueDataset) + (maxValueDataset * maxValueDataset) + 3) * sizeof(float)); // [0] --> Mutual information ::: NMI [1] --> entropyGen1 ::: [2] --> entropyGen2 ::: [3 - maxValueDataset] --> probMaps
		chunks[i] = ((sizeResultKendall + sizeResultsSpearman + fSpearmanRankG1 + fSpearmanRankG2 + fSpearmanStats + dNMIResults) / availableMemory) + 1;
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

