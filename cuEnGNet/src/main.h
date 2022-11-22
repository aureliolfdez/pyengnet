#ifndef EXTERN
#define EXTERN extern
#endif

#include <thread>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <sys/time.h>
#include <unistd.h>
#include <iterator>
#include <utility>
#include <string>
#include <cstdlib>
#include <time.h>
#include <set>
#include <map>
#include <unordered_set>
#include <mutex>
#include <cctype>
#include <string.h>

using namespace std;

EXTERN string sDataPath, sDelimiter;
EXTERN int iHeaderColumns, iDeviceCount;
EXTERN bool bOutput;
EXTERN double correctionThresholdSpearman, correctionThresholdKendall, correctionThresholdNMI;
EXTERN float *mData;
EXTERN int *mDataNormalized;
EXTERN ulong ulRowsData, ulColsData, patternSize,maxPairs, maxThreadsPerBlock, maxBlocksPerGrid, maxIteratorGPU, lastBlocksGrid;

// Preprocess
bool preprocess();
bool readerData();
void normalizedData();
void countColumnData(string sLine);
bool getParameters(int argc, char **argv);
void printHelp();

// CUDA
double runAlgorithm();