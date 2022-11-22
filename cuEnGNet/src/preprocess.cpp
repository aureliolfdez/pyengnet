#include "main.h"
#include <getopt.h>
using namespace std;

bool preprocess()
{
    bool bSuccess = readerData();
    if (bSuccess)
    {
        normalizedData();
    }
    return bSuccess;
}

bool readerData()
{
    bool bSuccess = false;
    ifstream inputFile(sDataPath.c_str());
    if (inputFile.is_open())
    {
        vector<string> rowsArray_Aux;
        string line;
        ulRowsData = ulColsData = 0;

        // 1) Remove first line
        getline(inputFile, line);

        // 2) Get number of rows and cols
        while (getline(inputFile, line))
        {
            if (ulRowsData == 0)
            {
                countColumnData(line);
            }
            rowsArray_Aux.push_back(line);
            ulRowsData++;
        }
        inputFile.close();

        // 3) Create input matrix
        mData = (float *)malloc(ulRowsData * ulColsData * sizeof(float));

        // 4) Fill input data matrix
        for (unsigned long ulCountRow = 0; ulCountRow < ulRowsData; ulCountRow++)
        {
            unsigned long ulCountCol = 1;
            size_t pos = 0;
            string line = rowsArray_Aux[ulCountRow];
            string token = line.substr(0, line.find(sDelimiter));
            line.erase(0, line.find(sDelimiter) + sDelimiter.length());

            do
            {
                pos = line.find(sDelimiter);
                token = line.substr(0, pos);
                line.erase(0, pos + sDelimiter.length());
                *(mData + ulCountRow * ulColsData + (ulCountCol - 1 - iHeaderColumns)) = stof(token);
                ulCountCol++;
            } while (pos != std::string::npos);
        }
        bSuccess = true;
    }
    return bSuccess;
}

void countColumnData(string sLine)
{
    size_t pos = 0;
    string token;
    int col = 0;

    do
    {
        pos = sLine.find(sDelimiter);
        token = sLine.substr(0, pos);
        sLine.erase(0, pos + sDelimiter.length());
        if (col > iHeaderColumns)
        {
            ulColsData++;
        }
        col++;
    } while (pos != std::string::npos);
}

void normalizedData()
{
    mDataNormalized = (int *)malloc(ulRowsData * (ulColsData + 1) * sizeof(int));
    for (unsigned long ulCountRow = 0; ulCountRow < ulRowsData; ulCountRow++)
    {
        int maxValue = (int)floor(*(mData + ulCountRow * ulColsData + 0));
        int minValue = (int)floor(*(mData + ulCountRow * ulColsData + 0));
        for (unsigned long ulCountCol = 0; ulCountCol < ulColsData; ulCountCol++)
        {
            int iExp = (int)floor(*(mData + ulCountRow * ulColsData + ulCountCol));
            *(mDataNormalized + ulCountRow * (ulColsData + 1) + ulCountCol) = iExp;
            if (iExp < minValue)
            {
                minValue = iExp;
            }

            if (iExp > maxValue)
            {
                maxValue = iExp;
            }
        }

        for (unsigned long ulCountCol = 0; ulCountCol < ulColsData; ulCountCol++)
        {
            *(mDataNormalized + ulCountRow * (ulColsData + 1) + ulCountCol) = *(mData + ulCountRow * ulColsData + ulCountCol) - minValue;
        }

        maxValue = maxValue - minValue + 1;
        *(mDataNormalized + ulCountRow * (ulColsData + 1) + ulColsData) =  maxValue;
    }
}

bool getParameters(int argc, char **argv)
{
    bool bError = false;
    int iFieldsRequired = 0;
    int c;

    sDelimiter = ",";
    bOutput = false;
    correctionThresholdSpearman = 0;
    correctionThresholdKendall = 0;
    correctionThresholdNMI = 0;
    iDeviceCount = 1;
    sDataPath = "";

    const char *const short_opts = "a:b:c:d:e:f:!:?";
    const option long_opts[] = {
        {"out", no_argument, nullptr, 'a'},   // Output
        {"cors", required_argument, nullptr, 'b'},   // Correction threshold Spearman
        {"gpu", required_argument, nullptr, 'c'},   // Device counts
        {"data", required_argument, nullptr, 'd'},  // Datapath
        {"cork", required_argument, nullptr, 'e'},   // Correction threshold Kendall
        {"corn", required_argument, nullptr, 'f'},   // Correction threshold NMI
        {"help", no_argument, nullptr, '!'},      // Help
        {nullptr, no_argument, nullptr, 0}};

    while ((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1)
    {
        switch (c)
        {
        case 'a':
            bOutput = true;
            break;
        case 'b':
            correctionThresholdSpearman = stof(optarg);
            break;
        case 'c':
            iDeviceCount = atoi(optarg);
            break;
        case 'd':
            sDataPath = optarg;
            iFieldsRequired++;
            break;
        case 'e':
            correctionThresholdKendall = stof(optarg);
            break;
        case 'f':
            correctionThresholdNMI = stof(optarg);
            break;
        case '!': // Help
        case ':':
            bError = true;
            break;
        case '?': // Illegal options
            cout << "Illegal option: -" << (char)optopt << endl;
            bError = true;
            break;
        }
    }

    if (bError)
    {
        printHelp();
    }
    else if (iFieldsRequired != 1)
    {
        cout << "Information must be provided for required parameters." << endl;
        cout << "--------------------------------------------------------" << endl
             << endl;
        printHelp();
        bError = true;
    }

    return bError;
}

void printHelp()
{
    cout << "Usage: ./cuEnGNet [OPTIONS]" << endl;
    cout << "--data <path>: Input dataset file [Required]" << endl;
    cout << "--cors <number>: Correction threshold Spearman value [Optional. Default: 0]" << endl;   
    cout << "--cork <number>: Correction threshold Kendall value [Optional. Default: 0]" << endl;   
    cout << "--corn <number>: Correction threshold NMI value [Optional. Default: 0]" << endl;   
    cout << "--gpu <number>: Number of GPUs to run [Optional. Default: 1]" << endl;   
    cout << "--out: Generate output file [Optional. Default: false]" << endl;
}
