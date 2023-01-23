#define EXTERN
#include "main.h"
using namespace std;

int main(int argc, char **argv)
{
    maxValueDataset = -1;
    if (!getParameters(argc, argv))
    {
        preprocess();
        if (maxValueDataset != -1)
        {
            if (bOutput)
            {
                mkdir("results", 0777);
            }
        }

        double elapsed = runAlgorithm();

        // Save results
        //cout << "Resume: " << endl
        //     << "========================" << endl;
        //cout << "Dataset: " << sDataPath << endl;
        //cout << "Output: " << bOutput << endl;
        //cout << "GPU Devices: " << iDeviceCount << endl;
        //cout << "Time measured (seconds): " << elapsed << endl;
    }

    return 0;
}