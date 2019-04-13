#include "mtsp_kmeans_ga.h"
#include "utils.h"
#include <iostream>

using namespace std;
int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cout << "Wrong command line parametes!" << endl;
        system("pause");
        return 0;
    }

    srand(unsigned(time(0)));
    time_t t_begin = clock();

    GetDataPath() = argv[1];
    GetDataFile() = argv[2];

    MtspKmeansGA();

    time_t t_end = clock();
    double running_time = double(t_end - t_begin) / CLOCKS_PER_SEC;
    cout << "\n"
         << "The running time is:  " << running_time << endl;
    return 0;
}