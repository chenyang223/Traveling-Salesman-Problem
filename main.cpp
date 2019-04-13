#include "ga.h"

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

string &GetDataFile()
{
    static string data_file;
    return data_file;
}

string &GetDataPath()
{
    static string data_path;
    return data_path;
}