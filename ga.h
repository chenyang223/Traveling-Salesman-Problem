#ifndef _GA_H_
#define _GA_H_
    #include <iostream>
    #include <cstdlib>
    #include <ctime>
    #include <algorithm>
    #include <vector>
    #include <cmath>
    #include <random> 
    #include <chrono>
    #include <iterator>
    #include <iomanip>
    #include <string>
    #include <sstream>
    #include <fstream>
    #include <cfloat>

    using namespace std;

    typedef struct
    {
        int num_nodes;
        vector<vector<double>> nodes_coordinate;
        vector<vector<double>> adjacent_matrix;
    } CityNetwork;

    //construct_network.cpp
    CityNetwork ConstructNetwork();
    void PrintNetwork(CityNetwork n);

    //tsp.cpp
    void TspGA();

    //mtsp.cpp
    void MtspMOGA();

    //mtsp_kmeans_moga.cpp
    void MtspKmeansMOGA();

    // tsp_kmeans_ga.cpp
    void MtspKmeansGA();
    
    //main.cpp
    string &GetDataFile();
    string &GetDataPath();
#endif