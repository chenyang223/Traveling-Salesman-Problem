
#ifndef TRAVEL_SALESMAN_PROBLEM_CONSTRUCT_NETWORK_H
#define TRAVEL_SALESMAN_PROBLEM_CONSTRUCT_NETWORK_H

#include <vector>
using namespace std;

typedef struct
{
    int num_nodes;
    vector<vector<double>> nodes_coordinate;
    vector<vector<double>> adjacent_matrix;
} CityNetwork;

CityNetwork ConstructNetwork();
void PrintNetwork(CityNetwork n);
void CreateNetworkReadData(CityNetwork &network);
void CreateNetworkRandom(CityNetwork &network);

#endif //TRAVEL_SALESMAN_PROBLEM_CONSTRUCT_NETWORK_H
