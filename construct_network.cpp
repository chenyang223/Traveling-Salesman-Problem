#include "utils.h"
#include "construct_network.h"
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
#include <sstream>
#include <fstream>
#include <cfloat>

using namespace std;
// For function: CreateNetworkRandom
const int kNumCity = 10;

CityNetwork ConstructNetwork()
{
    CityNetwork network;
    CreateNetworkReadData(network);
    //PrintNetwork(network);
    return network;
}

void PrintNetwork(CityNetwork n)
{
    cout << "num of cities: " << n.num_nodes << "\n";
    for (int i = 0; i < n.num_nodes; i++)
    {
        cout << setw(7) << setiosflags(ios::left);
        cout << i;
        cout << setw(7) << setiosflags(ios::left);
        cout << n.nodes_coordinate[i][0];
        cout << setw(7) << setiosflags(ios::left);
        cout << n.nodes_coordinate[i][1] << "\n";
    }
    cout << "\n";

    for (int i = 0; i < n.num_nodes; i++)
    {
        for (int j = 0; j < n.num_nodes; j++)
        {
            cout << setw(7) << setiosflags(ios::right);
            cout << n.adjacent_matrix[i][j] << " ";
        }
        cout << "\n";
    }

    cout << endl;
}

void CreateNetworkReadData(CityNetwork &network)
{
    ifstream read_in;
    read_in.open(GetDataPath() + GetDataFile());
    if (!read_in.is_open())
    {
        cout << "fail to open file." << endl;
        return;
    }
    string s;
    do
    {
        read_in >> s;
    } while (s != "DIMENSION:");
    read_in >> network.num_nodes;
    cout << "num of cities: " << network.num_nodes << endl;
    do
    {
        read_in >> s;
    } while (s != "EDGE_WEIGHT_TYPE:");
    read_in >> s;
    if (s != "EUC_2D")
    {
        cout << "edge weight type error: not euc_2d type." << endl;
        return;
    }
    do
    {
        read_in >> s;
    } while (s != "NODE_COORD_SECTION");
    network.nodes_coordinate.resize(network.num_nodes);
    network.adjacent_matrix.resize(network.num_nodes);
    for (int i = 0; i < network.num_nodes; i++)
    {
        int index;
        read_in >> index;
        if (index != i + 1)
        {
            cout << "mistake in reading coordinates" << endl;
            return;
        }
        double receiver;
        read_in >> receiver;
        (network.nodes_coordinate.at(i)).push_back(receiver);
        read_in >> receiver;
        (network.nodes_coordinate.at(i)).push_back(receiver);
        network.adjacent_matrix.at(i).resize(network.num_nodes);
        network.adjacent_matrix.at(i).at(i) = 0.0;
    }
    read_in >> s;
    if (s == "EOF")
        cout << "reading complete " << endl;

    for (int i = 0; i < network.num_nodes; i++)
        for (int j = 0; j < i; j++)
        {
            double distance = pow((network.nodes_coordinate[i][0] - network.nodes_coordinate[j][0]), 2) +
                              pow((network.nodes_coordinate[i][1] - network.nodes_coordinate[j][1]), 2);
            distance = sqrt(distance);
            network.adjacent_matrix[i][j] = distance;
            network.adjacent_matrix[j][i] = distance;
        }
}

void CreateNetworkRandom(CityNetwork &network)
{
    network.num_nodes = kNumCity;
    network.nodes_coordinate.resize(kNumCity);
    network.adjacent_matrix.resize(kNumCity);
    for (int i = 0; i < network.num_nodes; i++)
    {
        network.nodes_coordinate[i].push_back(rand() / double(RAND_MAX));
        network.nodes_coordinate[i].push_back(rand() / double(RAND_MAX));
        network.adjacent_matrix[i].resize(kNumCity);
        network.adjacent_matrix[i][i] = 0.0;
    }
    for (int i = 0; i < network.num_nodes; i++)
        for (int j = 0; j < i; j++)
        {
            double distance = pow((network.nodes_coordinate[i][0] - network.nodes_coordinate[j][0]), 2) +
                              pow((network.nodes_coordinate[i][1] - network.nodes_coordinate[j][1]), 2);
            distance = sqrt(distance);
            network.adjacent_matrix[i][j] = distance;
            network.adjacent_matrix[j][i] = distance;
        }
}
