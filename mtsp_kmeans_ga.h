#ifndef _GA_H_
#define _GA_H_

#include "construct_network.h"
#include <cfloat>
typedef struct
{
    double length;
    vector<int> path;
} TspChromosome;

void MtspKmeansGA();
const int kBaseNumInGroup = 20;
const double kDoubleLargeNumber = DBL_MAX;
const int kBaseGaSteps = 100;

vector<vector<int>> ReadKMeansCluster();
void InitializeGroup(CityNetwork &network, vector<TspChromosome> &group,
                            vector<vector<int>>::iterator base_path);
void PrintMtspResult(vector<vector<int>> v, double length);
void OutputBestSolution(vector<vector<int>> v, double length, int num_salesman);
void OutputInitSolution(vector<vector<int>> v, double length, int num_salesman);
void OutputBestInGeneration(int num_salesman);
vector<int>::iterator FindGreedy(CityNetwork &network,
                                        vector<int>::iterator begin,
                                        vector<int>::iterator end, int start_point);
void GAMain(CityNetwork &network, vector<TspChromosome> &group);
bool LKImprovePath(CityNetwork &network, vector<int> path, int depth);
void LocalSearchPath(CityNetwork &network, vector<int> &path);
double CalculatePathLength(CityNetwork &network, vector<int> path);
vector<int> DPXCrossover(CityNetwork &network, vector<int> father, vector<int> mother);
vector<int> &GetLKImprovedPath();
vector<double> &GetBestAtOneGeneration();
int &GetNumInGroup();
int &GetGASteps();

#endif
