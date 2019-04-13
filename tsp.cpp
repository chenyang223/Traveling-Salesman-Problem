#include "ga.h"

typedef struct
{
    double path_length;
    vector<int> path;
    double pr_roulette;
} TspSolution;

const int kNumInGroup = 30;
const int kNumCross = 15;

const double kPrMutation = 0.1;
const double kPrCrossover = 0.95;
const double kPrDisaster = 0.01;
const int kNumIteration = 1000;
const int kTspMode = 1;
// 0:for simple tsp  1: collaboration and average minimum

static void PrintTspSolutionVector(vector<TspSolution> s);
static void PrintTspSolution(TspSolution s);
static void CreateNetworkRandom(CityNetwork &network);
static void InitializeGroup(CityNetwork network, vector<TspSolution> &group);
static double CalculatePathLength(vector<int> path, CityNetwork network);
static void CalculatePrRoulette(vector<TspSolution> &group);
static void UpdatePathLength(CityNetwork n, vector<TspSolution> &g);
static int Evolution(CityNetwork network, vector<TspSolution> &group);
static int RouletteSelection(vector<TspSolution> group);
static void Crossover(TspSolution &father, TspSolution &mother, vector<TspSolution> &son);
static void Mutation(TspSolution &s);
static void UpdateGroup(vector<TspSolution> &group, vector<TspSolution> &son);
static bool TspSortFunction(TspSolution a, TspSolution b);
static int FindBest(vector<TspSolution> group);
static void OutputBestSolution(vector<TspSolution> group, int best);


void TspGA()
{
    CityNetwork network;
    vector<TspSolution> group(kNumInGroup);
    int best_solution_index;

    network = ConstructNetwork();
    InitializeGroup(network, group);
    PrintTspSolutionVector(group);
    best_solution_index = Evolution(network, group);
    cout << "the final result:" << endl;
    PrintTspSolution(group.at(best_solution_index));
    OutputBestSolution(group, best_solution_index);
}

static void OutputBestSolution(vector<TspSolution> group, int best)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "path_" << kTspMode << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    switch(kTspMode)
    {
    case 0:
        for (vector<int>::iterator i = group.at(best).path.begin();
             i != group.at(best).path.end(); i++)
            output << (*i) << " ";
        break;
    case 1:
        int middle = group.at(best).path.size()/2;
        for(int i = 0; i < middle; i++)
            output << group.at(best).path.at(i) << " ";
        output << "\n";
        for(int i = 0; i < middle; i++)
            output << group.at(best).path.at(i+middle) << " ";
            break;
    }
    output.close();
}

static void CalculatePrRoulette(vector<TspSolution> &group)
{
    double overall = 0.0;
    for (int i = 0; i < group.size(); i++)
        overall += 1.0 / group[i].path_length;
    for (int i = 0; i < group.size(); i++)
        group[i].pr_roulette = 1.0 / group[i].path_length / overall;
}

static bool TspSortFunction(TspSolution a, TspSolution b)
{
    return (a.path_length < b.path_length);
}

static void PrintTspSolutionVector(vector<TspSolution> s)
{
    cout << "num of members: " << s.size() << "\n"
         << endl;
    for (int i = 0; i < s.size(); i++)
    {
        cout << i << "  length: " << s[i].path_length << " pr: " << s[i].pr_roulette << endl;
        cout << "path: ";
        for (int j = 0; j < s[i].path.size(); j++)
            cout << s[i].path[j] << "->";
        cout << "\n"
             << endl;
    }
}


static void PrintVecInt(vector<int> v)
{
    for (vector<int>::iterator i = v.begin(); i != v.end(); i++)
        cout << *i << " ";
    cout << endl;
}

static void PrintTspSolution(TspSolution s)
{
    cout << "\n"
         << "path length: " << s.path_length;
    cout << "path: ";
    for (int j = 0; j < s.path.size(); j++)
        cout << s.path[j] << "->";
    cout << "\n"
         << endl;
}

static double CalculatePathLength(vector<int> path, CityNetwork network)
{
    double total_length = 0.0;
    switch (kTspMode)
    {
    case 0:
        for (int i = 0; i < path.size() - 1; i++)
            total_length += network.adjacent_matrix[path[i]][path[i + 1]];
        total_length += network.adjacent_matrix[path.back()][path.front()];
        break;
    case 1:
        for (int i = 0; i < path.size() - 1; i++)
            total_length += network.adjacent_matrix[path[i]][path[i + 1]];
        int middle = path.size() / 2;
        total_length -= network.adjacent_matrix[path.at(middle - 1)][path.at(middle)];
        break;
    }
    return total_length;
}

static void UpdatePathLength(CityNetwork n, vector<TspSolution> &g)
{
    for (vector<TspSolution>::iterator i = g.begin(); i != g.end(); i++)
        (*i).path_length = CalculatePathLength((*i).path, n);
}

static void InitializeGroup(CityNetwork network, vector<TspSolution> &group)
{
    double total_path_length = 0.0;
    vector<int> base_path;
    for (int i = 0; i < network.num_nodes; i++)
        base_path.push_back(i);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine rde(seed);

    for (int i = 0; i < group.size(); i++)
    {
        shuffle(base_path.begin(), base_path.end(), rde);
        group[i].path = base_path;
        group[i].path_length = CalculatePathLength(group[i].path, network);
    }
    CalculatePrRoulette(group);
}

static int Evolution(CityNetwork network, vector<TspSolution> &group)
{
    vector<TspSolution> son;
    
    for (int num_iter = 0; num_iter < kNumIteration; num_iter++)
    {
        for (int i = 0; i < kNumCross; i++)
        {
            double pr_cross = rand() / double(RAND_MAX);
            int father_index = RouletteSelection(group);
            int mother_index = RouletteSelection(group);
            if (pr_cross > kPrCrossover)
            {
                i--;
                continue;
            }
            TspSolution father = group[father_index];
            TspSolution mother = group[mother_index];
            Crossover(father, mother, son);
        }

        for (vector<TspSolution>::iterator i = son.begin(); i != son.end(); i++)
        {
            double pr_mutation = rand() / double(RAND_MAX);
            if (pr_mutation > kPrMutation)
                continue;
            Mutation(*i);
        }
        UpdatePathLength(network, son);
        UpdateGroup(group, son);
        CalculatePrRoulette(group);
        son.clear();
    }
    PrintTspSolutionVector(group);
    int best = FindBest(group);
    return best;
}

static int RouletteSelection(vector<TspSolution> group)
{
    double pr_selection = rand() / double(RAND_MAX + 1);
    double pr_accumulate = 0.0;
    for (int i = 0; i < group.size(); i++)
    {
        pr_accumulate += group[i].pr_roulette;
        if (pr_selection < pr_accumulate)
            return i;
    }
    cout << "error: no result from roulette selection" << endl;
    return 0;
}

static void Crossover(TspSolution &father, TspSolution &mother, vector<TspSolution> &son)
{
    int num_city = father.path.size();
    int cross_index1, cross_index2;
    do
    {
        cross_index1 = rand() % (num_city);
        cross_index2 = rand() % (num_city);
        if (cross_index1 > cross_index2)
            swap(cross_index1, cross_index2);
    } while (cross_index1 == 0 && cross_index2 == num_city - 1);

    vector<int> father_cross(father.path.begin() + cross_index1,
                             father.path.begin() + cross_index2 + 1);
    vector<int> mother_cross(mother.path.begin() + cross_index1,
                             mother.path.begin() + cross_index2 + 1);
    vector<int> common_cross;
    vector<int> fa_cross_sort;
    vector<int> mo_cross_sort;
    fa_cross_sort.assign(father_cross.begin(), father_cross.end());
    mo_cross_sort.assign(mother_cross.begin(), mother_cross.end());
    sort(fa_cross_sort.begin(), fa_cross_sort.end());
    sort(mo_cross_sort.begin(), mo_cross_sort.end());
    set_intersection(fa_cross_sort.begin(), fa_cross_sort.end(),
                     mo_cross_sort.begin(), mo_cross_sort.end(),
                     inserter(common_cross, common_cross.begin()));
    //PrintVecInt(common_cross);
    vector<int> father_exchange;
    vector<int> mother_exchange;
    set_difference(mo_cross_sort.begin(), mo_cross_sort.end(),
                   common_cross.begin(), common_cross.end(),
                   inserter(father_exchange, father_exchange.begin()));
    set_difference(fa_cross_sort.begin(), fa_cross_sort.end(),
                   common_cross.begin(), common_cross.end(),
                   inserter(mother_exchange, mother_exchange.begin()));
    vector<int>::iterator change;
    for (int i = 0; i < father_exchange.size(); i++)
    {
        change = find(father.path.begin(), father.path.end(), father_exchange[i]);
        *change = mother_exchange[i];
        change = find(mother.path.begin(), mother.path.end(), mother_exchange[i]);
        *change = father_exchange[i];
    }
    copy(mother_cross.begin(), mother_cross.end(), father.path.begin() + cross_index1);
    copy(father_cross.begin(), father_cross.end(), mother.path.begin() + cross_index1);
    //PrintVecInt(father.path);
    //PrintVecInt(mother.path);
    son.push_back(father);
    son.push_back(mother);
}

static void Mutation(TspSolution &s)
{
    int num_city = s.path.size();
    int cityindex1 = rand() % num_city;
    int cityindex2 = rand() % num_city;
    while (cityindex1 == cityindex2)
        cityindex2 = rand() % num_city;
    swap(s.path[cityindex1], s.path[cityindex2]);
}

static void UpdateGroup(vector<TspSolution> &group, vector<TspSolution> &son)
{
    if (group.size() > son.size())
        cout << "cross error: size of son is smaller than size of group";
    sort(son.begin(), son.end(), TspSortFunction);
    int original_best = FindBest(group);
    if (original_best != 0)
        swap(group.at(0), group.at(original_best));
    double pr_disaster = rand() / double(RAND_MAX);
    if (pr_disaster < kPrDisaster)
        copy(son.begin(), son.begin() + group.size(), group.begin());
    else
        copy(son.begin(), son.begin() + group.size() - 1, group.begin() + 1);
}

static int FindBest(vector<TspSolution> group)
{
    int best = 0;
    double best_path_length = group.at(best).path_length;
    for (int i = 1; i < group.size(); i++)
    {
        if (group.at(i).path_length < best_path_length)
        {
            best = i;
            best_path_length = group.at(i).path_length;
        }
    }
    return best;
}