#include "ga.h"

typedef struct
{
    double fitness_func;
    vector<int> path;
} TspChromosome;

const int kNumInGroup = 10;
const double kDoubleLargeNumber = DBL_MAX;
const int kDepotIndex = 0;

static vector<vector<int>> ReadKMeansCluster();
// no
static void InitializeGroup(CityNetwork network, vector<TspChromosome> &group,
                     vector<vector<int>>::iterator base_path);

static vector<int>::iterator FindGreedy(CityNetwork network,
                                 vector<int>::iterator begin,
                                 vector<int>::iterator end, int start_point);
static void MOGAMain(CityNetwork network, vector<TspChromosome> &group);
static void PrintTspChromosomeVector(vector<TspChromosome> v);
static vector<vector<int>> GenerateSwappingPath(int length);
static void TabuSearchPath(CityNetwork network, vector<TspChromosome> &group,
                    vector<TspChromosome>::iterator target, int fitness_mode,
                    vector<vector<int>>::iterator swap_begin,
                    vector<vector<int>>::iterator swap_end);
static double TabuSearchFitness(CityNetwork network, vector<TspChromosome> &group,
                         vector<TspChromosome>::iterator target,
                         vector<vector<int>>::iterator swapping,
                         int fitness_mode);
static double CalculatePathFitness(CityNetwork network,
                            vector<TspChromosome>::iterator target);
static void PrintMtspResult(vector<vector<int>> v, double length);
static void OutputBestSolution(vector<vector<int>> v, double length);

void MtspKmeansMOGA()
{
    CityNetwork network;
    network = ConstructNetwork();
    //PrintNetwork(network);
    vector<vector<int>> kmeans_cluster;
    kmeans_cluster = ReadKMeansCluster();
    vector<vector<int>> path_result;
    double total_path_length = 0.0;
    
    for (vector<vector<int>>::iterator i = kmeans_cluster.begin();
         i != kmeans_cluster.end(); i++)
    {
        vector<TspChromosome> group(kNumInGroup);
        InitializeGroup(network, group, i);
        MOGAMain(network, group);
        path_result.push_back(group.begin()->path);
        total_path_length += group.begin()->fitness_func;
    }
    PrintMtspResult(path_result,total_path_length);  
    OutputBestSolution(path_result,total_path_length);  
}

static void PrintMtspResult(vector<vector<int>> v, double length)
{
    cout << "num of salesman: " << v.size() << "\n";
    cout << "path length: " << length << "\n";
    for (vector<vector<int>>::iterator i = v.begin(); i != v.end(); i++)
    {
        for(vector<int>::iterator j = (*i).begin(); j != (*i).end();j++)
            cout << *j << " ";
        cout << "\n";
    }
    cout << endl;
}

static void OutputBestSolution(vector<vector<int>> v, double length)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "mtsp_path_kmeansMoga";
    filename_stream << "_" << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    output << "NUM_SALESMAN: " << v.size() << "\n";
    output << "PATH_LENGTH: " << length << "\n";
    for (vector<vector<int>>::iterator i = v.begin(); i != v.end(); i++)
    {
        output << kDepotIndex << " ";
        for(vector<int>::iterator j = (*i).begin(); j != (*i).end();j++)
            output << *j << " ";
        output << kDepotIndex << " ";
        output << "\n";
    }
    output.close();
}

static void InitializeGroup(CityNetwork network, vector<TspChromosome> &group,
                     vector<vector<int>>::iterator base_path_iter)
{
    int path_sequence_length = (*base_path_iter).size();
    for (vector<TspChromosome>::iterator i = group.begin();
         i != group.end(); i++)
    {
        vector<int> base_path = *base_path_iter;
        double fitness = 0;
        int start_point;
        vector<int>::iterator path_iterator = base_path.begin();
        vector<int>::iterator greedy_iterator;
        
        start_point = rand() % path_sequence_length;
        swap((*path_iterator), base_path.at(start_point));
        path_iterator++;
        while (path_iterator != base_path.end())
        {
            greedy_iterator = FindGreedy(network, path_iterator, base_path.end(),
                                         *(path_iterator - 1));
            swap(*greedy_iterator, *path_iterator);
            fitness += network.adjacent_matrix
                           [*(path_iterator - 1)][*path_iterator];
            path_iterator++;
        }
        fitness += network.adjacent_matrix
                       [*(path_iterator - 1)][start_point];
        i->path.assign(base_path.begin(), base_path.end());
        i->fitness_func = fitness;
    }
}

static vector<int>::iterator FindGreedy(CityNetwork network,
                                 vector<int>::iterator begin,
                                 vector<int>::iterator end, int start_point)
{
    vector<int>::iterator greedy_iterator = begin;
    double min_path_length = network.adjacent_matrix[start_point][*begin];
    for (vector<int>::iterator i = begin; i != end; i++)
    {
        double path_len = network.adjacent_matrix[start_point][*i];
        if (path_len < min_path_length)
        {
            min_path_length = path_len;
            greedy_iterator = i;
        }
    }
    return greedy_iterator;
}

static vector<vector<int>> ReadKMeansCluster()
{
    ifstream read;
    read.open(GetDataPath() + "clusters_" + GetDataFile());
    if (!read.is_open())
    {
        cout << "fail to open file." << endl;
    }
    string s;
    stringstream stream;
    int i;
    vector<int> vec;
    vector<vector<int>> kmeans_clusters;

    do
    {
        read >> s;
    } while (s != "Cluster:");
    read >> s;
    do
    {
        do
        {
            stream.clear();
            stream << s;
            stream >> i;
            vec.push_back(i);
            read >> s;
        } while (s != "Cluster:" && s != "EOF");
        kmeans_clusters.push_back(vec);
        vec.clear();
        read >> s;
    } while (s != "EOF");

    return kmeans_clusters;
}

static void PrintTspChromosomeVector(vector<TspChromosome> v)
{
    cout << "num of members: " << v.size() << "\n";
    for (int i = 0; i < v.size(); i++)
    {
        cout << i << "  fitness: " << v.at(i).fitness_func << "\n";
        cout << "path: "
             << "\n";
        for (int j = 0; j < v.at(i).path.size(); j++)
            cout << v.at(i).path[j] << "->";
        cout << "\n"
             << endl;
    }
}

static void MOGAMain(CityNetwork network, vector<TspChromosome> &group)
{
    struct MySort
    {
        bool operator()(TspChromosome i, TspChromosome j)
        {
            return (i.fitness_func < j.fitness_func);
        }
    } sort_object;

    int path_sequence_length = group.begin()->path.size();
    vector<vector<int>> swapping_path =
        GenerateSwappingPath(path_sequence_length);

    vector<double> prob_list;
    for (int i = 0; i < group.size(); i++)
        prob_list.push_back(i / (double)(group.size() - 1.0));

    int num_step_moga = 20;
    for (int i = 0; i < num_step_moga; i++)
    {
        sort(group.begin(), group.end(), sort_object);
        for (int j = 0; j < group.size(); j++)
        {
            int fitness_mode;
            if (rand() / (double)RAND_MAX > prob_list.at(j))
                fitness_mode = 0;
            else
                fitness_mode = 1;
            TabuSearchPath(network, group, group.begin() + j,
                           fitness_mode, swapping_path.begin(),
                           swapping_path.end());
            (group.begin() + j)->fitness_func =
                CalculatePathFitness(network, group.begin() + j);
        }
    }
    sort(group.begin(), group.end(), sort_object);
}

static double CalculatePathFitness(CityNetwork network,
                            vector<TspChromosome>::iterator target)
{
    double path_length = 0;
    vector<int>::iterator path_iter = target->path.begin();
    path_length += network.adjacent_matrix[kDepotIndex][*path_iter];
    path_iter++;
    while (path_iter != target->path.end())
    {
        path_length += network.adjacent_matrix
                           [*(path_iter - 1)][*(path_iter)];
        path_iter++;
    }
    path_length += network.adjacent_matrix[kDepotIndex][*(path_iter - 1)];
    return path_length;
}

static vector<vector<int>> GenerateSwappingPath(int length)
{
    vector<vector<int>> swapping(length * (length - 1) / 2);
    vector<vector<int>>::iterator swap_iterator = swapping.begin();
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < i; j++)
        {
            (*swap_iterator).resize(2);
            (*swap_iterator).at(0) = i;
            (*swap_iterator).at(1) = j;
            swap_iterator++;
        }
    }
    return swapping;
}

static void TabuSearchPath(CityNetwork network, vector<TspChromosome> &group,
                    vector<TspChromosome>::iterator target, int fitness_mode,
                    vector<vector<int>>::iterator swap_begin,
                    vector<vector<int>>::iterator swap_end)
{
    typedef struct
    {
        vector<int> swapping;
        double fitness;
    } Swapping;

    struct MySort
    {
        bool operator()(Swapping i, Swapping j)
        {
            return (i.fitness < j.fitness);
        }
    } sort_object;
    vector<Swapping> neighborhood;
    vector<vector<int>> tabu_list;
    double current_fitness = 0;
    double best_fitness = 0;
    int num_steps = 2;

    for (int i = 0; i < num_steps; i++)
    {
        for (vector<vector<int>>::iterator j = swap_begin;
             j != swap_end; j++)
        {
            Swapping temp;
            temp.swapping = (*j);
            temp.fitness = current_fitness +
                           TabuSearchFitness(network, group, target, j, fitness_mode);
            neighborhood.push_back(temp);
        }
        sort(neighborhood.begin(), neighborhood.end(), sort_object);
        if (neighborhood.begin()->fitness < best_fitness)
        {
            best_fitness = neighborhood.begin()->fitness;
            current_fitness = neighborhood.begin()->fitness;

            tabu_list.push_back(neighborhood.begin()->swapping);
        }
        else
        {
            for (vector<Swapping>::iterator i = neighborhood.begin();
                 i != neighborhood.end(); i++)
            {
                vector<vector<int>>::iterator temp;
                temp = find(tabu_list.begin(), tabu_list.end(), i->swapping);
                if (temp == tabu_list.end())
                {
                    current_fitness = i->fitness;

                    tabu_list.push_back(i->swapping);
                    break;
                }
            }
        }
        neighborhood.clear();
    }
}

static double TabuSearchFitness(CityNetwork network, vector<TspChromosome> &group,
                         vector<TspChromosome>::iterator target,
                         vector<vector<int>>::iterator swapping,
                         int fitness_mode)
{
    int position_a = (*swapping).at(0);
    int position_b = (*swapping).at(1);
    int path_sequence_length = target->path.size();
    switch (fitness_mode)
    {
    case 0:
    {
        int index_b = (*target).path.at(position_b);
        int index_a = (*target).path.at(position_a);
        int index_a_back;
        int index_a_forward;
        if (position_a == 0)
        {
            index_a_back = kDepotIndex;
            index_a_forward = (*target).path.at(position_a + 1);
        }
        else if (position_a == path_sequence_length - 1)
        {
            index_a_back = (*target).path.at(position_a - 1);
            index_a_forward = kDepotIndex;
        }
        else
        {
            index_a_back = (*target).path.at(position_a - 1);
            index_a_forward = (*target).path.at(position_a + 1);
        }

        int index_b_back;
        int index_b_forward;
        if (position_b == 0)
        {
            index_b_back = kDepotIndex;
            index_b_forward = (*target).path.at(position_b + 1);
        }
        else if (position_b == path_sequence_length - 1)
        {
            index_b_back = (*target).path.at(position_b - 1);
            index_b_forward = kDepotIndex;
        }
        else
        {
            index_b_back = (*target).path.at(position_b - 1);
            index_b_forward = (*target).path.at(position_b + 1);
        }
        double old_len, new_len;
        if (index_a_back == index_b)
        {
            old_len = network.adjacent_matrix[index_b_back][index_b] +
                      network.adjacent_matrix[index_a][index_a_forward];
            new_len = network.adjacent_matrix[index_b_back][index_a] +
                      network.adjacent_matrix[index_b][index_a_forward];
        }
        else if (index_b_back == index_a)
        {
            old_len = network.adjacent_matrix[index_a_back][index_a] +
                      network.adjacent_matrix[index_b][index_b_forward];
            new_len = network.adjacent_matrix[index_a_back][index_b] +
                      network.adjacent_matrix[index_a][index_b_forward];
        }
        else
        {
            old_len = network.adjacent_matrix[index_a_back][index_a] +
                      network.adjacent_matrix[index_a][index_a_forward] +
                      network.adjacent_matrix[index_b_back][index_b] +
                      network.adjacent_matrix[index_b][index_b_forward];
            new_len = network.adjacent_matrix[index_a_back][index_b] +
                      network.adjacent_matrix[index_b][index_a_forward] +
                      network.adjacent_matrix[index_b_back][index_a] +
                      network.adjacent_matrix[index_a][index_b_forward];
        }
        return new_len - old_len;
    }
    case 1:
    {
        int num_chromosome = group.size();
        double delta_hamming_dist = 0;
        double delta_fitness = 0;
        for (vector<TspChromosome>::iterator i = group.begin();
             i != group.end(); i++)
        {
            if (i == target)
                continue;
            int coeff = num_chromosome - (i - group.begin());
            delta_hamming_dist = 0;
            if ((*target).path.at(position_a) == (*i).path.at(position_a))
                delta_hamming_dist += 1;
            else if ((*target).path.at(position_a) == (*i).path.at(position_b))
                delta_hamming_dist -= 1;
            if ((*target).path.at(position_b) == (*i).path.at(position_b))
                delta_hamming_dist += 1;
            else if ((*target).path.at(position_b) == (*i).path.at(position_a))
                delta_hamming_dist -= 1;
            delta_fitness += coeff * delta_hamming_dist;
        }
        return delta_fitness;
    }
    default:
        cout << "Wrong fitness function mode for Tabu Search" << endl;
        return 0;
    }
}