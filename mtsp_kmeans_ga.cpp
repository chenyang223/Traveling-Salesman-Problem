//no depot mtsp

#include "mtsp_kmeans_ga.h"
#include "utils.h"
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

void MtspKmeansGA()
{
    CityNetwork network;
    network = ConstructNetwork();
//    PrintNetwork(network);
    vector<vector<int>> kmeans_cluster;
    kmeans_cluster = ReadKMeansCluster();
    int num_salesman = kmeans_cluster.size();
    GetNumInGroup() = kBaseNumInGroup * num_salesman;
    GetGASteps() = kBaseGaSteps * num_salesman;
    vector<vector<int>> path_result;
    vector<vector<int>> init_path_result;
    double total_path_length = 0.0;
    double init_total_path_length = 0.0;

    time_t time1 = clock();
    time_t time2 = clock();
    double init_time = 0;
    double ga_time = 0;

    for (auto i = kmeans_cluster.begin();
         i != kmeans_cluster.end(); i++)
    {
        cout << "current salesman : " << i - kmeans_cluster.begin() + 1
             << "/" << num_salesman << endl;
        vector<TspChromosome> group(GetNumInGroup());
        InitializeGroup(network, group, i);
        init_path_result.push_back(group.begin()->path);
        init_total_path_length += group.begin()->length;
        for (auto j = group.begin();
             j != group.end(); j++)
        {
            LocalSearchPath(network, j->path);
            j->length = CalculatePathLength(network, j->path);
        }

        time1 = clock();
        init_time += double(time1 - time2) / CLOCKS_PER_SEC;
        cout << "initializing time : " << double(time1 - time2) / CLOCKS_PER_SEC << endl;
        GAMain(network, group);
        LocalSearchPath(network, group.begin()->path);
        group.begin()->length = CalculatePathLength(network, group.begin()->path);
        time2 = clock();
        ga_time += double(time2 - time1) / CLOCKS_PER_SEC;
        path_result.push_back(group.begin()->path);
        total_path_length += group.begin()->length;
    }

    cout << "total initializing time : " << init_time << "\n";
    cout << "total ga time : " << ga_time << "\n";

    PrintMtspResult(path_result, total_path_length);
    OutputBestSolution(path_result, total_path_length, num_salesman);
    OutputInitSolution(init_path_result, init_total_path_length, num_salesman);
//    OutputBestInGeneration(num_salesman);
}

void GAMain(CityNetwork &network, vector<TspChromosome> &group)
{
    struct MySort
    {
        bool operator()(TspChromosome i, TspChromosome j)
        {
            return (i.length < j.length);
        }
    } sort_object;
    sort(group.begin(), group.end(), sort_object);

    time_t time1 = clock();
    time_t time2 = clock();
    double crossover_time = 0;
    double local_search_time = 0;

    for (int i = 0; i < GetGASteps(); i++)
    {
        // cout << "current generation : " << i << "\n";
        int father_index = rand() % GetNumInGroup();
        int mother_index;
        do
        {
            mother_index = rand() % GetNumInGroup();
        } while (mother_index == father_index);
        TspChromosome son;
        son.path = DPXCrossover(network, group.at(father_index).path,
                                group.at(mother_index).path);

        time1 = clock();
        crossover_time += double(time1 - time2) / CLOCKS_PER_SEC;
        // cout << "crossover time : " << double(time1 - time2) / CLOCKS_PER_SEC << "\n";
        LocalSearchPath(network, son.path);
        time2 = clock();
        son.length = CalculatePathLength(network, son.path);

        double old_best_path_length = group.begin()->length;

        int insert_or_not = 0;
        for (auto j = group.begin(); j != group.end(); j++)
        {
            if (son.length > j->length)
                continue;
            group.insert(j, son);
            insert_or_not = 1;
            break;
        }
        if (insert_or_not == 1)
            group.pop_back();
            
        double new_best_path_length = group.begin()->length;
        if (new_best_path_length < old_best_path_length)
        {
            LocalSearchPath(network, group.begin()->path);
            group.begin()->length = CalculatePathLength(network, group.begin()->path);
        }
        GetBestAtOneGeneration().at(i) += group.begin()->length;
        local_search_time += double(time2 - time1) / CLOCKS_PER_SEC;
        // cout << "local search time : " << double(time2 - time1) / CLOCKS_PER_SEC << "\n";
        // cout << "current father solution : " << group.at(father_index).length << "\n";
        // cout << "current mother solution : " << group.at(mother_index).length << "\n";
        // cout << "current son solution : " << son.length << "\n";
        // cout << "current best solution : " << group.begin()->length << "\n"
        //      << endl;
        if (group.begin()->length == (group.end() - 1)->length)
            break;
    }

    cout << "crossover time : " << crossover_time << "\n";
    cout << "local search time : " << local_search_time << "\n";
}

vector<vector<int>> ReadKMeansCluster()
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

vector<double> &GetBestAtOneGeneration()
{
    static vector<double> length(GetGASteps(), 0);
    return length;
}

vector<int> &GetLKImprovedPath()
{
    static vector<int> improved_path;
    return improved_path;
}

int &GetNumInGroup()
{
    static int num_in_group;
    return num_in_group;
}

int &GetGASteps()
{
    static int ga_steps;
    return ga_steps;
}

void PrintMtspResult(vector<vector<int>> v, double length)
{
    cout << "num of salesman: " << v.size() << "\n";
    cout << "path length: " << length << "\n";
    for (auto i = v.begin(); i != v.end(); i++)
    {
        for (auto j = (*i).begin(); j != (*i).end(); j++)
            cout << *j << " ";
        cout << "\n";
    }
    cout << endl;
}

void OutputBestInGeneration(int num_salesmans)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "length-generation_";
    filename_stream << num_salesmans;
    filename_stream << "_" << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    for (auto i = GetBestAtOneGeneration().begin();
         i != GetBestAtOneGeneration().end(); i++)
    {
        if (*i == 0)
            break;
        output << *i << " ";
    }
    output.close();
}

void OutputBestSolution(vector<vector<int>> v, double length, int num_salesman)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "mtsp_path_kmeansGA";
    filename_stream << "_" << num_salesman;
    filename_stream << "_" << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    output << "NUM_SALESMAN: " << v.size() << "\n";
    output << "PATH_LENGTH: " << length << "\n";
    for (auto i = v.begin(); i != v.end(); i++)
    {
        for (auto j = (*i).begin(); j != (*i).end(); j++)
            output << *j << " ";
        output << *((*i).begin()) << " ";
        output << "\n";
    }
    output.close();
}

void OutputInitSolution(vector<vector<int>> v, double length, int num_salesman)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "mtsp_init_path_kmeansGA";
    filename_stream << "_" << num_salesman;
    filename_stream << "_" << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    output << "NUM_SALESMAN: " << v.size() << "\n";
    output << "PATH_LENGTH: " << length << "\n";
    for (auto i = v.begin(); i != v.end(); i++)
    {
        for (auto j = (*i).begin(); j != (*i).end(); j++)
            output << *j << " ";
        output << *((*i).begin()) << " ";
        output << "\n";
    }
    output.close();
}

void InitializeGroup(CityNetwork &network, vector<TspChromosome> &group,
                            vector<vector<int>>::iterator base_path_iter)
{
    int path_sequence_length = (*base_path_iter).size();
    for (auto i = group.begin();
         i != group.end(); i++)
    {
        vector<int> base_path = *base_path_iter;
		if (path_sequence_length < 4)
		{
            i->path.assign(base_path.begin(), base_path.end());
            i->length = CalculatePathLength(network, base_path);
			continue;
		}
        double length = 0;
        int start_point;
        auto path_iterator = base_path.begin();
        vector<int>::iterator greedy_iterator;

        start_point = rand() % path_sequence_length;
        swap((*path_iterator), base_path.at(start_point));
        path_iterator++;
        while (path_iterator != base_path.end())
        {
            greedy_iterator = FindGreedy(network, path_iterator, base_path.end(),
                                         *(path_iterator - 1));
            swap(*greedy_iterator, *path_iterator);
            length += network.adjacent_matrix
                           [*(path_iterator - 1)][*path_iterator];
            path_iterator++;
        }
        length += network.adjacent_matrix
                       [*(path_iterator - 1)][start_point];
        i->path.assign(base_path.begin(), base_path.end());
        reverse(i->path.begin(), i->path.end());
        i->length = length;
    }
}

vector<int>::iterator FindGreedy(CityNetwork &network,
                                        vector<int>::iterator begin,
                                        vector<int>::iterator end, int start_point)
{
    auto greedy_iterator = begin;
    double min_path_length = network.adjacent_matrix[start_point][*begin];
    for (auto i = begin; i != end; i++)
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

double CalculatePathLength(CityNetwork &network, vector<int> path)
{
    double path_length = 0;
	if (path.size() == 1)
		return path_length;
    path_length += network.adjacent_matrix[*(path.begin())]
                                          [*(path.end() - 1)];
    for (auto i = path.begin(); i != path.end() - 1; i++)
        path_length += network.adjacent_matrix[*i][*(i + 1)];
    return path_length;
}

vector<int> DPXCrossover(CityNetwork &network, vector<int> father, vector<int> mother)
{
    vector<int> son;
	if (father.size() < 4)
	{
		son = father;
		return son;
	}

    vector<vector<int>> fragments;

    // Build two vectors to save the back forward relation of mother
    int link_vector_size = *max_element(mother.begin(), mother.end()) + 1;
    vector<int> mother_back_vector(link_vector_size);
    vector<int> mother_forward_vector(link_vector_size);
    for (auto i = mother.begin(); i != mother.end(); i++)
    {
        if (i == mother.end() - 1)
        {
            mother_forward_vector.at(*i) = *(mother.begin());
            mother_back_vector.at(*i) = *(i - 1);
        }
        else if (i == mother.begin())
        {
            mother_forward_vector.at(*i) = *(i + 1);
            mother_back_vector.at(*i) = *(mother.end() - 1);
        }
        else
        {
            mother_forward_vector.at(*i) = *(i + 1);
            mother_back_vector.at(*i) = *(i - 1);
        }
    }

    //Find all cut points in father
    vector<int> cut_point;
    for (auto i = father.begin(); i != father.end(); i++)
    {
        vector<int>::iterator i_forward;
        if (i == father.end() - 1)
            i_forward = father.begin();
        else
            i_forward = i + 1;

        if (*i_forward != mother_forward_vector.at(*i) &&
            *i_forward != mother_back_vector.at(*i))
            cut_point.push_back(i - father.begin());
    }

    //Find the fragments after cutting
    int num_fragments = cut_point.size();
    if (num_fragments == 0)
        return father;
    vector<int> local_fragment;
    //head part
    local_fragment.assign(father.begin(),
                          father.begin() + *(cut_point.begin()) + 1);
    //tail part
    local_fragment.insert(local_fragment.begin(),
                          father.begin() + 1 + cut_point.at(num_fragments - 1),
                          father.end());
    fragments.push_back(local_fragment);
    for (auto i = cut_point.begin() + 1;
         i != cut_point.end(); i++)
    {
        local_fragment.assign(father.begin() + *(i - 1) + 1,
                              father.begin() + *(i) + 1);
        fragments.push_back(local_fragment);
    }

    // Collect the start and end point of all fragments
    vector<int> start_point_vector;
    vector<int> end_point_vector;
    for (auto i = fragments.begin();
         i != fragments.end(); i++)
    {
        start_point_vector.push_back(*((*i).begin()));
        end_point_vector.push_back(*((*i).end() - 1));
    }

    // Rearange all fragments follow greedy algorithms
    int start_fragment_index = rand() % num_fragments;
    son = fragments.at(start_fragment_index);
    int back_end_point = *(son.end() - 1);
    start_point_vector.at(start_fragment_index) = -1;
    end_point_vector.at(start_fragment_index) = -1;
    for (int i = 0; i < num_fragments - 1; i++)
    {
        int best_start;
        double best_length_start = kDoubleLargeNumber;
        for (auto j = start_point_vector.begin();
             j != start_point_vector.end(); j++)
        {
            if (*j == -1)
                continue;
            double length = network.adjacent_matrix[back_end_point][*j];
            if (length < best_length_start)
            {
                best_start = *j;
                best_length_start = length;
            }
        }
        int best_end;
        double best_length_end = kDoubleLargeNumber;
        for (auto j = end_point_vector.begin();
             j != end_point_vector.end(); j++)
        {
            if (*j == -1)
                continue;
            double length = network.adjacent_matrix[back_end_point][*j];
            if (length < best_length_end)
            {
                best_end = *j;
                best_length_end = length;
            }
        }
        if (best_length_start < best_length_end)
        {
            auto iter = find(start_point_vector.begin(),
                                              start_point_vector.end(), best_start);
            int index = iter - start_point_vector.begin();
            start_point_vector.at(index) = -1;
            end_point_vector.at(index) = -1;
            son.insert(son.end(), (fragments.at(index)).begin(),
                       (fragments.at(index)).end());
        }
        else
        {
            auto iter = find(end_point_vector.begin(),
                                              end_point_vector.end(), best_end);
            int index = iter - end_point_vector.begin();
            start_point_vector.at(index) = -1;
            end_point_vector.at(index) = -1;
            son.insert(son.end(), (fragments.at(index)).rbegin(),
                       (fragments.at(index)).rend());
        }
        back_end_point = *(son.end() - 1);
    }
    reverse(son.begin(), son.end());
    return son;
}

// Recursion LK from one paper, not tested
bool LKImprovePath(CityNetwork &network, vector<int> path, int depth)
{
    if (path.size() < 4)
        return false;
    int depth_limit = 2;
    if (depth < depth_limit)
    {
        for (auto i = path.begin() + 2; i != path.end() - 1; i++)
        {
            double delta_length = network.adjacent_matrix[*(i - 1)][*(path.end() - 1)] -
                                  network.adjacent_matrix[*(i - 1)][*i];
            if (delta_length < 0)
            {
                double delta_circle = network.adjacent_matrix[*(path.begin())][*(i)] -
                                      network.adjacent_matrix[*(path.begin())][*(path.end() - 1)];
                reverse(i, path.end());
                if (delta_circle + delta_length < 0)
                {
                    GetLKImprovedPath() = path;
                    return true;
                }
                else if (LKImprovePath(network, path, depth + 1))
                    return true;
                else
                    return false;
            }
        }
        return false;
    }
    else
    {
        auto best_iter = path.begin() + 1;
        double best_delta_l = network.adjacent_matrix[*(best_iter - 1)][*(path.end() - 1)] -
                              network.adjacent_matrix[*(best_iter - 1)][*best_iter];
        for (auto i = path.begin() + 1; i != path.end() - 1; i++)
        {
            double current_delta_l = network.adjacent_matrix[*(i - 1)][*(path.end() - 1)] -
                                     network.adjacent_matrix[*(i - 1)][*i];
            if (current_delta_l < best_delta_l)
            {
                best_iter = i;
                best_delta_l = current_delta_l;
            }
        }
        if (best_delta_l < 0)
        {
            double delta_circle = network.adjacent_matrix[*(path.begin())][*(best_iter)] -
                                  network.adjacent_matrix[*(path.begin())][*(path.end() - 1)];
            reverse(best_iter, path.end());
            if (delta_circle + best_delta_l < 0)
            {
                GetLKImprovedPath() = path;
                return true;
            }
            else if (LKImprovePath(network, path, depth + 1))
                return true;
            else
                return false;
        }
        else
            return false;
    }
}

void LocalSearchPath(CityNetwork &network, vector<int> &path)
{
    if (path.size() < 4)
        return;
    for (auto j = path.begin() + 1; j != path.end() - 2; j++)
    {
        double delta_l = network.adjacent_matrix[*j][*(path.end() - 1)] +
                         network.adjacent_matrix[*(path.begin())][*(j + 1)] -
                         network.adjacent_matrix[*j][*(j + 1)] -
                         network.adjacent_matrix[*(path.begin())][*(path.end() - 1)];
        if (delta_l < 0)
            reverse(path.begin(), j + 1);
    }
    for (auto i = path.begin() + 1; i != path.end() - 2; i++)
    {
        for (auto j = i + 1; j != path.end() - 1; j++)
        {
            double delta_l = network.adjacent_matrix[*(i - 1)][*j] +
                             network.adjacent_matrix[*i][*(j + 1)] -
                             network.adjacent_matrix[*(i - 1)][*i] -
                             network.adjacent_matrix[*j][*(j + 1)];
            if (delta_l < 0)
                reverse(i, j + 1);
        }
    }
}