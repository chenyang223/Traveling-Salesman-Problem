#include "ga.h"

typedef struct
{
    double fitness_func;
    vector<int> path;
    vector<int> salesman;
} MtspChromosome;

const int kNumInGroup = 2;
const int kNumSalesman = 2;
const double kDoubleLargeNumber = DBL_MAX;

static void GreedyAlgorithm(CityNetwork network, vector<MtspChromosome> &group);
static void InitializeGroup(CityNetwork network, vector<MtspChromosome> &group);
static vector<int>::iterator FindGreedy(CityNetwork network,
                                        vector<int>::iterator begin,
                                        vector<int>::iterator end, int start_point);
static void PrintMtspChromosomeVector(vector<MtspChromosome> v);
static void MOGAMain(CityNetwork network, vector<MtspChromosome> &group);
static void TabuSearchPath(CityNetwork network, vector<MtspChromosome> &group,
                           vector<MtspChromosome>::iterator target, int fitness_mode,
                           vector<vector<int>>::iterator swap_begin,
                           vector<vector<int>>::iterator swap_end);
static void TabuSearchSalesman(CityNetwork network, vector<MtspChromosome> &group,
                               vector<MtspChromosome>::iterator target,
                               vector<vector<int>>::iterator swap_begin,
                               vector<vector<int>>::iterator swap_end);
static double TabuSearchFitness(CityNetwork network, vector<MtspChromosome> &group,
                                vector<MtspChromosome>::iterator target,
                                vector<vector<int>>::iterator swapping,
                                int fitness_mode);
static vector<vector<int>> GenerateSwappingPath(int length);
static vector<vector<int>> GenerateSwappingSalesman(int length);
static double CalculatePathFitness(CityNetwork network,
                                   vector<MtspChromosome>::iterator target);
static void OutputBestSolution(vector<MtspChromosome>::iterator best);

void MtspMOGA()
{
    CityNetwork network;
    network = ConstructNetwork();
    //PrintNetwork(network);

    vector<MtspChromosome> group(kNumInGroup);
    InitializeGroup(network, group);
    //PrintMtspChromosomeVector(group);

    MOGAMain(network, group);
    PrintMtspChromosomeVector(group);
    OutputBestSolution(group.begin());
}

static void OutputBestSolution(vector<MtspChromosome>::iterator best)
{
    stringstream filename_stream;
    string output_filename;
    filename_stream << GetDataPath() << "mtsp_path";
    filename_stream << "_" << GetDataFile();
    filename_stream >> output_filename;
    ofstream output;
    output.open(output_filename);

    int depot_index = best->path.size();
    output << "NUM_SALESMAN: " << best->salesman.size() << "\n";
    output << "PATH_LENGTH: " << best->fitness_func << "\n";
    vector<int>::iterator salesman_iter = best->salesman.begin();
    vector<int>::iterator path_iter = best->path.begin();

    for (; salesman_iter != best->salesman.end(); salesman_iter++)
    {
        int num = *salesman_iter;
        output << depot_index << " ";
        for (int i = 0; i < num; i++, path_iter++)
            output << *path_iter << " ";
        output << depot_index << "\n";
    }
    output.close();
}

static void InitializeGroup(CityNetwork network, vector<MtspChromosome> &group)
{
    int path_sequence_length = network.num_nodes - 1;
    for (vector<MtspChromosome>::iterator i = group.begin();
         i != group.end(); i++)
    {
        i->salesman.resize(kNumSalesman, 0);
        for (int j = 0; j < path_sequence_length; j++)
        {
            int temp = rand() % kNumSalesman;
            i->salesman.at(temp)++;
        }
        vector<int>::iterator min_iterator =
            min_element(i->salesman.begin(), i->salesman.end());
        while (*min_iterator == 0)
        {
            vector<int>::iterator max_iterator =
                max_element(i->salesman.begin(), i->salesman.end());
            (*min_iterator)++;
            (*max_iterator)--;
            min_iterator = min_element(i->salesman.begin(), i->salesman.end());
        }
    }
    GreedyAlgorithm(network, group);
}

static void GreedyAlgorithm(CityNetwork network, vector<MtspChromosome> &group)
{
    int depot_index = network.num_nodes - 1;
    int path_sequence_length = network.num_nodes - 1;

    for (vector<MtspChromosome>::iterator i = group.begin();
         i != group.end(); i++)
    {
        vector<int> base_path;
        for (int j = 0; j < path_sequence_length; j++)
            base_path.push_back(j);
        double fitness = 0;
        vector<int>::iterator salesman_iterator = i->salesman.begin();
        int counter = *salesman_iterator;
        vector<int>::iterator path_iterator = base_path.begin();
        vector<int>::iterator greedy_iterator;
        int start_point = rand() % path_sequence_length;
        swap((*path_iterator), base_path.at(start_point));
        fitness += network.adjacent_matrix[*path_iterator][depot_index];
        counter--;
        path_iterator++;
        while (path_iterator != base_path.end())
        {
            if (counter == 0)
            {
                fitness += network.adjacent_matrix
                               [*(path_iterator - 1)][depot_index];
                salesman_iterator++;
                counter = *salesman_iterator;
                greedy_iterator = FindGreedy(
                    network, path_iterator, base_path.end(), depot_index);
                swap(*greedy_iterator, *path_iterator);
                fitness += network.adjacent_matrix
                               [*path_iterator][depot_index];
                path_iterator++;
            }
            else
            {
                greedy_iterator = FindGreedy(
                    network, path_iterator, base_path.end(),
                    *(path_iterator - 1));
                swap(*greedy_iterator, *path_iterator);
                fitness += network.adjacent_matrix
                               [*(path_iterator - 1)][*path_iterator];
                path_iterator++;
            }
            counter--;
        }
        fitness += network.adjacent_matrix
                       [*(path_iterator - 1)][depot_index];
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

static void PrintMtspChromosomeVector(vector<MtspChromosome> v)
{
    cout << "num of members: " << v.size() << "\n";
    for (int i = 0; i < v.size(); i++)
    {
        cout << i << "  fitness: " << v.at(i).fitness_func << "\n";
        cout << "path: "
             << "\n";
        vector<int>::iterator p = v.at(i).path.begin();
        for (vector<int>::iterator s = v.at(i).salesman.begin();
             s != v.at(i).salesman.end(); s++)
        {
            cout << "num of cities: " << *s << " : ";
            for (int counter = *s; counter != 0; counter--, p++)
                cout << *p << "->";
            cout << "\n";
        }
    }
}

static void MOGAMain(CityNetwork network, vector<MtspChromosome> &group)
{
    struct MySort
    {
        bool operator()(MtspChromosome i, MtspChromosome j)
        {
            return (i.fitness_func < j.fitness_func);
        }
    } sort_object;

    vector<vector<int>> swapping_path =
        GenerateSwappingPath(network.num_nodes - 1);
    vector<vector<int>> swapping_salesman =
        GenerateSwappingSalesman(kNumSalesman);

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
            TabuSearchSalesman(network, group, group.begin() + j,
                               swapping_salesman.begin(),
                               swapping_salesman.end());
            TabuSearchPath(network, group, group.begin() + j,
                           fitness_mode, swapping_path.begin(),
                           swapping_path.end());
            (group.begin() + j)->fitness_func =
                CalculatePathFitness(network, group.begin() + j);
        }
    }
}

static double CalculatePathFitness(CityNetwork network,
                            vector<MtspChromosome>::iterator target)
{
    int depot_index = network.num_nodes - 1;
    double path_length = 0;
    vector<int>::iterator path_iter = target->path.begin();
    vector<int> salesman_copy = target->salesman;
    vector<int>::iterator salesman_iter = salesman_copy.begin();
    path_length += network.adjacent_matrix[depot_index][*path_iter];
    path_iter++;
    (*salesman_iter)--;
    while (path_iter != target->path.end())
    {
        if (*salesman_iter == 0)
        {
            path_length += network.adjacent_matrix
                               [depot_index][*(path_iter - 1)];
            path_length += network.adjacent_matrix
                               [depot_index][*(path_iter)];
            salesman_iter++;
        }
        else
        {
            path_length += network.adjacent_matrix
                               [*(path_iter - 1)][*(path_iter)];
        }
        path_iter++;
        (*salesman_iter)--;
    }
    path_length += network.adjacent_matrix[depot_index][*(path_iter - 1)];
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

static vector<vector<int>> GenerateSwappingSalesman(int length)
{
    vector<vector<int>> swapping(length * (length - 1));
    vector<vector<int>>::iterator swap_iterator = swapping.begin();
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            if (i == j)
                continue;

            (*swap_iterator).resize(2);
            (*swap_iterator).at(0) = i;
            (*swap_iterator).at(1) = j;
            swap_iterator++;
        }
    }
    return swapping;
}

static void TabuSearchPath(CityNetwork network, vector<MtspChromosome> &group,
                    vector<MtspChromosome>::iterator target, int fitness_mode,
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
    int num_steps = 3;

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

static void TabuSearchSalesman(CityNetwork network, vector<MtspChromosome> &group,
                        vector<MtspChromosome>::iterator target,
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
    double best_fitness = kDoubleLargeNumber;
    int num_steps = 2;

    for (int i = 0; i < num_steps; i++)
    {
        for (vector<vector<int>>::iterator j = swap_begin;
             j != swap_end; j++)
        {
            if (target->salesman.at((*j).at(0)) == 1)
                continue;
            Swapping temp;
            temp.swapping = (*j);
            (target->salesman.at((*j).at(0)))--;
            (target->salesman.at((*j).at(1)))++;
            temp.fitness = CalculatePathFitness(network, target);
            neighborhood.push_back(temp);
            (target->salesman.at((*j).at(0)))++;
            (target->salesman.at((*j).at(1)))--;
        }
        sort(neighborhood.begin(), neighborhood.end(), sort_object);
        if (neighborhood.begin()->fitness < best_fitness)
        {
            best_fitness = neighborhood.begin()->fitness;
            current_fitness = neighborhood.begin()->fitness;
            (target->salesman.at((neighborhood.begin()->swapping).at(0)))--;
            (target->salesman.at((neighborhood.begin()->swapping).at(1)))++;
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
                    (target->salesman.at((neighborhood.begin()->swapping).at(0)))--;
                    (target->salesman.at((neighborhood.begin()->swapping).at(1)))++;
                    tabu_list.push_back(i->swapping);
                    break;
                }
            }
        }
        neighborhood.clear();
    }
}

static double TabuSearchFitness(CityNetwork network, vector<MtspChromosome> &group,
                         vector<MtspChromosome>::iterator target,
                         vector<vector<int>>::iterator swapping,
                         int fitness_mode)
{
    int position_a = (*swapping).at(0);
    int position_b = (*swapping).at(1);
    int depot_index = network.num_nodes - 1;
    switch (fitness_mode)
    {
    case 0:
    {
        int index_b = (*target).path.at(position_b);
        int index_a = (*target).path.at(position_a);
        vector<int>::iterator salesman = (target->salesman).begin();
        int position_a_copy = position_a;
        int index_a_back;
        int index_a_forward;
        while (position_a_copy >= *salesman)
        {
            position_a_copy -= *salesman;
            salesman++;
        }
        if (position_a_copy == 0)
        {
            index_a_back = depot_index;
            if (*salesman == 1)
                index_a_forward = depot_index;
            else
                index_a_forward = (*target).path.at(position_a + 1);
        }
        else if (position_a_copy == *salesman - 1)
        {
            index_a_back = (*target).path.at(position_a - 1);
            index_a_forward = depot_index;
        }
        else
        {
            index_a_back = (*target).path.at(position_a - 1);
            index_a_forward = (*target).path.at(position_a + 1);
        }

        salesman = (target->salesman).begin();
        int position_b_copy = position_b;
        int index_b_back;
        int index_b_forward;
        while (position_b_copy >= *salesman)
        {
            position_b_copy -= *salesman;
            salesman++;
        }
        if (position_b_copy == 0)
        {
            index_b_back = depot_index;
            if (*salesman == 1)
                index_b_forward = depot_index;
            else
                index_b_forward = (*target).path.at(position_b + 1);
        }
        else if (position_b_copy == *salesman - 1)
        {
            index_b_back = (*target).path.at(position_b - 1);
            index_b_forward = depot_index;
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
        for (vector<MtspChromosome>::iterator i = group.begin();
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
