#include <mpi.h>
#include <omp.h>

#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sys/time.h>

#define FILENAME "att48_d.txt"

using namespace std;


int map[100][100];
vector< vector<int> > gene, gene_child;
vector<int> cost, cost_child, child1, child2;
vector<double> fitness, fitness_child;
vector<double> posibility;
int elite1, elite2;
int my_rank, comm_sz;



int cities = 48;
int population = 40;
int generations = 1000;
double crossover_p = 0.8;
double mutation_p = 0.01;
int immigration = 2;
int numThreads = 1;


void test(int argc, char* argv[]);
void init();
void readFile();
void firstGeneration();
void findElite();
bool sameCost(int);
void calculatePosibility();
int whellSelection();
void crossover();
vector<int> CrossoverTour(vector<int> a, vector<int> b);
void mutation();
void immigrate();



int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0)
    {
        cout << "0";
        for (int i = 1; i < 9; i++)
            cout << " \t" << i;
    }
    
    
    for (int i = 0; i < 10; i++)
    {
        
        generations = 300 * (i + 1) / comm_sz;
        if (my_rank == 0)
        {
            cout << endl;
            cout << generations * comm_sz;
        }
        for (int j = 1; j < 9; j++)
        {
            numThreads = j;
            init();
            readFile();
            firstGeneration();
            findElite();
            calculatePosibility();
            struct timeval time0, time1;
            gettimeofday(&time0, NULL);
            test(argc, argv);
            gettimeofday(&time1, NULL);
            if(my_rank == 0)
                cout << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << "\t";
        }
    }
    MPI_Finalize();
    return 0;
}

void test(int argc, char* argv[])
{

    

    srand(time(NULL) + my_rank);


    int i, j;
    for (j = 0; j < generations; j++) {
        for (i = 0; i < (population - immigration) / 2 - 1; i++) {
            child1.clear();
            child2.clear();
            crossover();
            mutation();
        }
        gene_child.push_back(gene[elite1]);
        gene_child.push_back(gene[elite2]);
        cost_child.push_back(cost[elite1]);
        cost_child.push_back(cost[elite2]);
        fitness_child.push_back(fitness[elite1]);
        fitness_child.push_back(fitness[elite2]);

        gene = gene_child;
        cost = cost_child;
        fitness = fitness_child;

        gene_child.clear();
        cost_child.clear();
        fitness_child.clear();
        calculatePosibility();
        
        if(comm_sz != 1)
        {
            for (i = 0; i < immigration; i++) {
            immigrate();
            }
        }
        findElite();

    }


    if (my_rank == 0) {
        for (int i = 1; i < comm_sz; i++) {
            int temp_cost;
            MPI_Recv(&temp_cost, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (temp_cost < cost[elite1])
                cost[elite1] = temp_cost;
        }
        
    }
    else {
        MPI_Send(&cost[elite1], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    
}

void init()
{
    gene.clear();
    gene_child.clear();
    cost.clear();
    cost_child.clear();
    child1.clear();
    child2.clear();
    fitness.clear();
    fitness_child.clear();
    posibility.clear();
}


void readFile() {
    ifstream f(FILENAME);
    for (int i = 0; i < cities; i++) {
        for (int j = 0; j < cities; j++) {
            f >> map[i][j];
        }
    }
    f.close();
}


void firstGeneration() {
    vector<int> randomnum;
    for (int i = 0; i < cities; i++)
        randomnum.push_back(i);
    for (int i = 0; i < population; i++) {
        random_shuffle(randomnum.begin(), randomnum.end());
        gene.push_back(randomnum);
        int cost_sum = 0;
        for (int j = 0; j < cities - 1; j++)
            cost_sum += map[randomnum[j]][randomnum[j + 1]];
        cost_sum += map[randomnum[0]][randomnum[cities - 1]];
        cost.push_back(cost_sum);
        fitness.push_back((double)1 / (double)cost_sum);
    }
}


void findElite() {
    int min = 2147483647, sec = 2147483647;
    elite1 = 0, elite2 = 0;

#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < cost.size(); i++) {
        if (cost[i] < min) {
            elite1 = i;
            min = cost[i];
        }
    }

#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < cost.size(); i++) {
        if (cost[i] < sec && cost[i] > min) {
            elite2 = i;
            sec = cost[i];
        }
    }
}

void calculatePosibility() {
    posibility.clear();
    double fitness_sum = 0.0;
    for (int i = 0; i < fitness.size(); i++)
        fitness_sum += fitness[i];
    for (int i = 0; i < fitness.size(); i++)
        posibility.push_back(fitness[i] / fitness_sum);
}


int whellSelection() {
    double ran = (double)rand() / RAND_MAX;
    int i = 0;
    while (ran > 0) {
        for (i = 0; i < posibility.size(); i++) {
            ran -= posibility[i];
            if (ran < 0) return i;
        }
    }
    return i;
}


void crossover() {
    int dad = whellSelection();
    int mom = whellSelection();
    while (dad == mom)
        mom = whellSelection();

    double ran = (double)rand() / RAND_MAX;

    if (ran < crossover_p) {
        int start = rand() % (cities - 1);
        int end = rand() % cities;
        while (start >= end)
            end = rand() % cities;

        vector<int> dad_buf, mom_buf;
        dad_buf = gene[dad];
        mom_buf = gene[mom];

        child1 = CrossoverTour(dad_buf, mom_buf);
        child2 = CrossoverTour(mom_buf, dad_buf);

    }
    else {
        child1 = gene[dad];
        child2 = gene[mom];
    }


    int cost_sum = 0;

    for (int j = 0; j < cities - 1; j++)
        cost_sum += map[child1[j]][child1[j + 1]];
    cost_sum += map[child1[0]][child1[cities - 1]];
    cost_child.push_back(cost_sum);
    fitness_child.push_back((double)1 / (double)cost_sum);

    cost_sum = 0;
    for (int j = 0; j < cities - 1; j++)
        cost_sum += map[child2[j]][child2[j + 1]];
    cost_sum += map[child2[0]][child2[cities - 1]];
    cost_child.push_back(cost_sum);
    fitness_child.push_back((double)1 / (double)cost_sum);

    gene_child.push_back(child1);
    gene_child.push_back(child2);
}


vector<int> CrossoverTour(vector<int> a, vector<int> b) {
    vector<int> result(cities);

    bool contained[100] = { 0 };
    int start_index = rand() % cities;
    int end_index = rand() % cities;

    while (start_index == end_index)
        end_index = rand() % cities;

    if (start_index > end_index)
        swap(start_index, end_index);

    for (int i = start_index; i <= end_index; ++i)
        contained[result[i] = a[i]] = true;

    for (int i = 0; i < cities; ++i) {
        if (!contained[b[i]]) {
            for (int k = 0; k < cities; ++k) {
                if (!result[k]) {
                    contained[result[k] = b[i]] = true;
                    break;
                }
            }
        }
    }

    return result;
}


void mutation() {
    double ran = (double)rand() / RAND_MAX;
    if (ran < mutation_p) {
        int index1, index2, tmp;
        index1 = rand() % cities;
        index2 = rand() % cities;
        tmp = child1[index1];
        child1[index1] = child1[index2];
        child1[index2] = tmp;
        tmp = child2[index1];
        child2[index1] = child2[index2];
        child2[index2] = tmp;
    }
}


void immigrate() {
    int index = whellSelection();
    vector<int> recv, send;
    send = gene[index];
    recv.resize(send.size());

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank % 2) {
        MPI_Send(&send[0], cities, MPI_INT, (my_rank + 1) % comm_sz, 0, MPI_COMM_WORLD);
        MPI_Recv(&recv[0], cities, MPI_INT, (my_rank - 1 + comm_sz) % comm_sz, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(&recv[0], cities, MPI_INT, (my_rank - 1 + comm_sz) % comm_sz, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&send[0], cities, MPI_INT, (my_rank + 1) % comm_sz, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    gene_child.push_back(recv);
    int cost_sum = 0;
    for (int j = 0; j < cities - 1; j++)
        cost_sum += map[recv[j]][recv[j + 1]];
    cost_sum += map[recv[0]][recv[cities - 1]];
    cost_child.push_back(cost_sum);
    fitness_child.push_back((double)1 / (double)cost_sum);
}
