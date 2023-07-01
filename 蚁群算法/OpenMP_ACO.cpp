#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <Windows.h>
#include "ants.h"
#include <omp.h>

#define NUM_THREAD 8

using namespace std;

struct Ant
{
    float tourLength;
    bool tabu[MAX_NODES];
    int path[MAX_NODES];
    int curNode, nextNode, pathIndex;
};

Ant antColony[MAX_ANTS];
double pheromone[MAX_NODES][MAX_NODES];
float best = (float)MAX_TOUR;
int bestIndex;
EdgeMatrix *dist;

void init()
{
    for (int from = 0; from < MAX_NODES; from++)
    {
        for (int to = 0; to < MAX_NODES; to++)
        {
            pheromone[from][to] = INIT_PHER;
        }
    }

    for (int ant = 0; ant < MAX_ANTS; ant++)
    {
        antColony[ant].curNode = rand() % MAX_NODES;
        for (int to = 0; to < MAX_NODES; to++)
        {
            antColony[ant].tabu[to] = 0;
        }
        antColony[ant].pathIndex = 1;
        antColony[ant].path[0] = antColony[ant].curNode;
        antColony[ant].nextNode = -1;
        antColony[ant].tourLength = 0;
        antColony[ant].tabu[antColony[ant].curNode] = 1;
    }
}

void restartAnts()
{
#pragma omp parallel for num_threads(NUM_THREAD)
    for (int ant = 0; ant < MAX_ANTS; ant++)
    {
        antColony[ant].nextNode = -1;
        antColony[ant].tourLength = 0.0;
        for (int i = 0; i < MAX_NODES; i++)
        {
            antColony[ant].tabu[i] = 0;
            antColony[ant].path[i] = -1;
        }
        antColony[ant].curNode = rand() % MAX_NODES;
        antColony[ant].pathIndex = 1;
        antColony[ant].path[0] = antColony[ant].curNode;
        antColony[ant].tabu[antColony[ant].curNode] = 1;
    }
}
long double calculateProbability(int from, int to)
{
    long double a = pow(pheromone[from][to], ALPHA) * pow((1.0 / (*dist)[from][to]), BETA);
    return (long double)(pow(pheromone[from][to], ALPHA) * pow((1.0 / (*dist)[from][to]), BETA));
}

int selectNextCity(int ant)
{
    int from = antColony[ant].curNode;
    long double sum = 0.0;
    for (int to = 0; to < MAX_NODES; to++)
    {
        if (antColony[ant].tabu[to] == 0)
        {
            sum += calculateProbability(from, to);
        }
    }

    int lastBestIndex = 0.0;
    srand((unsigned)time(NULL));
    long double luckyNumber = (long double)rand() / RAND_MAX;

    for (int to = 0; to < MAX_NODES; to++)
    {
        if (antColony[ant].tabu[to] == 0)
        {
            long double product = calculateProbability(from, to) / sum;
            if (product > 0)
            {
                luckyNumber -= product;
                lastBestIndex = to;

                if (luckyNumber <= 0.0)
                {
                    return to;
                }
            }
        }
    }
    return lastBestIndex;
}

void updatePheromone()
{
    int from, to, i, ant;
#pragma omp parallel for num_threads(NUM_THREAD)
    for (from = 0; from < MAX_NODES; from++)
    {
        for (to = 0; to < from; to++)
        {
            pheromone[from][to] *= 1.0 - RHO;

            if (pheromone[from][to] < 0.0)
            {
                pheromone[from][to] = INIT_PHER;
            }
            pheromone[to][from] = pheromone[from][to];
        }
    }

#pragma omp parallel for num_threads(NUM_THREAD)
    for (ant = 0; ant < MAX_ANTS; ant++)
    {
        for (i = 0; i < MAX_NODES; i++)
        {
            from = antColony[ant].path[i];
            if (i < MAX_NODES - 1)
            {
                to = antColony[ant].path[i + 1];
            }
            else
            {
                to = antColony[ant].path[0];
            }
            pheromone[from][to] += (QVAL / antColony[ant].tourLength);
            pheromone[to][from] = pheromone[from][to];
        }
    }
}

float euclideanDistance(int x1, int x2, int y1, int y2)
{
    int xd = x1 - x2;
    int yd = y1 - y2;
    return (float)(sqrt(xd * xd + yd * yd) + 0.001);
}

float pseudoEuclideanDistance(int x1, int x2, int y1, int y2)
{
    int xd = x1 - x2;
    int yd = y1 - y2;
    float rij = sqrt((xd * xd + yd * yd) / 10.0);
    return ceil(rij);
}

void constructTSP(std::string graph, cityType *cities, EdgeMatrix *dist)
{
    std::ifstream infile(("instances/" + graph + ".tsp").c_str());
    std::string line;
    bool euclidean = true;

    int city;
    float x, y;
    bool reading_nodes = false;
    while (std::getline(infile, line))
    {
        istringstream iss(line);
        string word;
        if (!reading_nodes)
        {
            iss >> word;
            if (word.compare("EDGE_WEIGHT_TYPE") == 0)
            {
                iss >> word >> word;
                euclidean = !word.compare("EUC_2D");
            }
            else if (word.compare("NODE_COORD_SECTION") == 0)
            {
                reading_nodes = true;
            }
        }
        else if (iss >> city >> x >> y)
        {
            cities[city - 1].x = x;
            cities[city - 1].y = y;
        }
    }
    infile.close();

    for (int from = 0; from < MAX_NODES; from++)
    {
        (*dist)[from][from] = 0.0;

        float edge_dist;
        for (int to = from + 1; to < MAX_NODES; to++)
        {
            if (euclidean)
            {
                edge_dist = euclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
            }
            else
            {
                edge_dist = pseudoEuclideanDistance(cities[from].x, cities[to].x, cities[from].y, cities[to].y);
            }
            if (edge_dist == 0)
            {
                edge_dist = 1.0;
            }
            (*dist)[from][to] = edge_dist;
            (*dist)[to][from] = edge_dist;
        }
    }
}

int bestSolution()
{
    for (int ant = 0; ant < MAX_ANTS; ant++)
    {
        if (antColony[ant].tourLength < best)
        {
            best = antColony[ant].tourLength;
            bestIndex = ant;
        }
    }
    printf("new best: %1.f\n", best);
    return best;
}

int main()
{
    cityType cities[MAX_NODES];
    string graph = "lin105";
    dist = new EdgeMatrix();
    constructTSP(graph, cities, dist);

    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
    QueryPerformanceCounter((LARGE_INTEGER *)&time0);
    init();
    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
#pragma omp parallel for num_threads(NUM_THREAD)
        for (int ant = 0; ant < MAX_ANTS; ant++)
        {
            while (antColony[ant].pathIndex < MAX_NODES)
            {
                antColony[ant].nextNode = selectNextCity(ant);
                antColony[ant].tabu[antColony[ant].nextNode] = 1;
                antColony[ant].path[antColony[ant].pathIndex++] = antColony[ant].nextNode;
                antColony[ant].tourLength += (*dist)[antColony[ant].curNode][antColony[ant].nextNode];
                antColony[ant].curNode = antColony[ant].nextNode;
            }
            antColony[ant].tourLength += (*dist)[antColony[ant].path[MAX_NODES - 1]][antColony[ant].path[0]];
        }
        bestSolution();
        updatePheromone();
        if (i != MAX_ITERATIONS - 1)
        {
            restartAnts();
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&time1);
    cout << "CPU计算时间为" << (time1 - time0) * 1000.0 / freq << "ms" << endl;
}
