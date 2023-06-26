#include <assert.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <random>
#include <ctime>
#include<Windows.h>

using namespace std;

#define ANTS 1024
#define ALPHA 2
#define BETA 10
#define RHO 0.5
#define Q 50
#define MAX_ITERATIONS 10

#define NODES 105
#define DIST 10000
#define PHERO_INITIAL (1.0 / NODES)
#define TOTAL_DIST (DIST * NODES)

struct ant {
    int curNode, nextNode, pathIndex;
    int tabu[NODES];
    int solution[NODES];
    float solutionLen;
};

struct nodeTSP {
    float x, y;
};

float heuristic[NODES * NODES];
double phero[NODES * NODES];
ant antColony[ANTS];
float bestSol[ANTS];
float globalBest = TOTAL_DIST;

float euclideanDistance(float x1, float x2, float y1, float y2) {
    float xd = x1 - x2;
    float yd = y1 - y2;
    return (float)(sqrt(xd * xd + yd * yd));
}

void constructTSP(string graph, nodeTSP* nodes) {
    ifstream infile(("instances/" + graph + ".tsp").c_str());
    string line;
    bool euclidean = true;
    int node;
    float x, y;
    bool reading_nodes = false;

    while (getline(infile, line)) {
        istringstream iss(line);
        string word;
        if (!reading_nodes) {
            iss >> word;
            if (word.compare("EDGE_WEIGHT_TYPE") == 0) {
                iss >> word >> word;
                cout << "edge type: " << word << endl;
                euclidean = !word.compare("EUC_2D");
            }
            else if (word.compare("NODE_COORD_SECTION") == 0) {
                reading_nodes = true;
            }
        }
        else if (iss >> node >> x >> y) {
            nodes[node - 1].x = x;
            nodes[node - 1].y = y;
        }
    }
    infile.close();
    for (int from = 0; from < NODES; from++) {
        for (int to = from + 1; to < NODES; to++) {
            float edge_weight;
            if (euclidean) {
                edge_weight = euclideanDistance(nodes[from].x, nodes[to].x,
                    nodes[from].y, nodes[to].y);
            }

            if (edge_weight == 0) {
                edge_weight = 1.0;
            }
            heuristic[from + to * NODES] = edge_weight;
            heuristic[to + from * NODES] = edge_weight;
            phero[from + to * NODES] = PHERO_INITIAL;
            phero[to + from * NODES] = PHERO_INITIAL;
        }
    }
}

double probFunctionProduct(int from, int to) {
    double result;
    result = pow(phero[from + to * NODES], ALPHA) *
        pow(1 / (heuristic[from + to * NODES]), BETA);
    if (!isnan(result)) {
        return (double)((result));
    }
    else {
        return 0;
    }
}

int NextNode(int pos) {
    int to, from;
    double denom = 0.00000001;
    from = antColony[pos].curNode;
    for (to = 0; to < NODES; to++) {
        if (antColony[pos].tabu[to] == 0) {
            denom += probFunctionProduct(from, to);
        }
    }
    assert(denom != 0.0);
    to++;
    int count = NODES - antColony[pos].pathIndex;
    do {
        double p;
        to++;
        if (to >= NODES)
            to = 0;
        if (antColony[pos].tabu[to] == 0) { 
            p = probFunctionProduct(from, to) / denom;
            srand((unsigned)time(NULL));
            double x = (double)(rand() % 1000000000) / 1000000000.0;
            if (x < p) {
                break;
            }
            count--;
            if (count == 0) {
                break;
            }
        }
    } while (1);
    return to;
}

void initializeAnts() {
    for (int i = 0; i < ANTS; i++) {
        for (int node = 0; node < NODES; node++) {
            antColony[i].tabu[node] = 0;
            antColony[i].solution[node] = -1;
        }
        bestSol[i] = (float)TOTAL_DIST;
        srand((unsigned)time(NULL));
        antColony[i].curNode = rand() % NODES;
        antColony[i].nextNode = -1;
        antColony[i].pathIndex = 1;
        antColony[i].solutionLen = 0.0;

        antColony[i].tabu[antColony[i].curNode] = 1;
        antColony[i].solution[0] = antColony[i].curNode;
    }
}

void updateAnts() {
    for (int i = 0; i < ANTS; i++) {
        int node = 0;

        while (node++ < NODES) {
            if (antColony[i].pathIndex < NODES) {
                antColony[i].nextNode =
                    NextNode(i);
                antColony[i].tabu[antColony[i].nextNode] = 1;
                antColony[i].solution[antColony[i].pathIndex++] =
                    antColony[i].nextNode;
                antColony[i].solutionLen +=
                    heuristic[antColony[i].curNode +
                    (antColony[i].nextNode * NODES)];
                if (antColony[i].pathIndex == NODES) {
                    antColony[i].solutionLen +=
                        heuristic[antColony[i].solution[NODES - 1] +
                        (antColony[i].solution[0] * NODES)];
                }
                antColony[i].curNode = antColony[i].nextNode;
            }
        }
    }
}

void updatePheromone() {
    for (int from = 0; from < NODES; from++) {
        for (int to = 0; to < NODES; to++) {
            if (from != to) {
                phero[from + to * NODES] *= (1.0 - RHO);
                if (phero[from + to * NODES] < 0.0)
                    phero[from + to * NODES] = PHERO_INITIAL;
            }
        }
    }

    for (int i = 0; i < ANTS; i++) {
        for (int j = 0; j < NODES; j++) {
            if (j < NODES - 1) {
                int from = antColony[i].solution[j];
                int to = antColony[i].solution[j + 1];
                phero[from + to * NODES] += (Q / antColony[i].solutionLen);
                phero[to + from * NODES] += (Q / antColony[i].solutionLen);
            }
            else {
                int from = antColony[i].solution[j];
                int to = antColony[i].solution[0];
                phero[from + to * NODES] += (Q / antColony[i].solutionLen);
                phero[to + from * NODES] += (Q / antColony[i].solutionLen);
            }
        }
    }

    for (int i = 0; i < ANTS; i++) {
            if (bestSol[i] < globalBest) {
                globalBest = bestSol[i];
            }
        }
}

void restartAnts()
{
    for (int i = 0; i < ANTS; i++)
    {
        for (int node = 0; node < NODES; node++) {
            antColony[i].tabu[node] = 0;
            antColony[i].solution[node] = 1;
        }
        if (antColony[i].solutionLen < bestSol[i] &&
            antColony[i].solutionLen > 0) {
            bestSol[i] = antColony[i].solutionLen;
        }
        srand((unsigned)time(NULL));
        antColony[i].curNode = rand() % NODES;
        antColony[i].solution[0] = antColony[i].curNode;
        antColony[i].tabu[antColony[i].curNode] = 1;
        antColony[i].nextNode = -1;
        antColony[i].solutionLen = 0;
        antColony[i].pathIndex = 1;
    }
}

void search() {
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        updateAnts();
        updatePheromone();
        restartAnts();
    }
}

int main() {
    srand(time(NULL));

    nodeTSP nodes[NODES];
    ant* ants = new ant[ANTS];
    float* pheromone = new float[NODES * NODES];
    constructTSP("lin105", nodes);
    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&time0);
    initializeAnts();
    search();
    QueryPerformanceCounter((LARGE_INTEGER*)&time1);
    cout << "串行算法时间为" << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    float bestLen = FLT_MAX;
    int bestIndex = -1;

    for (int i = 0; i < ANTS; i++) {
        if (ants[i].solutionLen < bestLen) {
            bestLen = ants[i].solutionLen;
            bestIndex = i;
        }
    }

    cout << "Best length: " << bestLen << endl;
    cout << "Best tour: ";

    for (int i = 0; i < NODES; i++) {
        cout << ants[bestIndex].solution[i] << " ";
    }

    cout << endl;

    delete[] ants;
    delete[] pheromone;

    return 0;
}
