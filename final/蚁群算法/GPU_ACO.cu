#include <assert.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include<Windows.h>

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

using namespace std;
//内存数据
float* heuristic;
double* phero;
struct ant antColony[ANTS];
float bestSol[ANTS];
float globalBest = TOTAL_DIST;
curandState state[ANTS];
const size_t heuristic_size = sizeof(float) * size_t(NODES * NODES);
const size_t phero_size = sizeof(double) * size_t(NODES * NODES);
//GPU内存数据
float* heuristic_d;
double* phero_d;
struct ant* antColony_d;
float* bestSol_d;
curandState* state_d;
int BLOCKS, THREADS;

__global__ void initializeAnts(struct ant* antColony_d, curandState* state_d,
    float* bestSol_d);
__global__ void setuCurandStates(curandState* stated_d, unsigned long t,
    float* bestSol_d);
__global__ void restartAnts(struct ant* antColony_d, curandState* state_d,
    float* bestSol_d);
__global__ void constructSolution(struct ant* antColony_d, curandState* state_d, float* heuristic_d, double* phero_d);
__global__ void atomicUpdate(struct ant* antColony_d, double* phero_d);
__device__ double probFunctionProduct(int from, int to, double* phero_d, float* heuristic_d);
__device__ int NextNode(struct ant* antColony_d, int pos, float* heuristic_d, double* phero_d, curandState* state_d);

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

__global__ void setupCurandStates(curandState* state_d, unsigned long t) {
    int gid = blockDim.x * blockIdx.x + threadIdx.x;
    curand_init(t, gid, 0, &state_d[gid]);
}

__global__ void initializeAnts(struct ant* antColony_d, curandState* state_d, float* bestSol_d) {

    int ant_id = blockDim.x * blockIdx.x + threadIdx.x;
    for (int node = 0; node < NODES; node++) {

        antColony_d[ant_id].tabu[node] =
            0;
        antColony_d[ant_id].solution[node] =
            -1;
    }
    bestSol_d[ant_id] = (float)TOTAL_DIST;
    antColony_d[ant_id].curNode = curand(&state_d[ant_id]) % NODES;
    antColony_d[ant_id].solution[0] = antColony_d[ant_id].curNode;
    antColony_d[ant_id].tabu[antColony_d[ant_id].curNode] =
        1;
    antColony_d[ant_id].nextNode = -1;
    antColony_d[ant_id].solutionLen = 0;
    antColony_d[ant_id].pathIndex = 1;
}

__global__ void restartAnts(struct ant* antColony_d, curandState* state_d,
    float* bestSol_d) {

    int ant_id = blockDim.x * blockIdx.x + threadIdx.x;

    for (int node = 0; node < NODES; node++) {
        antColony_d[ant_id].tabu[node] =
            0;
        antColony_d[ant_id].solution[node] =
            -1;
    }
    if (antColony_d[ant_id].solutionLen < bestSol_d[ant_id] &&
        antColony_d[ant_id].solutionLen > 0) {
        bestSol_d[ant_id] = antColony_d[ant_id].solutionLen;

    }
    antColony_d[ant_id].curNode = curand(&state_d[ant_id]) % NODES;
    antColony_d[ant_id].solution[0] = antColony_d[ant_id].curNode;
    antColony_d[ant_id].tabu[antColony_d[ant_id].curNode] =
        1; 
    antColony_d[ant_id].nextNode = -1;
    antColony_d[ant_id].solutionLen = 0;
    antColony_d[ant_id].pathIndex = 1;
}

void acoSolve() {
    int iteration = 0;
    while (iteration++ < MAX_ITERATIONS) {
        constructSolution << <BLOCKS, THREADS >> > (antColony_d, state_d, heuristic_d,
            phero_d);

        cudaDeviceSynchronize();
        // Move solution back to Host
        cudaMemcpy(antColony, antColony_d, sizeof(antColony),
            cudaMemcpyDeviceToHost);
        for (int from = 0; from < NODES; from++) {
            for (int to = 0; to < NODES; to++) {
                if (from != to) {
                    phero[from + to * NODES] *= (1.0 - RHO);
                    if (phero[from + to * NODES] < 0.0) {
                        phero[from + to * NODES] = PHERO_INITIAL;
                    }
                }
            }
        } 

        cudaMemcpy(phero_d, phero, phero_size, cudaMemcpyHostToDevice);
        cudaMemcpy(bestSol, bestSol_d, sizeof(bestSol), cudaMemcpyDeviceToHost);
        atomicUpdate << <BLOCKS, THREADS >> > (antColony_d, phero_d);

        for (int i = 0; i < ANTS; i++) {
            if (bestSol[i] < globalBest) {
                globalBest = bestSol[i];
            }
        }

        restartAnts << <BLOCKS, THREADS >> > (antColony_d, state_d, bestSol_d);
        cudaDeviceSynchronize();

    }

    printf("Best Solution %f ", globalBest);
}

__global__ void atomicUpdate(struct ant* antColony_d, double* phero_d) {

    int ant_id = blockDim.x * blockIdx.x + threadIdx.x;
    int from, to;
    for (int i = 0; i < NODES; i++) {
        from = antColony_d[ant_id].solution[i];
        if (i > NODES - 1) {
            to = antColony_d[ant_id].solution[i + 1];
        }
        else {
            to = antColony_d[ant_id].solution[0];
        }
        atomicAdd(&phero_d[from + to * NODES], Q / antColony_d[ant_id].solutionLen * RHO);
        atomicAdd(&phero_d[from + to * NODES], Q / antColony_d[ant_id].solutionLen * RHO);
    }
}
__global__ void constructSolution(struct ant* antColony_d, curandState* state_d,
    float* heuristic_d, double* phero_d) {

    int ant_id = blockDim.x * blockIdx.x + threadIdx.x;
    int node = 0;

    while (node++ < NODES) {
        if (antColony_d[ant_id].pathIndex < NODES) {
            antColony_d[ant_id].nextNode =
                NextNode(antColony_d, ant_id, heuristic_d, phero_d, state_d);
            antColony_d[ant_id].tabu[antColony_d[ant_id].nextNode] = 1;
            antColony_d[ant_id].solution[antColony_d[ant_id].pathIndex++] =
                antColony_d[ant_id].nextNode;
            antColony_d[ant_id].solutionLen +=
                heuristic_d[antColony_d[ant_id].curNode +
                (antColony_d[ant_id].nextNode * NODES)];

            if (antColony_d[ant_id].pathIndex == NODES) {
                antColony_d[ant_id].solutionLen +=
                    heuristic_d[antColony_d[ant_id].solution[NODES - 1] +
                    (antColony_d[ant_id].solution[0] * NODES)];
            }
            antColony_d[ant_id].curNode = antColony_d[ant_id].nextNode;
        }
    }
}

__device__ double probFunctionProduct(int from, int to, double* phero_d,
    float* heuristic_d) {
    double result;
    result = pow(phero_d[from + to * NODES], ALPHA) *
        pow(1 / (heuristic_d[from + to * NODES]), BETA);
    if (!isnan(result)) {
        return (double)((result));
    }
    else {
        return 0;
    }
}

__device__ int NextNode(struct ant* antColony_d, int pos, float* heuristic_d,
    double* phero_d, curandState* state_d) {
    int to, from;
    double denom = 0.00000001;
    from = antColony_d[pos].curNode;
    for (to = 0; to < NODES; to++) {
        if (antColony_d[pos].tabu[to] == 0) {
            denom += probFunctionProduct(from, to, phero_d, heuristic_d);
        }
    }
    assert(denom != 0.0);
    to++;
    int count = NODES - antColony_d[pos].pathIndex;
    do {
        double p;
        to++;
        if (to >= NODES)
            to = 0;
        if (antColony_d[pos].tabu[to] ==
            0) {
            p = probFunctionProduct(from, to, phero_d, heuristic_d) / denom;
            double x = (double)(curand(&state_d[pos]) % 1000000000) / 1000000000.0;
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

int main() {
    float exec_time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    if (ANTS <= 1024) {
        BLOCKS = 1;
        THREADS = ANTS;
    }
    else {
        THREADS = 1024;
        BLOCKS = ceil(ANTS / (float)THREADS);
    }
    heuristic = (float*)malloc(NODES * NODES * sizeof(float));
    phero = (double*)malloc(NODES * NODES * sizeof(double));

    nodeTSP nodes[NODES];
    constructTSP("lin105", nodes);
    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&time0);
    cudaMalloc((void**)&antColony_d, sizeof(antColony));
    cudaMalloc((void**)&state_d, sizeof(state));
    cudaMalloc((void**)&bestSol_d, sizeof(bestSol));
    cudaMalloc((void**)&heuristic_d, heuristic_size);
    cudaMalloc((void**)&phero_d, phero_size);

    cudaMemcpy(heuristic_d, heuristic, heuristic_size, cudaMemcpyHostToDevice);
    cudaMemcpy(phero_d, phero, phero_size, cudaMemcpyHostToDevice);

    time_t t;
    time(&t);
    setupCurandStates << <BLOCKS, THREADS >> > (state_d, (unsigned long)t);
    cudaDeviceSynchronize();
    initializeAnts << <BLOCKS, THREADS >> > (antColony_d, state_d, bestSol_d);
    cudaDeviceSynchronize();

    cudaEventRecord(start, 0);
    acoSolve();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&exec_time, start, stop);

    printf("%5.5f \n", exec_time / 1000);
    free(phero);
    free(heuristic);

    cudaFree(antColony_d);
    cudaFree(heuristic_d);
    cudaFree(phero_d);
    cudaFree(state_d);
    cudaFree(bestSol_d);
    QueryPerformanceCounter((LARGE_INTEGER*)&time1);
    cout << "GPU加速算法时间为" << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    return 0;
}
