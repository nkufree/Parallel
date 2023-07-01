#define MAX_NODES 1889 
#define MAX_DIST 10000
#define MAX_TOUR (MAX_NODES * MAX_DIST)
#define MAX_ANTS 300
#define NUM_EDGES ((MAX_NODES * MAX_NODES - MAX_NODES) / 2)

struct cityType {
  float x,y;
};

class EdgeMatrix {
  float *dist;
public:
  EdgeMatrix() {
    dist = new float[MAX_NODES * MAX_NODES];
  }
  ~EdgeMatrix() {
    delete dist;
  }
  float* operator[](unsigned int i) {
   return &dist[MAX_NODES * i];
 }

  float *get_array(){
    return dist;
  }
};

#define ALPHA 2 
#define BETA 10 
#define RHO 0.5 
#define QVAL 10 
#define MAX_ITERATIONS 10
#define MAX_TIME (MAX_ITERATIONS * MAX_NODES)
#define INIT_PHER (1.0 / MAX_NODES)
