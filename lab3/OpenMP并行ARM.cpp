#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include <queue>
#include <omp.h>
#include <sys/time.h>

using namespace std;


struct node {
    int num, x, y;
};


struct solution {
    vector<int> path;
    double sum, F, P;
    int gen;
};

int taskNum = 1;
int num_threads = 2;
int curGen = 1; 
bool isFound = false;  
const int maxGen = 10000; 
const int start = 1;    
int packNum = 20;    
double mutationP = 0.04;  
double crossP = 0.7;   
double curCost = 0; 
int Hayflick = 6;   
vector<node> city;  
vector<int> num;    
vector<solution> pack;  
double dis[1000][1000] = {};   
queue<double> record;  

regex pattern("^(\\d+) (\\d+) (\\d+)$");


double distanceBetween(node &a, node &b) {
    int x = a.x - b.x;
    int y = a.y - b.y;
    return ceil(sqrt((x * x + y * y) / 10.0));
}


double sum(vector<int> &v) {
    double sum = 0;
    for (auto it = v.begin() + 1; it != v.end(); it++) {
        sum += dis[*(it - 1)][*it];
    }
    return sum;
}


bool cmp(const solution &a, const solution &b) {
    return a.sum < b.sum;
}


void initData() {
    /*
    ifstream f("att48.tsp", ios::in);


    if (!f) {
        cout << "File Not Found" << endl;
    }
    string s;
    while (getline(f, s)) {

        if (regex_match(s, pattern)) {
            sregex_iterator it(s.begin(), s.end(), pattern);
            city.push_back(node{stoi(it->str(1)), stoi(it->str(2)), stoi(it->str(3))});
        }
    }
    f.close();*/
    for(int i = 1; i < 49; i++)
    {
        city.push_back(node{i,i*100,i*200});
    }

    for (int i = 0; i < city.size(); i++)
    {
        for (int j = 0; j < city.size(); j++)
        {
            int x = city[i].x - city[j].x;
            int y = city[i].y - city[j].y;
            dis[i + 1][j + 1] = ceil(sqrt((x * x + y * y) / 10.0));

        }
    }
    for (auto it = city.begin(); it != city.end(); it++) {
        if (it->num != start) {
            num.push_back(it->num);
        }
    }
}


void initPack(int gen) {
    for (int i = 0; i < packNum; i++) {

        solution temp;
        temp.path.push_back(start);
        random_shuffle(num.begin(), num.end());
        for (auto it = num.begin(); it != num.end(); it++) {
            temp.path.push_back(*it);
        }
        temp.path.push_back(start);


        temp.gen = gen;
        temp.sum = sum(temp.path);
        pack.push_back(temp);
    }
    if (gen == 1) {
        solution greed;
        unordered_map<int, bool> isVisited;
        greed.path.push_back(start);
        isVisited[start] = true;
        for (int i = 0; i < num.size(); i++) {
            int lastCity = greed.path.back();
            int min = -1;
            double minDistance = 999999;
            for (int i = 1; i <= city.size(); i++) {
                if (isVisited.count(i) == 0 && (min == -1 || dis[lastCity][i] < minDistance)) {
                    min = i;
                    minDistance = dis[lastCity][i];
                }
            }
            greed.path.push_back(min);
            isVisited[min] = true;
        }
        greed.path.push_back(start);
        greed.gen = gen;
        greed.sum = sum(greed.path);
        pack[0] = greed;
    }
}


void passOn() {
    record.push(pack[0].sum);

    if (record.size() > 500) {
        record.pop();

        if (curCost != pack[0].sum && record.front() == record.back()) {

            vector<int> Prometheus = pack[0].path;
            int gen = pack[0].gen;
            double sum = pack[0].sum;
            pack.clear();

            if (Hayflick > 0) {
                Hayflick--;
                packNum *= log(packNum) / log(10);
            }

            initPack(gen);
            pack[0].path = Prometheus;
            pack[0].sum = sum;
            sort(pack.begin(), pack.end(), cmp);
            curCost = sum;
            mutationP += 0.1;

            while (!record.empty()) {
                record.pop();
            }
        }
    }
}


solution cross(solution &firstParent, solution &secondParent) {

    int length = int((rand() % 1000 / 1000.0) * city.size());
    int off = rand() % city.size() + 1;
    vector<int> nextGen(firstParent.path.size());
    unordered_map<int, bool> selected;
    nextGen[0] = start;
    for (int i = off; i < nextGen.size() - 1 && i < off + length; i++) {
        nextGen[i] = firstParent.path[i];
        selected[nextGen[i]] = true;
    }
    for (int i = 1, j = 1; i < nextGen.size(); i++) {
        if (nextGen[i] == 0) {
            for (; j < secondParent.path.size(); j++) {
                if (selected.count(secondParent.path[j]) == 0) {
                    nextGen[i] = secondParent.path[j];
                    selected[secondParent.path[j]] = true;
                    break;
                }
            }
        }
    }
    return solution{nextGen, sum(nextGen), 0, 0, firstParent.gen + 1};
}


void mutation(solution &cur) {

    int length = int((rand() % 1000 / 1000.0) * city.size());
    int off = rand() % city.size() + 1;
    vector<int> m;
    unordered_map<int, bool> selected;
    m.push_back(start);
    for (int i = off; i < cur.path.size() - 1 && i < off + length; i++) {
        m.push_back(cur.path[i]);
        selected[cur.path[i]] = true;
    }
    for (int i = 1; i < cur.path.size(); i++) {
        if (selected.count(cur.path[i]) == 0) {
            m.push_back(cur.path[i]);
        }
    }
    for (int i = 0; i < m.size(); i++) {
        cur.path[i] = m[i];
    }
    cur.sum = sum(cur.path);
}


void gmutation(solution &cur) {

    int selected = rand() % (city.size() - 4) + 2, min = 1;
    int selectedCity = cur.path[selected];
    int begin = 0, end = 0;
    for (int i = 1; i <= city.size(); i++)
    {
        if (i != selectedCity && dis[selectedCity][i] < dis[selectedCity][min])
        {
            min = i;
        }
    }
    for (int i = 1; i < cur.path.size() - 1; i++) {
        if (cur.path[i] == min) {
            if (i > selected + 1) {
                begin = selected + 1;
                end = i;
            } else if (i < selected - 1) {
                begin = i;
                end = selected - 1;
            }
            break;
        }
    }
    vector<int> stack;
    for (int i = begin; i <= end; i++) {
        stack.push_back(cur.path[i]);
    }
    for (int i = begin; i <= end; i++) {
        cur.path[i] = stack.back();
        stack.pop_back();
    }
    cur.sum = sum(cur.path);
}


vector<solution> process() {
    vector<solution> nextGenPack;
    sort(pack.begin(), pack.end() - 1, cmp);

    if (pack[0].sum == 10628) {
        isFound = true;
    }
    passOn();   

    nextGenPack.push_back(pack[0]);
    nextGenPack[0].gen++;

    int cycleNum = pack.size() - 1;
    #pragma omp parallel for num_threads(num_threads) schedule(static,taskNum)
    for (int i = 0; i < cycleNum; i++) {
        int tid = omp_get_thread_num();
        int start = tid * (cycleNum) / (num_threads * taskNum);
        int end = min((tid + 1) * cycleNum / (num_threads * taskNum),cycleNum);

        double total = 0;

        for (int i = start; i < end; i++) {
            pack[i].F = 10000 / pack[i].sum;
            pack[i].P = (i == start ? 0 : pack[i - 1].P) + pack[i].F;
            total += pack[i].F;
        }
        if (rand() % 10000 / 10000.0 < crossP) {
            double selected = (rand() % 10000 / 10000.0) * total;
            
            for (int j = start; j < end - 1; j++) {
                if (selected < pack[j].P) {
                    solution temp = cross(pack[i], pack[j]);
                    #pragma omp critical
                    {
                        nextGenPack.push_back(temp);
                    }
                    break;
                }
            }
        }
        else
        {
            pack[i].gen++;
            #pragma omp critical
            {
                nextGenPack.push_back(pack[i]);
            }
            
        }
        if (rand() % 10000 / 10000.0 < mutationP) {
            #pragma omp critical
            {
                Hayflick < 6 ? gmutation(nextGenPack.back()) : mutation(nextGenPack.back());
            }
        }
    }

    return nextGenPack;
}

void ceshi()
{
    curGen = 1; 
    isFound = false;  
    curCost = 0; 
    Hayflick = 6;  
    city.clear(); 
    num.clear();    
    pack.clear();  
    srand(unsigned(time(NULL)));    
    initData();     
    initPack(1);    
    struct timeval time0, time1;
    gettimeofday(&time0, NULL);
    pack = process();
    curGen++;
    gettimeofday(&time1, NULL);
    cout << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << "\t";
}

int main() {
    int thread[7] = { 2,3,4,5,6,7,8 };
    int a[10] = { 100,400,700,1000,1300,1600,1900,2100,2400,2700 };
    for (int i = 0; i < 10; i++)
    {
        packNum = a[i];
        cout << packNum << "\t";
        for (int j = 0; j < 7; j++)
        {
            num_threads = thread[j];
            ceshi();
        }
        cout << endl;
    }
    return 0;
}
