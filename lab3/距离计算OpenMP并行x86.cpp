#include <immintrin.h>
#include <iostream>
#include<Windows.h>
#include<cmath>
#include <omp.h>

using namespace std;

int cityx[1024] = {};
int cityy[1024] = {};
int dis[1024][1024] = {};
int ssize = 0;
int num = 0;
int num_threads = 0;

void jisuan(int start, int end)
{
    for (int j = start; j < end; j++)
    {
        __m256i x1 = _mm256_set1_epi32(cityx[j]);
        __m256i y1 = _mm256_set1_epi32(cityy[j]);
        for (int i = 0; i < num; i += 8)
        {
            __m256i x2 = _mm256_loadu_si256((__m256i*) & cityx[i]);
            __m256i y2 = _mm256_loadu_si256((__m256i*) & cityy[i]);

            __m256i xdiff = _mm256_sub_epi32(x2, x1);
            __m256i ydiff = _mm256_sub_epi32(y2, y1);
            __m256i ds_vec = _mm256_mullo_epi32(xdiff, xdiff);
            ds_vec = _mm256_add_epi32(ds_vec, _mm256_mullo_epi32(ydiff, ydiff));
            __m256i ds = _mm256_cvtps_epi32(_mm256_sqrt_ps(_mm256_cvtepi32_ps(ds_vec)));
            _mm256_storeu_si256((__m256i*) & (dis[j + 1][i + 1]), ds);
        }
    }
}

int main() {
    for (int i = 0; i < 1024; i++)
    {
        cityx[i] = i;
        cityy[i] = 1024 - i;
    }
    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    int thread[7] = { 2,3,4,5,6,7,8 };
    int a[10] = { 100,200,300,400,500,600,700,800,900,1000 };
    for (int i = 0; i < 10; i++)
    {
        num = a[i];
        cout << num << "\t";
        QueryPerformanceCounter((LARGE_INTEGER*)&time0);
        for (int k = 0; k < 100; k++)
            jisuan(0, num);
        QueryPerformanceCounter((LARGE_INTEGER*)&time1);
        cout << (time1 - time0) * 1000.0 / freq << "ms" << "\t";
        for (int j = 0; j < 7; j++)
        {
            num_threads = thread[j];
            ssize = num / num_threads;
            QueryPerformanceCounter((LARGE_INTEGER*)&time0);
            #pragma omp parallel num_threads(num_threads)
            int tid = omp_get_thread_num();
            int start = tid * num / num_threads;
            int end = (tid + 1) * num / num_threads  > num ? num : (tid + 1) * num / num_threads ;
            for (int k = 0; k < 100; k++)
                jisuan(start, end);
            QueryPerformanceCounter((LARGE_INTEGER*)&time1);
            cout << (time1 - time0) * 1000.0 / freq << "ms" << "\t";
        }
        cout << endl;
    }
    return 0;
}

