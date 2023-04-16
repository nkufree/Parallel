#include <emmintrin.h>
#include<Windows.h>
#include<cmath>
#include<iostream>

using namespace std;

alignas(32) int dis[1024][1024] = {};
alignas(32) int cityx[1024] = {};
alignas(32) int cityy[1024] = {};



void jisuan(int num)
{
    for (int i = 0; i < num; i++)
    {
        __m128i x1 = _mm_set1_epi32(cityx[i]);
        __m128i y1 = _mm_set1_epi32(cityy[i]);
        for (int j = 0; j < num; j += 4)
        {
            __m128i x2 = _mm_load_si128((__m128i*)&cityx[j]);
            __m128i y2 = _mm_load_si128((__m128i*)&cityy[j]);

            __m128i xdiff = _mm_sub_epi32(x1, x2);
            __m128i ydiff = _mm_sub_epi32(y1, y2);

            __m128i xdiff_squared = _mm_madd_epi16(xdiff, xdiff);
            __m128i ydiff_squared = _mm_madd_epi16(ydiff, ydiff);

            __m128i sum = _mm_add_epi32(xdiff_squared, ydiff_squared);

            __m128 sqrt_result = _mm_sqrt_ps(_mm_castsi128_ps(sum));

            __m128i int_result = _mm_cvtps_epi32(sqrt_result);
            _mm_store_si128((__m128i*) & dis[i + 1][j + 1], int_result);
        }
    }
}

void chuan(int num)
{
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < num; j++)
        {
            int dx = cityx[i] - cityx[j];
            int dy = cityy[i] - cityy[j];
            dis[i][j] = floor(sqrt(dx * dx + dy * dy));
        }
    }
}

int main()
{
    int num[10] = { 100,200,300,400,500,600,700,800,900,1000 };
    for (int i = 0; i < 1024; i++)
    {
        cityx[i] = i;
        cityy[i] = 1024 - i;
    }
    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    for (int i = 0; i < 10; i++)
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&time0);
        for (int loop = 0; loop < 1000; loop++)
            jisuan(num[i]);
        QueryPerformanceCounter((LARGE_INTEGER*)&time1);
        cout << num[i] << " " << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    }
    
    cout << endl << endl;
    for (int i = 0; i < 10; i++)
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&time0);
        for (int loop = 0; loop < 1000; loop++)
            chuan(num[i]);
        QueryPerformanceCounter((LARGE_INTEGER*)&time1);
        cout << num[i] << " " << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    }
    

    return 0;
}
