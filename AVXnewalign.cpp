
#include <immintrin.h>
#include <iostream>
#include<Windows.h>
#include<cmath>

using namespace std;

alignas(32) float cityx[1024] = {};
alignas(32) float cityy[1024] = {};
alignas(32) float dis[1024][1024] = {};

void jisuan(int num)
{
    __m256 lat1, lat2, lon1, lon2, dlat, dlon, a, b, c, d;

    for (int i = 0; i < num; i++) {
        lat1 = _mm256_load_ps(cityy + i);
        lon1 = _mm256_load_ps(cityx + i);

        for (int j = 0; j < num; j += 4) {
            lat2 = _mm256_load_ps(cityy + j);
            lon2 = _mm256_load_ps(cityx + j);

            dlat = _mm256_sub_ps(lat1, lat2);
            dlon = _mm256_sub_ps(lon1, lon2);

            a = _mm256_sin_ps(_mm256_div_ps(dlat, _mm256_set1_ps(2)));
            a = _mm256_mul_ps(a, a);
            b = _mm256_mul_ps(_mm256_cos_ps(lat1), _mm256_cos_ps(lat2));
            c = _mm256_sin_ps(_mm256_div_ps(dlon, _mm256_set1_ps(2)));
            c = _mm256_mul_ps(c, c);
            c = _mm256_fmadd_ps(b, c, a);
            c = _mm256_sqrt_ps(c);
            d = _mm256_asin_ps(c);
            d = _mm256_mul_ps(_mm256_set1_ps(12742), d); // 地球半径为6371km,两倍为12742km

            _mm256_store_ps(dis[i] + j, d);
        }
    }
}

void chuan(int num)
{
    for (int i = 0; i < num; i++)
    {
        float lat1 = cityx[i] / (2 * 3.1415926);
        float lon1 = cityy[i] / (2 * 3.1415926);
        for (int j = 0; j < num; j++)
        {
            float lat2 = cityx[j] / (2 * 3.1415926);
            float lon2 = cityy[j] / (2 * 3.1415926);
            float dlat = lat1 - lat2;
            float dlon = lon1 - lon2;
            float a = sin(dlat / 2);
            a = a * a;
            float b = cos(cityy[i]) * cos(cityy[j]);
            float c = sin(dlon / 2);
            c = sqrt(c * c);
            b = b * c;
            dis[i][j] = 12742 * asin(a+b);
        }
    }
}

int main()
{
    int num[10] = { 100,200,300,400,500,600,700,800,900,1000 };
    for (int i = 0; i < 1024; i++)
    {
        cityx[i] = i / 10.0;
        cityy[i] = (1024 - i) / 10.0;
    }
    long long time0, time1, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    for (int i = 0; i < 10; i++)
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&time0);
        for (int loop = 0; loop < 100; loop++)
            jisuan(num[i]);
        QueryPerformanceCounter((LARGE_INTEGER*)&time1);
        cout << num[i] << " " << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    }

    cout << endl << endl;
    for (int i = 0; i < 10; i++)
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&time0);
        for (int loop = 0; loop < 100; loop++)
            chuan(num[i]);
        QueryPerformanceCounter((LARGE_INTEGER*)&time1);
        cout << num[i] << " " << (time1 - time0) * 1000.0 / freq << "ms" << endl;
    }

}
