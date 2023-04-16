#include <sys/time.h>
#include <arm_neon.h>
#include <iostream>
#include <cmath>

using namespace std;

alignas(32) int cityx[1024] = {};
alignas(32) int cityy[1024] = {};
alignas(32) int dis[1024][1024] = {};

void jisuan(int num)
{
    for (int i = 0; i < num; i++)
    {
        int32x4_t x1 = vld1q_s32(cityx + i);
        int32x4_t y1 = vld1q_s32(cityy + i);
        for (int j = 0; j < num; j += 4)
        {
            int32x4_t x2 = vld1q_s32(cityx + j);
            int32x4_t y2 = vld1q_s32(cityy + j);
            int32x4_t dx = vsubq_s32(x1, x2);
            int32x4_t dy = vsubq_s32(y1, y2);
            int32x4_t xx = vmulq_s32(dx, dx);
            int32x4_t yy = vmulq_s32(dy, dy);
            int32x4_t sum = vaddq_s32(xx, yy);
            float32x4_t distance = vrsqrteq_f32(vreinterpretq_f32_s32(sum));
            int32x4_t uint_distance = vcvtq_s32_f32(distance);
            int32_t* dis_ptr = &dis[i + 1][j + 1];
            vst1q_s32(dis_ptr, uint_distance);
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
            dis[i][j] = ceil(sqrt(dx * dx + dy * dy));
        }
    }
}

int main()
{
    int num[10] = {100,200,300,400,500,600,700,800,900,1000};
    for (int i = 0; i < 1024; i++)
    {
        cityx[i] = i;
        cityy[i] = 1024 - i;
    }
    struct timeval time0, time1;

    for (int i = 0; i < 10; i++)
    {
        gettimeofday(&time0, NULL);
        for (int loop = 0; loop < 1000; loop++)
            jisuan(num[i]);
        gettimeofday(&time1, NULL);
        cout << num[i] << '\t' << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << endl;
    }
    cout << endl << endl;
    

    for (int i = 0; i < 10; i++)
    {
        gettimeofday(&time0, NULL);
        for (int loop = 0; loop < 1000; loop++)
            chuan(num[i]);
        gettimeofday(&time1, NULL);
        cout << num[i] << '\t' << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << endl;
    }
    return 0;
}