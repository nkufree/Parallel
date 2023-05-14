#include <arm_neon.h>
#include <iostream>
#include<cmath>
#include <pthread.h>
#include <sys/time.h>

using namespace std;

int cityx[1024] = {};
int cityy[1024] = {};
int dis[1024][1024] = {};
int ssize = 0;
int num = 0;
int num_threads = 0;

void jisuan(int start, int end)
{

    for (int i = start; i < end; i++)
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

void* worker_thread(void* arg) {
    pthread_t thread_id = pthread_self();
    int* a = (int*)arg;
    int start = ssize * *a;
    int end = ssize * (*a + 1) > num ? num : ssize * (*a + 1);
    delete a;
    for (int k = 0; k < 100; k++)
        jisuan(start,end);
    return NULL;
}

void para()
{
    pthread_t* threads = new pthread_t[num_threads];
    for (int i = 0; i < num_threads; i++) {
        int* a = new int;
        *a = i;
        pthread_create(&threads[i], NULL, worker_thread, a);
    }
    int* m = new int;
    *m = num_threads;
    worker_thread(m);
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    delete[] threads;
}

int main() {
    for (int i = 0; i < 1024; i++)
    {
        cityx[i] = i;
        cityy[i] = 1024 - i;
    }
    struct timeval time0, time1;
    int thread[7] = {1,2,3,4,5,6,7 };
    int a[10] = { 100,200,300,400,500,600,700,800,900,1000 };
    for (int i = 0; i < 10; i++)
    {
        num = a[i];
        cout << num << "\t";
        gettimeofday(&time0, NULL);
        for (int k = 0; k < 100; k++)
            jisuan(0,num);
        gettimeofday(&time1, NULL);
        cout << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << "\t";
        for (int j = 0; j < 7; j++)
        {
            num_threads = thread[j];
            ssize = num / (num_threads+1);
            gettimeofday(&time0, NULL);
            para();
            gettimeofday(&time1, NULL);
            cout << (time1.tv_sec - time0.tv_sec) * 1000000 + time1.tv_usec - time0.tv_usec << "us" << "\t";
        }
        cout << endl;
    }
    return 0;
}

