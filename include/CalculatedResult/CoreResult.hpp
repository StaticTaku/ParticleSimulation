#pragma once
#include <Settings.hpp>

struct CoreResult
{
    real _accel[N][DIM];
    real (*accel)[DIM];

    CoreResult():accel(_accel) {}

    void ClearAccel(int number)
    {
        #pragma omp parallel
        {
        #pragma omp for
            for(int i = 0;i<number;++i)
                for(int d = 0;d<DIM;++d)
                    accel[i][d] = 0;
        }
    }
};