#pragma once
#include <Settings.hpp>
#include <cmath>
#include <PhysicsValue/CoreGridData.hpp>

struct RectangleGridData:public CoreGridData
{
    int number[DIM];
    real dr[DIM];

    RectangleGridData(int* num)
    {
        for(int d = 0;d<DIM;++d)
            number[d] = num[d];
    }
};