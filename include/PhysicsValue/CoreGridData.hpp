#pragma once
#include <Settings.hpp>

struct CoreGridData
{
    int number;
    real position[N][DIM];
    real velocity[N][DIM];
    real density[N];
    real pressure[N];

    CoreGridData(int num):number(num) {}
};