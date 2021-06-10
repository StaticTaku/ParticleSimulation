#pragma once
#include <Settings.hpp>
#include <cmath>

struct CoreGridData
{
    int number;
    real position[N][DIM];
    real density[N];
    real momentum[N][DIM];
    real energy[N];
    CoreGridData(int num):number(num) {}

    real GetVelocity(int i, int dim) const
    {   
        return momentum[i][dim-1]/density[i];
    }

    real GetPressure(int i, int dim) const 
    {
        return (energy[i] - momentum[i][dim-1]*momentum[i][dim-1]/(2*density[i]))*(heatCapRatio-1);
    }

    real GetSpeedOfSound(int i) const
    {
        return std::sqrt(heatCapRatio*GetPressure(i,1)/density[i]);
    }
};