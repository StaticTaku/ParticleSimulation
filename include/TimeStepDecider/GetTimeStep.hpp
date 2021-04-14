#pragma once 
#include <cmath>
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>

constexpr real MIN = 1e-4;
inline real GetNextTimeStep(const SphDataWithGamma& data)
{
    real MinDt = 1e5;
#pragma omp parallel
    {
        real tempThreadMinDt = 1e5;
        int i;
#pragma omp for
        for(i=0;i<data.number;++i)
            tempThreadMinDt = std::min(tempThreadMinDt,sqrt(data.heatCapRatio*data.pressure[i]/data.density[i]));

#pragma omp critical
        {
            MinDt = std::min(MinDt,tempThreadMinDt);
        }
    }

    return std::min(Ccfl*data.h/MinDt,MIN);
}