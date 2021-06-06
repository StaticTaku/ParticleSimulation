#pragma once 
#include <cmath>
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <CalculatedResult/FluidResult.hpp>

inline real GetNextTimeStep(FluidResult& result,const SphDataWithGamma& data)
{
    real delta_tf = 1e5;
    real delta_tcv = 1e5;
#pragma omp parallel
    {
        real tempThread_delta_tf = 1e5;
        real tempThread_delta_tcv = 1e5;
        real ci;
        real fi = 0;
        int i;
#pragma omp for
        for(i=0;i<data.number;++i)
        {
            ci = sqrt(heatCapRatio*data.pressure[i]/data.density[i]);
            for(int d = 0;d<DIM;++d)
                fi = data.accel[i][d]*data.accel[i][d];
            tempThread_delta_tf = std::min(tempThread_delta_tf,data.h/std::abs(std::sqrt(fi)));
            tempThread_delta_tcv = std::min(tempThread_delta_tcv,data.h/(ci + 0.6*(alpha*ci + cbeta*result.max_ui[i])));
            result.max_ui[i] = 0;
        }

#pragma omp critical
        {
            delta_tf = std::min(delta_tf,tempThread_delta_tf);
            delta_tcv = std::min(delta_tcv,tempThread_delta_tcv);
        }
    }

    return Ccfl*std::min(delta_tf,delta_tcv);
}