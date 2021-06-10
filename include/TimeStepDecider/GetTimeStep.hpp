#pragma once 
#include <cmath>
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/CoreGridData.hpp>
#include <CalculatedResult/FluidResult.hpp>

#if defined(SphMethod)
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
#elif defined(GodunovMethod)
inline real GetNextTimeStep(real dx, const CoreGridData& data)
{
    real max_u_minus_c_and_u_plus_c = 1e-5;
    real dt;
#pragma omp parallel
    {
        real tempThread_max_u_minus_c_and_u_plus_c = 1e-5;
        real ci;
        real ui;
        int i;
#pragma omp for
        for(i=0;i<data.number;++i)
        {
            ci = data.GetSpeedOfSound(i);
            ui = data.GetVelocity(i,1);
            tempThread_max_u_minus_c_and_u_plus_c = std::max(tempThread_max_u_minus_c_and_u_plus_c, std::abs(ui-ci));
            tempThread_max_u_minus_c_and_u_plus_c = std::max(tempThread_max_u_minus_c_and_u_plus_c, std::abs(ui+ci));
        }

#pragma omp critical
        {
            max_u_minus_c_and_u_plus_c = std::max(max_u_minus_c_and_u_plus_c,tempThread_max_u_minus_c_and_u_plus_c);
        }
    }

    return Ccfl*dx/max_u_minus_c_and_u_plus_c;
}
#endif
