#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/SphCoreData.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <Advancer/VelocityVerlet.hpp>

class FluidValueAdvancer:public VelocityVerlet
{
public:
    void UpdatePosition(SphCoreDataWithFixedH& data, const CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                    data.position[i][j] += (data.velocity[i][j] + 0.5 * data.accel[i][j] * dt )* dt;
            }
        }
    #pragma omp barrier
    }

    void UpdateVelocity(SphCoreDataWithFixedH& data, CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                {
                    data.velocity[i][j] += 0.5 * (data.accel[i][j] + result.accel[i][j]) * dt;
                }
            }
        }
    #pragma omp barrier
    }
    
    void UpdateDensity(SphCoreDataWithFixedH& data, const FluidResult& result)
    {
    #pragma omp parallel
        {
            real h = data.h;
    #pragma omp for
            for(int i = 0;i<data.number;++i)
            {
                data.density[i] = 0;
                for(int j = 0;j<data.number;++j)
                    data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],h);
            }
        }
    }

    void PredictInternalEnergyAndPressure(SphDataWithGamma& data, const FluidResult& result, real dt)
    {
    #pragma omp parallel
        {
            real heatCapRatio_1 = data.heatCapRatio-1;
            real predictedInternalEnergy;
    #pragma omp for
            for(int i = 0;i<data.number;++i)
            {
                predictedInternalEnergy = data.internalEnergy[i] + result.internalEnergyDif[i] * dt;
                data.pressure[i] = predictedInternalEnergy*(heatCapRatio_1*data.density[i]);
            }
        }
    }

    void UpdateInternalEnergyAndPressure(SphDataWithGamma& data, const FluidResult& result, real dt)
    {
    #pragma omp parallel
        {
            real heatCapRatio_1 = data.heatCapRatio-1;
    #pragma omp for
            for(int i = 0;i<data.number;++i)
            {
                data.internalEnergy[i] += 0.5 * (result.internalEnergyDif[i] + result.pastInternalEnergyDif[i]) * dt;//修正オイラー法
                data.pressure[i] = data.internalEnergy[i]*(heatCapRatio_1*data.density[i]);
            }
        }
    }

};