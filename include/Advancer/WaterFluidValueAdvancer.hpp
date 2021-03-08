#pragma once
#include <Settings.hpp>
#include <cmath>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/SphCoreData.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <Advancer/VelocityVerlet.hpp>

constexpr real B = 3000;
constexpr real dens0 = 997;
constexpr int n = 7;

class WaterFluidValueAdvancer:public VelocityVerlet
{
public:
    void UpdateDensityAndPressure(SphCoreDataWithFixedH& data, const FluidResult& result)
    {
    #pragma omp parallel
        {
            real h = data.h;
    #pragma omp for
            for(int i = 0;i<data.number;++i)
            {
                data.density[i] = 0;
                for(int j = 0;j<data.number;++j)
                {
                    data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],h);
                    data.pressure[i] = B*(std::pow(data.density[i]/dens0,n)-1);
                }
            }
        }
    }

    void UpdateInternalEnergy(SphCoreDataWithFixedH& data, const FluidResult& result, real dt)
    {
    #pragma omp parallel
        {
    #pragma omp for
            for(int i = 0;i<data.number;++i)
                data.internalEnergy[i] += 0.5 * (result.internalEnergyDif[i] + result.pastInternalEnergyDif[i]) * dt;
        }
    }
};