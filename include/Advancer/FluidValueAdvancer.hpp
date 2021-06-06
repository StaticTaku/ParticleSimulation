#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/CpiData.hpp>
#include <PhysicsValue/SphCoreData.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <CalculatedResult/CpiCoreResult.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <Advancer/LeapFrog.hpp>

namespace SPH
{
    class FluidValueAdvancer_VelocityVerlet:public VelocityVerlet
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
                        data.position[i][j] += (data.actualVelocity[i][j] + 0.5 * data.accel[i][j] * dt )* dt;
                }
            }
        #pragma omp barrier
        }

        void PredictVelocity(SphCoreDataWithFixedH& data, const CoreResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i, j;
        #pragma omp for
                for (i = 0; i < data.number; ++i)
                {
                    for (j = 0; j < DIM; ++j)
                        data.velocity[i][j] = data.actualVelocity[i][j] + data.accel[i][j] * dt;
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
                        data.actualVelocity[i][j] += 0.5 * (data.accel[i][j] + result.accel[i][j]) * dt;
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
                    data.internalEnergy[i] = data.actualInternalEnergy[i] + result.internalEnergyDif[i] * dt;
                    data.pressure[i] = data.internalEnergy[i]*(heatCapRatio_1*data.density[i]);
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
                    data.actualInternalEnergy[i] += 0.5 * (result.internalEnergyDif[i] + result.pastInternalEnergyDif[i]) * dt;//修正オイラー法
                    data.pressure[i] = data.actualInternalEnergy[i]*(heatCapRatio_1*data.density[i]);
                }
            }
        }

    };

    class FluidValueAdvancer_LeapFrog:LeapFrog
    {
    public:
        void PredictVelocity(SphCoreDataWithFixedH& data, CoreResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i,j;
        #pragma omp for
                for(i = 0;i<data.number;++i)
                {
                    for(j = 0; j<DIM;++j)
                        data.velocity[i][j] = data.actualVelocity[i][j] + data.accel[i][j]*dt;
                }
            }            
        }

        void HalfUpdateVelocity(CoreData& data, const CoreResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i, j;
        #pragma omp for
                for (i = 0; i < data.number; ++i)
                {
                    for (j = 0; j < DIM; ++j)
                        data.actualVelocity[i][j] += data.accel[i][j]*dt/2.0;
                }
            }
        #pragma omp barrier
        }

        void UpdatePosition(CoreData& data, const CoreResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i, j;
        #pragma omp for
                for (i = 0; i < data.number; ++i)
                {
                    for (j = 0; j < DIM; ++j)
                        data.position[i][j] += data.actualVelocity[i][j]*dt;
                }
            }
        #pragma omp barrier
        }


        void PredictInternalEnergy(SphCoreDataWithFixedH& data, FluidResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i,j;
        #pragma omp for
                for(i = 0;i<data.number;++i)
                {
                    data.internalEnergy[i] = data.actualInternalEnergy[i] + result.internalEnergyDif[i]*dt;
                }
            }            

        }
        void HalfUpdateInternalEnergy(SphCoreDataWithFixedH& data, FluidResult& result, real dt)
        {
        #pragma omp parallel
            {
                int i,j;
        #pragma omp for
                for(i = 0;i<data.number;++i)
                {
                    data.actualInternalEnergy[i] +=  result.internalEnergyDif[i]*dt/2;
                }
            }            
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

        void UpdatePressure(SphDataWithGamma& data, const FluidResult& result, real dt)
        {
        #pragma omp parallel
            {
                real heatCapRatio_1 = data.heatCapRatio-1;
        #pragma omp for
                for(int i = 0;i<data.number;++i)
                    data.pressure[i] = data.internalEnergy[i]*(heatCapRatio_1*data.density[i]);
            }
        }


    };
}

namespace CIP
{
    class FluidValueAdvancer
    {
    public:
        void UpdateVelocityDensityPressure(CpiData& data, const CpiCoreResult& result, real dt, real dx)
        {
        #pragma omp parallel
            {
                int i, j;
        #pragma omp for
                for (i = 1; i < data.number-1; ++i)
                {
                    const real qi = VonNeumannArtificialViscosity(result.density[i],result.pressure[i],result.velocity_xDif[i][0],data.gamma,dx);
                    const real qi_1 = VonNeumannArtificialViscosity(result.density[i+1],result.pressure[i+1],result.velocity_xDif[i+1][0],data.gamma,dx);
                    for (j = 0; j < DIM; ++j)
                    {
                        data.velocity[i][j] = result.velocity[i][0] - 2/(result.density[i]+result.density[i+1]) * (data.pressure[i+1]+ qi_1 - result.pressure[i] - qi) * dt/dx;
                    }

                    data.density[i] = result.density[i] - result.density[i]*(result.velocity[i+1][0]-result.velocity[i][0])*dt/dx;
                    data.pressure[i] = result.pressure[i] - (data.gamma*result.pressure[i] + (data.gamma-1)*qi)*(result.velocity[i+1][0]-result.velocity[i][0])*dt/dx;
                }
            }
        #pragma omp barrier
        }
        
        void UpdateVelocityDensityPressure_xDif(CpiData& data, const CpiCoreResult& result, real dt, real dx)
        {
        #pragma omp parallel
            {
                int i, j;
        #pragma omp for
                for (i = 1; i < data.number-1; ++i)
                {
                    for (j = 0; j < DIM; ++j)
                    {
                        data.velocity_xDif[i][j] = result.velocity[i][0] + ((data.velocity[i+1][0] - result.velocity[i+1][0])-(data.velocity[i-1][0]-result.velocity[i-1][0]))/(2*dx) - result.velocity_xDif[i][0]*(data.velocity[i+1][0]-data.velocity[i-1][0])*dt/(2*dx);
                    }

                    data.density_xDif[i] = result.density[i] + ((data.density[i+1] - result.density[i+1])-(data.density[i-1]-result.density[i-1]))/(2*dx) - result.density_xDif[i]*(data.density[i+1]-data.density[i-1])*dt/(2*dx);
                    data.pressure_xDif[i] = result.pressure[i] + ((data.pressure[i+1] - result.pressure[i+1])-(data.pressure[i-1]-result.pressure[i-1]))/(2*dx) - result.pressure_xDif[i]*(data.pressure[i+1]-data.pressure[i-1])*dt/(2*dx);
                }
            }
        #pragma omp barrier
        }

        real VonNeumannArtificialViscosity(real dens,real pressure,real u_dx,real gamma,real lamda)
        {
            if(u_dx >= 0)
                return 0;
            
            const real Cs = std::sqrt(gamma*pressure/dens);

            return alpha*(-dens*Cs*u_dx*lamda + (gamma+1)/2.0*u_dx*u_dx*lamda*lamda);
        }

    };
}