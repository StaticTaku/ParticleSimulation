#pragma once
#include <Settings.hpp>
#include <utility>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/CpiData.hpp>
#include <PhysicsValue/SphCoreData.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <CalculatedResult/CpiCoreResult.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <Advancer/LeapFrog.hpp>

#if defined(SphMethod) 
namespace SPH
{
    #if defined(VelocityVerletMethod)
    namespace VelocityVerlet
    {
        class FluidValueAdvancer:public VelocityVerlet
        {
        public:
            FluidValueAdvancer():actualVelocity(_actualVelocity),actualInternalEnergy(_actualInternalEnergy),pastInternalEnergyDif(_pastInternalEnergyDif),pastAccel(_pastAccel) {}
            real _actualVelocity[N][DIM];
            real (*actualVelocity)[DIM];

            real _actualInternalEnergy[N];
            real (*actualInternalEnergy);

            real _pastInternalEnergyDif[N];
            real (*pastInternalEnergyDif);

            real _pastAccel[N][DIM];
            real (*pastAccel)[DIM];

            void SetActual(SphCoreDataWithFixedH& data)
            {
                for(int i = 0;i<data.number;++i)
                {
                    for(int d = 0;d<DIM;++d)
                        actualVelocity[i][d] = data.velocity[i][d];
                    actualInternalEnergy[i] = data.internalEnergy[i];
                }
            }

            void UpdatePosition(SphCoreDataWithFixedH& data, const CoreResult& result, real dt)
            {
            #pragma omp parallel
                {
                    int i, j;
            #pragma omp for
                    for (i = 0; i < data.number; ++i)
                    {
                        for (j = 0; j < DIM; ++j)
                            data.position[i][j] += (actualVelocity[i][j] + 0.5 * data.accel[i][j] * dt )* dt;
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
                            data.velocity[i][j] = actualVelocity[i][j] + data.accel[i][j] * dt;
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
                            actualVelocity[i][j] += 0.5 * (data.accel[i][j] + pastAccel[i][j]) * dt;
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
            #pragma omp for
                    for(int i = 0;i<data.number;++i)
                    {
                        data.internalEnergy[i] = actualInternalEnergy[i] + pastInternalEnergyDif[i] * dt;
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
                        actualInternalEnergy[i] += 0.5 * (result.internalEnergyDif[i] + pastInternalEnergyDif[i]) * dt;//修正オイラー法
                        data.pressure[i] = actualInternalEnergy[i]*(heatCapRatio_1*data.density[i]);
                    }
                }
            }

        };
    }

    #elif defind(LeapFrogMethod)
    namespace LeapFrog
    {
        class FluidValueAdvancer:public LeapFrog
        {
        public:
            real _actualVelocity[N][DIM];
            real (*actualVelocity)[DIM];//updated by dt/2

            real _actualInternalEnergy[N];
            real (*actualInternalEnergy);//updated by dt/2

            FluidValueAdvancer():actualVelocity(_actualVelocity),actualInternalEnergy(_actualInternalEnergy) {}

            void SetActual(const SphCoreDataWithFixedH& data)
            {
                for(int i = 0;i<data.number;++i)
                {
                    for(int d = 0;d<DIM;++d)
                        actualVelocity[i][d] = data.velocity[i][d];
                    actualInternalEnergy[i] = data.internalEnergy[i]; 
                }
            }

            void PredictVelocity(SphCoreDataWithFixedH& data, CoreResult& result, real dt)
            {
            #pragma omp parallel
                {
                    int i,j;
            #pragma omp for
                    for(i = 0;i<data.number;++i)
                    {
                        for(j = 0; j<DIM;++j)
                            data.velocity[i][j] = actualVelocity[i][j] + data.accel[i][j]*dt;
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
                            actualVelocity[i][j] += data.accel[i][j]*dt/2.0;
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
                            data.position[i][j] += actualVelocity[i][j]*dt;
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
                        data.internalEnergy[i] = actualInternalEnergy[i] + result.internalEnergyDif[i]*dt;
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
                        actualInternalEnergy[i] +=  result.internalEnergyDif[i]*dt/2;
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
    #endif
}
#endif

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