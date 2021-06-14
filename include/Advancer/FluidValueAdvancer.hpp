#pragma once
#include <Settings.hpp>
#include <utility>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/SphCoreData.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <PhysicsValue/RectangleGridData.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <Advancer/LeapFrog.hpp>
#include <Utility/RiemanSolver.hpp>

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

    #elif defined(LeapFrogMethod)
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
#elif defined(GodunovMethod)
namespace Godunov
{
    class FluidValueAdvancer
    {
    private:
        RiemanSolver riemanSolver;
    public:
        void UpdateDensity_Momentum_Energy(const real* dr, const real dt,RectangleGridData& data)
        {
        #pragma omp parallel
            {
                real dens_vel_pres[3];
                real Flux_i_minusHalf[3];
                real Flux_i_plusHalf[3];
        #pragma omp for
                for(int i = 0;i<data.AllNum();++i)
                {
                    riemanSolver.GetShockTubeAnswerAtZero(data.density[i-1],data.momentum[i-1][0]/data.density[i-1],(data.energy[i-1]-data.momentum[i-1][0]*data.momentum[i-1][0]/(2*data.density[i-1]))*(heatCapRatio-1)
                                                        ,data.density[i],data.momentum[i][0]/data.density[i],(data.energy[i]-data.momentum[i][0]*data.momentum[i][0]/(2*data.density[i]))*(heatCapRatio-1)
                                                        ,dens_vel_pres);

                    Flux_i_minusHalf[0] = dens_vel_pres[0]*dens_vel_pres[1];                
                    Flux_i_minusHalf[1] = dens_vel_pres[2]+Flux_i_minusHalf[0]*dens_vel_pres[1];
                    Flux_i_minusHalf[2] = (heatCapRatio/(heatCapRatio-1)*dens_vel_pres[2]+Flux_i_minusHalf[0]*dens_vel_pres[1]/2)*dens_vel_pres[1];

                    riemanSolver.GetShockTubeAnswerAtZero(data.density[i],data.momentum[i][0]/data.density[i],(data.energy[i]-data.momentum[i][0]*data.momentum[i][0]/(2*data.density[i]))*(heatCapRatio-1)
                                                        ,data.density[i+1],data.momentum[i+1][0]/data.density[i+1],(data.energy[i+1]-data.momentum[i+1][0]*data.momentum[i+1][0]/(2*data.density[i+1]))*(heatCapRatio-1)
                                                        ,dens_vel_pres);

                    Flux_i_plusHalf[0] = dens_vel_pres[0]*dens_vel_pres[1];                
                    Flux_i_plusHalf[1] = dens_vel_pres[2]+Flux_i_plusHalf[0]*dens_vel_pres[1];
                    Flux_i_plusHalf[2] = (heatCapRatio/(heatCapRatio-1)*dens_vel_pres[2]+Flux_i_plusHalf[0]*dens_vel_pres[1]/2)*dens_vel_pres[1];

                    data.density[i] += -dt/dx * (Flux_i_plusHalf[0] - Flux_i_minusHalf[0]);
                    data.momentum[i][0] += -dt/dx * (Flux_i_plusHalf[1] - Flux_i_minusHalf[1]);
                    data.energy[i] += -dt/dx * (Flux_i_plusHalf[2] - Flux_i_minusHalf[2]);


                }
            }
        }
    };
}
#endif