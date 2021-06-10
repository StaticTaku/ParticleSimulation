#define ShockTube1D
#define GodunovMethod
#include <Settings.hpp>

#include <iostream>
#include <fstream>
#include <PhysicsValue/CoreGridData.hpp>
#include <Advancer/FluidValueAdvancer.hpp>
#include <InitValueSetter/SetShockTube.hpp>
#include <TimeStepDecider/GetTimeStep.hpp>
#include <omp.h>


using namespace std;

ofstream fs("Data/ShockTube/Godunov/1D/fluidValue.csv");

int main()
{
    omp_set_num_threads(12);
    CoreGridData data(N);
    Godunov::FluidValueAdvancer advancer;
    real dx = Grid::SetShockTube1D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,length,data);

    real t = 0;
    constexpr real endTime = 0.14154;
    bool isNotReached = true;

    while(isNotReached)
    {
        real dt = GetNextTimeStep(dx,data);
        if(t + dt > endTime)
        {
            dt = endTime - t;
            isNotReached = false;
        }
        t += dt;//時間をdt進める

        advancer.UpdateDensity_Momentum_Energy(dx,dt,data);

        cout << t << "\n";
    }

    for(int i = 0;i<data.number;++i)
        fs << data.position[i][0] << "," << (data.energy[i]-data.momentum[i][0]*data.momentum[i][0]/(2*data.density[i]))*(heatCapRatio-1) << "," << data.density[i] << "," << data.momentum[i][0]/data.density[i] << "\n";
}