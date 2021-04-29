#define ShockTube
#include <Settings.hpp>

#include <iostream>
#include <fstream>
#include <KernelFunction/CubicSpline.hpp>
#include <Calculator/CalcFluidValue.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <Advancer/FluidValueAdvancer.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <InitValueSetter/SetShockTube.hpp>
#include <TimeStepDecider/GetTimeStep.hpp>
#include <omp.h>

using namespace std;

ofstream fs("Data/ShockTube/fluidValue.csv");

int main()
{   
    omp_set_num_threads(12);
    SphDataWithGamma data(KernelW,GradKernelW,h,heatCapRatio,N);
    FluidValueAdvancer fluidAdvancer;
    FluidResult result;
    CalcFluidValue CalcFluid(alpha,cbeta);

    SetShockTube(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,length,data);

    result.ClearAccelAndIEDif(data.number);
    CalcFluid(data,result);
    std::swap(data.accel,result.accel);//初期時刻におけるsph粒子のaccelをdata構造体に設定

    double t = 0;

    while(t < 10)
    {
        for(int j = 0;j<data.number;++j)
        {
            fs << data.position[j][0] << "," << data.pressure[j] << "," << data.density[j] << "," << data.internalEnergy[j] << "," << data.velocity[j][0] << "\n";
        }
        fs << "\n\n";
        
        cout << t << "\n";

        dt = 0.002;
        t += dt;

        fluidAdvancer.UpdatePosition(data,result,dt);
        fluidAdvancer.UpdateDensity(data,result);
        fluidAdvancer.PredictInternalEnergyAndPressure(data,result,dt);

        std::swap(result.internalEnergyDif,result.pastInternalEnergyDif);//時刻tの内部エネルギーの時間微分(internalEnergyDif)をpastInternalEnergyDifへ移動

        result.ClearAccelAndIEDif(data.number);
        CalcFluid(data,result);//内部エネルギーの時間微分(result.internalEnergyDif)と加速度(result.accel)を求める

        fluidAdvancer.UpdateInternalEnergyAndPressure(data,result,dt);
        fluidAdvancer.UpdateVelocity(data,result,dt);
        std::swap(data.accel,result.accel);
    }  

    cout << t << "\n";
}