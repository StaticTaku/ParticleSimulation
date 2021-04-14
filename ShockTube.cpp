#include <iostream>
#include <fstream>
#include <Settings.hpp>
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

constexpr int number = N;
constexpr real h = 0.05;
constexpr real heatCapRatio = 1.4;
constexpr real alpha = 1.0; //人口粘性の強さを決める係数
constexpr real cbeta = 2*alpha; //人口粘性の強さを決める係数
real dt = 0;

constexpr double dens_L = 1.0;
constexpr double vel_L = 0.75;
constexpr double pres_L = 1.0;

constexpr double dens_R = 0.125;
constexpr double vel_R = 0.0;
constexpr double pres_R = 0.1;  
constexpr real length = 12;

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

    while(t < 5)
    {
        for(int j = 0;j<data.number;++j)
        {
            fs << data.position[j][0] << "," << data.pressure[j] << "," << data.density[j] << "," << data.internalEnergy[j] << "," << data.velocity[j][0] << "\n";
        }
        fs << "\n\n";
        
        cout << t << "\n";

        dt = GetNextTimeStep(data);
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