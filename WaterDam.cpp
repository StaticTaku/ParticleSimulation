#include <iostream>
#include <fstream>
#include <Settings.hpp>
#include <KernelFunction/CubicSpline.hpp>
#include <Calculator/CalcFluidValue.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <Advancer/WaterFluidValueAdvancer.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <InitValueSetter/SetRandomValue.hpp>
#include <TimeStepDecider/GetTimeStep.hpp>
#include <omp.h>
using namespace std;

constexpr int number = N;
constexpr real h = 0.05;
constexpr real heatCapRatio = 1.4;
constexpr real alpha = 0.5; //人口粘性の強さを決める係数
constexpr real cbeta = 2*alpha; //人口粘性の強さを決める係数
real dt = 0;

ofstream fs("Data/FluidWithSelfGravity/fluidValue.csv");

int main()
{   
    omp_set_num_threads(12);
    SphDataWithGamma data(KernelW,GradKernelW,h,heatCapRatio,N);
    WaterFluidValueAdvancer fluidAdvancer;
    FluidResult result;
    CalcFluidValue CalcFluid(alpha,cbeta);

    /*TO DO:ダムの初期条件設定*/

    result.ClearAccelAndIEDif(data.number);
    CalcFluid(data,result);
    std::swap(data.accel,result.accel);

    double t = 0;

    while(t < 10)
    {
        /*TO DO:水のsph粒子の位置をセーブデータへ記録*/

        cout << t << "\n";

        dt = GetNextTimeStep(data);
        t += dt;

        fluidAdvancer.UpdatePosition(data,result,dt);
        fluidAdvancer.UpdateDensityAndPressure(data,result);

        std::swap(result.internalEnergyDif,result.pastInternalEnergyDif);

        result.ClearAccelAndIEDif(data.number);
        CalcFluid(data,result);//内部エネルギーの時間微分と加速度を求める

        fluidAdvancer.UpdateInternalEnergy(data,result,dt);
        fluidAdvancer.UpdateVelocity(data,result,dt);
        std::swap(data.accel,result.accel);
    }
    cout << t << "\n";  
}