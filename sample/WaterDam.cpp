#define WaterDam
#include <Settings.hpp>

#include <iostream>
#include <fstream>
#include <KernelFunction/CubicSpline.hpp>
#include <Calculator/CalcFluidValue.hpp>
#include <Calculator/AddGravity.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <Advancer/WaterFluidValueAdvancer.hpp>
#include <Advancer/VelocityVerlet_Boundary.hpp>
#include <InitValueSetter/SetWaterDam.hpp>
#include <TimeStepDecider/GetTimeStep.hpp>
#include <omp.h>

using namespace std;

ofstream fs("Data/WaterDam/fluidValue.csv");
ofstream fs2("Data/WaterDam/density.csv");

int main()
{   
    omp_set_num_threads(12);
    SphDataWithGamma data(KernelW,GradKernelW,h,heatCapRatio,N);
    WaterFluidValueAdvancer fluidAdvancer;
    FluidResult result;
    CalcFluidValue CalcFluid(alpha,cbeta);

    int boundaryCondition = SetWaterDam(30,data);
    fluidAdvancer.UpdateDensityAndPressure(data,result);
    result.ClearAccelAndIEDif(data.number);
    CalcFluid(data,result);
    AddGravity(data.number, result);
    
    std::swap(data.accel,result.accel);

    double t = 0;


    for(int i = 0;i<1000;++i)
    {
        /*TO DO:水のsph粒子の位置をセーブデータへ記録*/
        for(int i = boundaryCondition;i<data.number;++i)
            fs << data.position[i][0] << "," << data.position[i][1] << "\n";
        for(int i = boundaryCondition;i<data.number;++i)
            fs2 << data.density[i] << "," <<  data.pressure[i] <<  "\n";

        fs << "\n\n";
        fs2 << "\n\n";

        cout << i << "\n";

        fluidAdvancer.UpdatePosition(data,result,boundaryCondition,dt);
        fluidAdvancer.UpdateDensityAndPressure(data,result);

        result.ClearAccelAndIEDif(data.number);
        CalcFluid(data,result);//内部エネルギーの時間微分と加速度を求める
        if(std::isnan(result.accel[i][0]))
        {
            std::cout << "122221212122";
        }
        AddGravity(data.number, result);

        fluidAdvancer.UpdateVelocity(data,result,dt);
        std::swap(data.accel,result.accel);
    }
    cout << t << "\n";  
}