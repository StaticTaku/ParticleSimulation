#define ShockTube2D
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
#include <Utility/Data.hpp>
#include <omp.h>

using namespace std;

#if defined(ShockTube)
    ofstream fs("Data/ShockTube/1D/fluidValue.csv");
#elif defined(ShockTube2D)
    ofstream fs("Data/ShockTube/2D/fluidValue.csv");
#endif

int main()
{   
    omp_set_num_threads(12);
    SphDataWithGamma data(KernelW,GradKernelW,h,heatCapRatio,N); //現時点でのposition,velocity,accel,mass,density,pressure,internalEnergy,h,heatCapRatioを持った構造体
    FluidValueAdvancer fluidAdvancer; //現時点でのdensiy,pressure,internalEnergy,position,velocityを求め、dataに格納するクラス
    FluidResult result; //accel,internalEnergyDifの計算結果を格納する構造体
    CalcFluidValue CalcFluid(alpha,cbeta); //現時点でのaccel,internalEnergyDifを求め、resultに格納するクラス

    #if defined(ShockTube)
        SetShockTube1D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,length,data);//初期条件を設定する。position,velocity,mass,density,pressure,internalEnergyを決定する
    #elif defined(ShockTube2D)
        SetShockTube2D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,x_length,y_length,x_num,y_num,data);//初期条件を設定する。position,velocity,mass,density,pressure,internalEnergyを決定する
    #endif

    result.ClearAccelAndIEDif(data.number);
    CalcFluid(data,result); //現時点でのaccel,internalEnergyDifを求める
    std::swap(data.accel,result.accel);//初期時刻におけるsph粒子のaccelをdata構造体に設定

    real t = 0;
    constexpr real endTime = 0.14154;
    dt = 0.00002;
    for(int _ = 0;_ < endTime/dt;++_)
    {
        cout << t << "\n";

        t += dt;//時間をdt進める

        //この時点でdata構造体に格納されているデータはすべてdtだけ前のデータ

        fluidAdvancer.UpdatePosition(data,result,dt); //現在のpositionを求める
        fluidAdvancer.UpdateDensity(data,result); //現在のdensityを求める
        fluidAdvancer.PredictInternalEnergyAndPressure(data,result,dt); //暫定的に現在のpressureを求める

        std::swap(result.internalEnergyDif,result.pastInternalEnergyDif);//dt前の内部エネルギーの時間微分が格納されている(internalEnergyDifに格納されている)をpastInternalEnergyDifへ移動

        result.ClearAccelAndIEDif(data.number);
        CalcFluid(data,result);//現在の内部エネルギーの時間微分(result.internalEnergyDif)と加速度(result.accel)を求める

        fluidAdvancer.UpdateInternalEnergyAndPressure(data,result,dt);//現在のinternalEnergyとpressureを求める
        fluidAdvancer.UpdateVelocity(data,result,dt); //現在のvelocityを求める
        std::swap(data.accel,result.accel); //現時点での加速度(result.accel)をdata.accelに移動

        //この時点でdata構造体に格納されている物理量はすべて最新の状態に
    }

    dt = endTime-t;
    
    t += dt;//t=endTime
    cout << t << "\n";
    //この時点でdata構造体に格納されているデータはすべてdtだけ前のデータ

    fluidAdvancer.UpdatePosition(data,result,dt); //現在のpositionを求める
    fluidAdvancer.UpdateDensity(data,result); //現在のdensityを求める
    fluidAdvancer.PredictInternalEnergyAndPressure(data,result,dt); //暫定的に現在のpressureを求める

    std::swap(result.internalEnergyDif,result.pastInternalEnergyDif);//dt前の内部エネルギーの時間微分が格納されている(internalEnergyDifに格納されている)をpastInternalEnergyDifへ移動

    result.ClearAccelAndIEDif(data.number);
    CalcFluid(data,result);//現在の内部エネルギーの時間微分(result.internalEnergyDif)と加速度(result.accel)を求める

    fluidAdvancer.UpdateInternalEnergyAndPressure(data,result,dt);//現在のinternalEnergyとpressureを求める
    fluidAdvancer.UpdateVelocity(data,result,dt); //現在のvelocityを求める
    std::swap(data.accel,result.accel); //現時点での加速度(result.accel)をdata.accelに移動

    //この時点でdata構造体に格納されている物理量はすべて最新の状態に

    #if defined(Shocktube)
        for(int i = 0;i<data.number;++i)
            fs << data.position[i][0] << "," << data.pressure[i] << "," << data.density[i] << "," << data.velocity[i][0] << "," << data.internalEnergy[i] << "\n";
    #elif defined(ShockTube2D)
        for(int x = 0;x<x_num;++x)
            fs << data.position[(y_num/2)*x_num+x][0] << "," << data.pressure[(y_num/2)*x_num+x] << "," << data.density[(y_num/2)*x_num+x] << "," << data.velocity[(y_num/2)*x_num+x][0] << "," << data.internalEnergy[(y_num/2)*x_num+x] << "\n";
    #endif
}