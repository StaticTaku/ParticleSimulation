#define ShockTube1D
#define SphMethod
#define VelocityVerletMethod
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

#if defined(ShockTube1D)
    ofstream fs("Data/ShockTube/SPH/1D/fluidValue.csv");
#elif defined(ShockTube2D)
    ofstream fs("Data/ShockTube/SPH/2D/fluidValue.csv");
#elif defined(ShockTube3D)
    ofstream fs("Data/ShockTube/SPH/3D/fluidValue.csv");
#endif

int main()
{   
    omp_set_num_threads(12);
    SphDataWithGamma data(KernelW,GradKernelW,h,heatCapRatio,N); //現時点でのposition,velocity,accel,mass,density,pressure,internalEnergy,h,heatCapRatioを持った構造体
    #if defined(LeapFrogMethod)
        SPH::LeapFrog::FluidValueAdvancer fluidAdvancer; //現時点でのdensiy,pressure,internalEnergy,position,velocityを求め、dataに格納するクラス
    #elif defined(VelocityVerletMethod)
        SPH::VelocityVerlet::FluidValueAdvancer fluidAdvancer; //この時点でdata,result構造体に格納されている物理量はすべて最新の状態に
    #endif
    FluidResult result; //internalEnergyDifの計算結果を格納する構造体
    SPH::CalcFluidValue CalcFluid; //現時点でのaccel,internalEnergyDifを求め、data.accel,result.internalEnergyDifに格納するクラス

    #if defined(ShockTube1D)
        SPH::SetShockTube1D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,length,data);//初期条件を設定する。position,velocity,mass,density,pressure,internalEnergyを決定する
    #elif defined(ShockTube2D)
        SPH::SetShockTube2D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,x_length,y_length,x_num,y_num,data);//初期条件を設定する。position,velocity,mass,density,pressure,internalEnergyを決定する
    #elif defined(ShockTube3D)
        cerr << "3次元衝撃波管の初期条件作成プログラムがまだ実装できていない\n";
        return 1;
    #endif
    
    fluidAdvancer.SetActual(data);//使用する積分スキームで必要なactualVelocity,actualInternalEnergyの初期値を設定

    result.ClearIEDif(data.number);
    data.ClearAccel();
    CalcFluid(data,result); //現時点でのaccel(data.accel),internalEnergyDif(result.internalEnergyDif)を求める

    real t = 0;
    constexpr real endTime = 0.14154;
    bool isNotReached = true;

    while(isNotReached)
    {
        dt = GetNextTimeStep(result,data);
        if(t + dt > endTime)
        {
            dt = endTime - t;
            isNotReached = false;
        }
        t += dt;//時間をdt進める

    
    #if defined(LeapFrogMethod)
        //この時点でdata構造体に格納されているデータはすべてdtだけ前のデータ
        //現在:t+dt

        fluidAdvancer.PredictVelocity(data,result,dt); //t+dtでの加速度を求めるため、t+dtでの速度を予測
        fluidAdvancer.HalfUpdateVelocity(data,result,dt); //t+dt/2での速度を求める
        fluidAdvancer.UpdatePosition(data,result,dt); //t+dtでの位置を求める。
        fluidAdvancer.PredictInternalEnergy(data,result,dt); //t+dtでの加速度を求めるため、t+dtでの内部エネルギーを予測
        fluidAdvancer.HalfUpdateInternalEnergy(data,result,dt); //t+dt/2での内部エネルギーを求める
        
        fluidAdvancer.UpdateDensity(data,result); //t+dtでの密度を求める 
        fluidAdvancer.UpdatePressure(data,result,dt); //t+dtでの加速度を求めるため、t+dtでの圧力を予測

        result.ClearIEDif(data.number);
        data.ClearAccel();
        CalcFluid(data,result);//現在(t+dt)の内部エネルギーの時間微分(result.internalEnergyDif)と加速度(data.accel)を求める

        fluidAdvancer.HalfUpdateVelocity(data,result,dt);
        fluidAdvancer.HalfUpdateInternalEnergy(data,result,dt);

        //この時点でdata,result構造体に格納されている物理量はすべて最新(t+dt)の状態に
        
    #elif defined(VelocityVerletMethod)
        //この時点でdata,result構造体に格納されているデータはすべてdtだけ前のデータ

        std::swap(fluidAdvancer.pastAccel,data.accel);
        std::swap(fluidAdvancer.pastInternalEnergyDif,result.internalEnergyDif);
        fluidAdvancer.UpdatePosition(data,result,dt); //現在のpositionを求める
        fluidAdvancer.UpdateDensity(data,result); //現在のdensityを求める
        fluidAdvancer.PredictVelocity(data,result,dt); //加速度の計算のため、オイラー法で現在の速度を予想
        fluidAdvancer.PredictInternalEnergyAndPressure(data,result,dt); //加速度の計算のため、オイラー法で現在の圧力を予想

        result.ClearIEDif(data.number);
        data.ClearAccel();
        CalcFluid(data,result);//現在の内部エネルギーの時間微分(result.internalEnergyDif)と加速度(result.accel)を求める

        fluidAdvancer.UpdateInternalEnergyAndPressure(data,result,dt);//現在のinternalEnergyとpressureを求める
        fluidAdvancer.UpdateVelocity(data,result,dt); //現在のvelocityを求める*/

        //この時点でdata,result構造体に格納されている物理量はすべて最新の状態に
    
    #endif

        cout << t << "\n";
    }

    for(int i = 0;i<data.number;++i)
        fs << data.position[i][0] << "," << data.pressure[i] << "," << data.density[i] << "," << data.velocity[i][0] << "," << data.internalEnergy[i] << "\n";
    
    #if defined(Shocktube1D)
        for(int i = 0;i<data.number;++i)
            fs << data.position[i][0] << "," << data.pressure[i] << "," << data.density[i] << "," << data.velocity[i][0] << "," << data.internalEnergy[i] << "\n";
    #elif defined(ShockTube2D)
        for(int x = 0;x<x_num;++x)
        fs << data.position[(y_num/2)*x_num+x][0] << "," << data.position[(y_num/2)*x_num+x][1] << "," << data.pressure[(y_num/2)*x_num+x] << "," << data.density[(y_num/2)*x_num+x] << "," << data.velocity[(y_num/2)*x_num+x][0] << "," << data.velocity[(y_num/2)*x_num+x][1] << "," << data.internalEnergy[(y_num/2)*x_num+x] << "\n";
    
        fs << "\n\n";

        for(int y = 0;y<y_num;++y)
            fs << data.position[(y)*x_num+x_num/4][0] << "," << data.position[(y)*x_num+x_num/4][1] << "," << data.pressure[(y)*x_num+x_num/4] << "," << data.density[(y)*x_num+x_num/4] << "," << data.velocity[(y)*x_num+x_num/4][0] << "," << data.velocity[(y)*x_num+x_num/4][1] << "," << data.internalEnergy[(y)*x_num+x_num/4] << "\n";
    #elif defined(ShockTube3D)
    #endif
}