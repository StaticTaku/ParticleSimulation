#define _ShockTube
#include <Settings.hpp>

#include <iostream>
#include <fstream>

#include <Calculator/CalcFluidValue.hpp>
#include <PhysicsValue/CpiData.hpp>
#include <CalculatedResult/CpiCoreResult.hpp>
#include <Advancer/FluidValueAdvancer.hpp>
#include <InitValueSetter/SetShockTube.hpp>
#include <Utility/Data.hpp>
#include <omp.h>

using namespace std;

int main()
{
      
    omp_set_num_threads(12);
    CpiData data(heatCapRatio,N); //現時点でのposition,velocity,accel,mass,density,pressure,internalEnergy,h,heatCapRatioを持った構造体
    CIP::FluidValueAdvancer fluidAdvancer; //現時点でのdensiy,pressure,internalEnergy,position,velocityを求め、dataに格納するクラス
    CpiCoreResult result; //accel,internalEnergyDifの計算結果を格納する構造体
    CPI::CalcFluidValue calcFluid; //現時点でのaccel,internalEnergyDifを求め、resultに格納するクラス

    ofstream fs("Data/ShockTube/CIP/1D/fluidValue.csv");

    //const real dx = Grid::SetShockTube1D(dens_L,pres_L,vel_L,dens_R,pres_R,vel_R,length,data);
    const real dx = Grid::SetShockTube1D(dens_L,pres_L,3,0,0,0,length,data);
    cout << dx << "\n";
    for(int i = 0;i<data.number;++i)
    {
        data.density_xDif[i] = 0;
        data.pressure_xDif[i] = 0;
        data.velocity_xDif[i][0] = 0;
    }

    real dt = 0.001;

    real t = 0;
    constexpr real endTime = 1;

    for(int _ = 0;_ < endTime/dt;++_)
    {
        cout << t << "\n";

        t += dt;//時間をdt進める

        calcFluid.CalcIntervalValue(dx,dt,data,result); //現時点でのaccel,internalEnergyDifを求める
        for(int i = 0;i<data.number;++i)
        {
            data.density[i] = result.density[i];
            data.velocity[i][0] = result.velocity[i][0];
            data.pressure[i] = result.pressure[i];
        }
        //fluidAdvancer.UpdateVelocityDensityPressure(data,result,dt,dx);
        //fluidAdvancer.UpdateVelocityDensityPressure_xDif(data,result,dt,dx);

    }

    dt = endTime-t;
    
    t += dt;//t=endTime
    cout << t << "\n";
    //この時点でdata構造体に格納されているデータはすべてdtだけ前のデータ

    calcFluid.CalcIntervalValue(dx,dt,data,result); //現時点でのaccel,internalEnergyDifを求める
    for(int i = 0;i<data.number;++i)
    {
        data.density[i] = result.density[i];
        data.velocity[i][0] = result.velocity[i][0];
        data.pressure[i] = result.pressure[i];
    }
    for(int i = 0;i<data.number;++i)
        fs << data.position[i][0] << "," << data.pressure[i] << "," << data.density[i] << "," << data.velocity[i][0] <<  "\n";

}