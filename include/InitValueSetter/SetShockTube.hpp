#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>

inline void SetShockTube(double dens_L,double pres_L,double vel_L,double dens_R,double pres_R,double vel_R,double length, SphDataWithGamma& data)
{
    for(int i = 0;i<data.number/2;++i)
    {
        data.pressure[i] = pres_L;
        data.velocity[i][0] = vel_L;
        data.position[i][0] = -length/2+length*double(i)/data.number;
        data.mass[i] = 0.01*dens_L;
    }

    for(int i = data.number/2;i<data.number;++i)
    {
        data.pressure[i] = pres_R;
        data.velocity[i][0] = vel_R;
        data.position[i][0] = -length/2+length*double(i)/data.number;
        data.mass[i] = 0.01*dens_R;
    }

    for(int i = 0;i<data.number;++i)
    {
        data.density[i] = 0;
        for(int j = 0;j<data.number;++j)
            data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);

        data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);           
    }    
}

inline void SetShockTube(double dens_L, double pres_L, double* vel_L, double dens_R, double pres_R, double* vel_R, double length, SphDataWithGamma& data)
{
    int num = std::pow(data.number,1.0/3);
    int last_num = data.number/num*num;

    data.number = last_num*num*num;

    for(int i = 0;i<data.number/2;++i)
    {
        data.pressure[i] = pres_L;

        for(int d = 0;d<DIM;++d)
            data.velocity[i][d] = vel_L[d];
        data.mass[i] = 0.01*dens_L;
    }

    for(int x = 0;x < num;++x)
    {
        for(int y = 0; y< num;++y)
        {
            for(int z = 0;z < last_num;++last_num)
            {
                data.position[]
            }
        }
    }
    for(int i = 0;i<num2;++i)
        data.position[i][DIM-1] = -length/2+length*double(i)/num2;

    for(int i = data.number/2;i<data.number;++i)
    {
        data.pressure[i] = pres_R;
        data.velocity[i][0] = vel_R;
        data.position[i][0] = -length/2+length*double(i)/data.number;
        data.mass[i] = 0.01*dens_R;
    }

    for(int i = 0;i<data.number;++i)
    {
        data.density[i] = 0;
        for(int j = 0;j<data.number;++j)
            data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);

        data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);           
    } 
}