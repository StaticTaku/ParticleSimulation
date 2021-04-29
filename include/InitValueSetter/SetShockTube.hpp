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
    }

    real savedDens_L = data.density[data.number/4];
    real savedDens_R = data.density[3*data.number/4];

    for(int i = 0;i<data.number/2;++i)
    {
        data.density[i] /= savedDens_L;
        data.mass[i] /= savedDens_L;
    }

    for(int i = data.number/2;i<data.number;++i)
    {
        data.density[i] = data.density[i]/savedDens_R * dens_R;
        data.mass[i] = data.mass[i]/savedDens_R * dens_R;
    }
    for(int i = 0;i<data.number;++i)
        data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);               
}

void SetPos(int* state, int num, int dim, double length, SphDataWithGamma& data)
{
    if(dim < DIM)
    {
        for(int i = 0;i<num;++i)
        {
            state[dim] = i;
            SetPos(state, num, dim+1, length, data);
        }
    }
    else
    {
        int id = 0;
        for(int d = 0;d<DIM;++d)
            id += std::pow(num,d)*state[d];
        for(int d = 0;d<DIM;++d)
            data.position[id][d] = -length/2+length*double(state[d])/(num-1);

        return;
    }
}

void SetValue(int* state, int num,int dim, double length, SphDataWithGamma& data, double dens_L, double pres_L, double* vel_L, double dens_R, double pres_R, double* vel_R)
{
    if(dim < DIM-1)
    {
        for(int i = 0;i<num;++i)
        {
            state[dim+1] = i;
            SetValue(state,num,dim+1,length,data,dens_L,pres_L,vel_L,dens_R,pres_R,vel_R);
        }
    }
    else
    {
        int tempID = 0;
        for(int d = 1;d<DIM;++d)
            tempID += std::pow(num,d)*state[d];

        for(int i = 0;i<num/2;++i)
        {
            int id = tempID + i;
            for(int d = 0;d<DIM;++d)
                data.velocity[id][d] = vel_L[d];

            data.pressure[id] = pres_L;
            data.mass[id] = 0.01*dens_L;
        }

        for(int i = num/2;i<num;++i)
        {
            int id = tempID + i;
            for(int d = 0;d<DIM;++d)
                data.velocity[id][d] = vel_R[d];

            data.pressure[id] = pres_R;
            data.mass[id] = 0.01*dens_R;
        }

    }
}

inline void SetShockTube(double dens_L, double pres_L, double* vel_L, double dens_R, double pres_R, double* vel_R, double length, SphDataWithGamma& data)
{
    int num = std::pow(data.number,1.0/DIM);
    int state[DIM];

    data.number = 1;
    for(int d = 0;d<DIM;++d)
        data.number *= num; 

    SetPos(state,num,0,length,data);
    
    for(int d = 0;d<DIM;++d)
        state[d] = 0;
    
    SetValue(state,num,0,length,data,dens_L,pres_L,vel_L,dens_R,pres_R,vel_R);

    for(int i = 0;i<data.number;++i)
    {
        data.density[i] = 0;
        for(int j = 0;j<data.number;++j)
            data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);

        data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);           
    } 
}

