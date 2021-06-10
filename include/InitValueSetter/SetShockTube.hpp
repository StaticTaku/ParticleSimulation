#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/CoreGridData.hpp>
#include <iostream>

#if defined(SphMethod)
namespace SPH
{
    inline void SetShockTube1D(real dens_L,real pres_L,real vel_L,real dens_R,real pres_R,real vel_R,real length, SphDataWithGamma& data)
    {
        for(int i = 0;i<data.number/2;++i)
        {
            data.pressure[i] = pres_L;
            data.velocity[i][0] = vel_L;
            data.position[i][0] = -length/2+length*real(i)/data.number;
            data.mass[i] = mass_coef*dens_L;
        }

        for(int i = data.number/2;i<data.number;++i)
        {
            data.pressure[i] = pres_R;
            data.velocity[i][0] = vel_R;
            data.position[i][0] = -length/2+length*real(i)/data.number;
            data.mass[i] = mass_coef*dens_R;
        }
        
        for(int i = 0;i<data.number;++i)
        {
            data.density[i] = 0;
            for(int j = 0;j<data.number;++j)
                data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);
            data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]); 
        }
    }

    inline void SetShockTube2D(real dens_L,real pres_L,real vel_L,real dens_R,real pres_R,real vel_R,real x_length,real y_length, int x_num, int y_num,SphDataWithGamma& data)
    {
        if(x_num*y_num != data.number)
        {
            std::cerr << "2次元衝撃波管問題では,衝撃波波のx軸方向のサンプリング粒子数とy軸方向のサンプリング粒子数の積を、全体の粒子数と等しくなるようにしてください。\n";
            exit(0);
        }
        for(int y = 0;y < y_num; ++y)
        {
            double y_pos =  -y_length/2+y_length*real(y)/y_num;

            for(int x = 0;x < x_num/2;++x)
            {
                int id = y*x_num + x;
                data.position[id][0] = -x_length/2+x_length*real(x)/x_num;
                data.position[id][1] = y_pos;

                data.velocity[id][0] = vel_L;
                data.velocity[id][1] = 0;

                data.pressure[id] = pres_L;
                data.mass[id] = mass_coef*dens_L;
            }
            for(int x = x_num/2;x < x_num; ++x)
            {
                int id = y*x_num + x;
                data.position[id][0] = -x_length/2+x_length*real(x)/x_num;
                data.position[id][1] = y_pos;

                data.velocity[id][0] = vel_R;
                data.velocity[id][1] = 0;

                data.pressure[id] = pres_R;
                data.mass[id] = mass_coef*dens_R;
            }
        }

        for(int i = 0;i<data.number;++i)
        {
            data.density[i] = 0;
            for(int j = 0;j<data.number;++j)
                data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);
            data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);
        }
    }

    /*inline void SetShockTube3D(real dens_L,real pres_L,real vel_L,real dens_R,real pres_R,real vel_R,real x_length,real y_length, real z_length, int x_num, int y_num, int z_num, SphDataWithGamma& data)
    {
        if(x_num*y_num*z_num != data.number)
        {
            std::cerr << "3次元衝撃波管問題では,衝撃波波のx軸方向のサンプリング粒子数とy軸方向のサンプリング粒子数とz軸方向のサンプリング粒子数の積を、全体の粒子数と等しくなるようにしてください。\n";
            exit(0);
        }
        for(int z = 0;z < z_num; ++z)
        {
            double z_pos = -z_length/2+z_length*real(z)/z_num;
            for(int y = 0;y < y_num; ++y)
            {
                double y_pos =  -y_length/2+y_length*real(y)/y_num;

                for(int x = 0;x < x_num/2;++x)
                {
                    int id = z*(x_num*y_num) + y*x_num + x;
                    data.position[id][0] = -x_length/2+x_length*real(x)/x_num;
                    data.position[id][1] = y_pos;
                    data.position[id][2] = z_pos;

                    data.velocity[id][0] = vel_L;
                    data.velocity[id][1] = 0;
                    data.velocity[id][2] = 0;

                    data.pressure[id] = pres_L;
                    data.mass[id] = 0.01*dens_L;
                }
                for(int x = x_num/2;x < x_num; ++x)
                {
                    int id = y*x_num + x;
                    data.position[id][0] = -x_length/2+x_length*real(x)/x_num;
                    data.position[id][1] = y_pos;
                    data.position[id][2] = z_pos;

                    data.velocity[id][0] = vel_R;
                    data.velocity[id][1] = 0;
                    data.velocity[id][2] = 0;

                    data.pressure[id] = pres_R;
                    data.mass[id] = 0.01*dens_R;
                }
            }
        }

        for(int i = 0;i<data.number;++i)
        {
            data.density[i] = 0;
            for(int j = 0;j<data.number;++j)
                data.density[i] += data.mass[j]*data.kernelW(data.position[i],data.position[j],data.h);
        }

        real savedDens_L = data.density[z_num/2*(x_num*y_num) + y_num/2*x_num + x_num/4];
        real savedDens_R = data.density[z_num/2*(x_num*y_num) + y_num/2*x_num + 3*x_num/4];

        for(int z = 0;z < z_num; ++z)
        {
            for(int y = 0;y < y_num; ++y)
            {
                for(int x = 0;x < x_num/2;++x)
                {
                    data.density[x] = data.density[x]/savedDens_L * dens_L;
                    data.mass[x] = data.mass[x]/savedDens_L * dens_L;
                }
                for(int x = x_num/2;x < x_num; ++x)
                {
                    data.density[x] = data.density[x]/savedDens_R * dens_R;
                    data.mass[x] = data.mass[x]/savedDens_R * dens_R;
                }
            }
        }
        for(int i = 0;i<data.number;++i)
            data.internalEnergy[i] = data.pressure[i]/((data.heatCapRatio-1)*data.density[i]);
    }*/
}

#elif defined(GodunovMethod)
namespace Grid
{
    inline real SetShockTube1D(real dens_L,real pres_L,real vel_L,real dens_R,real pres_R,real vel_R,real length,CoreGridData& data)
    {   
        real dx = length/(data.number-1);

        for(int i = 0;i<data.number;++i)
            data.position[i][0] = -length/2 + dx*i;

        for(int i = 0;i<data.number/2 + 1;++i)
        {
            data.density[i] = dens_L;
            data.momentum[i][0] = dens_L*vel_L;
            data.energy[i] = pres_L/(heatCapRatio-1) + dens_L/2 * vel_L*vel_L;
        }
        for(int i = data.number/2 + 1;i<data.number;++i)
        {
            data.density[i] = dens_R;
            data.momentum[i][0] = dens_R*vel_R;
            data.energy[i] = pres_R/(heatCapRatio-1) + dens_R/2 * vel_R*vel_R;
        }

        return dx;
    }
}
#endif