#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <PhysicsValue/CpiData.hpp>
#include <CalculatedResult/CpiCoreResult.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <cmath>

#if defined(ShockTube1D) || defined(ShockTube2D)
namespace SPH
{
    class CalcFluidValue
    {
    public:
        real alpha = 0.5; //人口粘性の強さを決める係数
        real beta = 2*alpha; //人口粘性の強さを決める係数

        CalcFluidValue(real _alpha, real _beta):alpha(_alpha),beta(_beta) {}
        void operator()(const SphDataWithGamma& data, FluidResult& result)
        {
        #pragma omp parallel
            {
                real accel[DIM];
                real pi,pj;
                real hi,hj;
                real densi,densj;
                real pos_i[DIM];
                real vel_i[DIM];
                real pos_j[DIM];
                real vel_j[DIM];
                real grad_Wij_hi[DIM];
                real grad_Wij_hj[DIM];

                real _internalEnergyDif;
                real pi_ij;

                int i,j,k;
        #pragma omp for
                for(i = 0;i<data.number;++i)
                {    
                    pi = data.pressure[i];
                    densi = data.density[i];
                    hi = data.h;
                    for(k=0;k<DIM;++k)
                    {
                        pos_i[k] = data.position[i][k];
                        vel_i[k] = data.velocity[i][k];
                    }

                    for(k = 0;k<DIM;++k)
                        accel[k] = 0;
                    
                    for(j = 0;j<data.number;++j)
                    {
                        pj = data.pressure[j];
                        densj = data.density[j];
                        hj = data.h;

                        for(k=0;k<DIM;++k)
                        {
                            pos_j[k] = data.position[j][k];
                            vel_j[k] = data.velocity[j][k];
                        }
                            
                        data.gradKernelW(pos_i,pos_j,hi,grad_Wij_hi);
                        data.gradKernelW(pos_i,pos_j,hj,grad_Wij_hj);

                        _internalEnergyDif = 0;
                        //オイラー方程式から導かれるsph粒子の加速度と内部エネルギーの微分
                        for(k = 0;k<DIM;++k)
                        {
                            accel[k] += (-data.mass[j]*(pi/(densi*densi)*grad_Wij_hi[k]+pj/(densj*densj)*grad_Wij_hj[k]));
                            _internalEnergyDif += (data.velocity[i][k]-data.velocity[j][k])*grad_Wij_hi[k];
                        }
                        result.internalEnergyDif[i] += pi/densi*data.mass[j]*_internalEnergyDif;

                        
                        //人口粘性による追加の加速度と内部エネルギーの微分の値
                        pi_ij = PI_ij(pos_i,pos_j,vel_i,vel_j,densi,densj,hi,hj,pi,pj,data.heatCapRatio);

                        if(!(abs(pi_ij) <= ZERO))
                        {
                            _internalEnergyDif = 0;
                            for(k = 0;k<DIM;++k)
                            {
                                _internalEnergyDif += (vel_i[k]-vel_j[k])*(grad_Wij_hi[k]+grad_Wij_hj[k]);
                                accel[k] += (-data.mass[j]*pi_ij*(grad_Wij_hi[k]+grad_Wij_hj[k])/2.0);
                            }
                            result.internalEnergyDif[i] += 0.5*data.mass[j]*pi_ij*0.5*_internalEnergyDif;
                        }
                    }
                    for(k = 0;k<DIM;++k)
                        result.accel[i][k] += accel[k];
                }

        #pragma omp barrier
            }    
        }

        double PI_ij(double* ri,double* rj,double* vi,double* vj,double densi, double densj, double hi, double hj, double presi, double presj, double heatCapRatio)
        {
            real rij[DIM];
            real vij[DIM];
            double rij_dot_vij = 0;

            for(int d = 0;d<DIM;++d)
            {
                rij[d] = (ri[d]-rj[d]);
                vij[d] = (vi[d]-vj[d]);
                rij_dot_vij += rij[d]*vij[d];
            }
            if(rij_dot_vij >= 0)
            {
                return 0;
            }

            double hij = (hi + hj)/2.0;
            double densij = (densi + densj)/2.0;
            double cij = (sqrt(heatCapRatio*presi/densi)+sqrt(heatCapRatio*presj/densj))/2.0;
            
            double _uij = 0;
            for(int d = 0;d<DIM;++d)
            {
                _uij += (rij[d]*rij[d]);
            }
            
            double uij = hij*rij_dot_vij/(_uij+eta*hij*hij);
            
            return (-alpha*cij*uij+beta*uij*uij)/densij;
        }   
    };
}

#elif defined(_ShockTube)
namespace CPI
{
    class CalcFluidValue
    {
    public:
        
        void CalcIntervalValue(real dx, real dt, const CpiData& data, CpiCoreResult& result)
        {
        #pragma omp parallel
            {
                real accel[DIM];
                real pi,pj;
                real densi,densj;
                real pos_i[DIM];
                real vel_i[DIM];
                int i,j,k;
                real udt;

                real a,b,c,d;
        #pragma omp for
                for(i = 1;i<data.number-1;++i)
                {    
                    pi = data.pressure[i];
                    densi = data.density[i];
                    for(k=0;k<DIM;++k)
                    {
                        pos_i[k] = data.position[i][k];
                        vel_i[k] = data.velocity[i][k];
                    }

                    udt = vel_i[i]*dt;
                    
                   if(vel_i[i] >= 0)
                   {
                        a = a_plus(i,dx,data.density,data.density_xDif);
                        b = b_plus(i,dx,data.density,data.density_xDif);
                        c = data.density_xDif[i];
                        d = data.density[i];

                        result.density[i] = F(a,b,c,d,-udt);
                        result.density_xDif[i] = F_dx(a,b,c,-udt);
                    
                        a = a_plus(i,dx,data.pressure,data.pressure_xDif);
                        b = b_plus(i,dx,data.pressure,data.pressure_xDif);
                        c = data.pressure_xDif[i];
                        d = data.pressure[i];

                        result.pressure[i] = F(a,b,c,d,-udt);
                        result.pressure_xDif[i] = F_dx(a,b,c,-udt);

                        a = a_plus(i,dx,data.velocity,data.velocity_xDif);
                        b = b_plus(i,dx,data.velocity,data.velocity_xDif);
                        c = data.velocity_xDif[i][0];
                        d = data.velocity[i][0];

                        result.velocity[i][0] = F(a,b,c,d,-udt);
                        result.velocity_xDif[i][0] = F_dx(a,b,c,-udt);
                   }
                    
                    if(vel_i[i] < 0)
                   {
                        a = a_minus(i,dx,data.density,data.density_xDif);
                        b = b_minus(i,dx,data.density,data.density_xDif);
                        c = data.density_xDif[i];
                        d = data.density[i];

                        result.density[i] = F(a,b,c,d,-udt);
                        result.density_xDif[i] = F_dx(a,b,c,-udt);
                    
                        a = a_minus(i,dx,data.pressure,data.pressure_xDif);
                        b = b_minus(i,dx,data.pressure,data.pressure_xDif);
                        c = data.pressure_xDif[i];
                        d = data.pressure[i];

                        result.pressure[i] = F(a,b,c,d,-udt);
                        result.pressure_xDif[i] = F_dx(a,b,c,-udt);

                        a = a_minus(i,dx,data.velocity,data.velocity_xDif);
                        b = b_minus(i,dx,data.velocity,data.velocity_xDif);
                        c = data.velocity_xDif[i][0];
                        d = data.velocity[i][0];

                        result.velocity[i][0] = F(a,b,c,d,-udt);
                        result.velocity_xDif[i][0] = F_dx(a,b,c,-udt);
                   }



                }

        #pragma omp barrier
            }    
        }

        real a_plus(int i, real dx, const real* f, const real* f_dx)
        {
            return (f_dx[i] + f_dx[i-1])/(dx*dx) - 2*(f[i]-f[i-1])/(dx*dx*dx);
        }

        real a_plus(int i, real dx, const real (*f)[DIM], const real (*f_dx)[DIM])
        {
            return (f_dx[i][0] + f_dx[i-1][0])/(dx*dx) - 2*(f[i][0]-f[i-1][0])/(dx*dx*dx);
        }

        real a_minus(int i, real dx, const real* f, const real* f_dx)
        {
            return (f_dx[i] + f_dx[i+1])/(dx*dx) + 2*(f[i]-f[i+1])/(dx*dx*dx);
        }

        real a_minus(int i, real dx, const real (*f)[DIM], const real (*f_dx)[DIM])
        {
            return (f_dx[i][0] + f_dx[i+1][0])/(dx*dx) + 2*(f[i][0]-f[i+1][0])/(dx*dx*dx);
        }

        real b_plus(int i, real dx, const real* f, const real* f_dx)
        {
            return 3*(f[i-1] - f[i])/(dx*dx) + 2*(f_dx[i] + f_dx[i-1])/dx;
        }

        real b_plus(int i, real dx, const real (*f)[DIM], const real (*f_dx)[DIM])
        {
            return 3*(f[i-1][0] - f[i][0])/(dx*dx) + 2*(f_dx[i][0] + f_dx[i-1][0])/dx;
        }

        real b_minus(int i, real dx, const real* f, const real* f_dx)
        {
            return 3*(f[i+1] - f[i])/(dx*dx) - 2*(f_dx[i] + f_dx[i+1])/dx;
        }

        real b_minus(int i, real dx, const real (*f)[DIM], const real (*f_dx)[DIM])
        {
            return 3*(f[i+1][0] - f[i][0])/(dx*dx) - 2*(f_dx[i][0] + f_dx[i+1][0])/dx;
        }

        real F(real a,real b, real c, real d, real x)
        {
            return a*x*x*x + b*x*x + c*x + d;
        }

        real F_dx(real a, real b, real c, real x)
        {
            return 3*a*x*x + 2*b*x + c;
        }
    };
}
#endif