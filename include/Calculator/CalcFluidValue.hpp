#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphDataWithGamma.hpp>
#include <CalculatedResult/FluidResult.hpp>
#include <cmath>


class CalcFluidValue
{
public:
    real alpha = 0.5; //人口粘性の強さを決める係数
	real beta = 2*alpha; //人口粘性の強さを決める係数
    real dt;

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