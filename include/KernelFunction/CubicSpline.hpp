#pragma once
#include <Settings.hpp>
#include <cmath>


constexpr real sigma[3] = {2.0/3, 10.0 / (7.0 * M_PI), 1.0 / M_PI};
constexpr real sigma2[3] = {3.0/2, 45.0/(14 * M_PI), 9.0/(4*M_PI)};

real KernelW(const real ri[],const real rj[], const real hi)
{
    real length_ij_2 = 0;
    for(int d = 0;d<DIM;++d)
        length_ij_2 += (ri[d]-rj[d])*(ri[d]-rj[d]);
    const real length_ij =  sqrt(length_ij_2);

    const real q = length_ij / hi;
    if(0 <= q && q <= 1)
        return sigma[DIM-1]/pow(hi,DIM) * (1-1.5*q*q + 0.75*q*q*q);
    else if(1 <= q && q <= 2)
        return sigma[DIM-1]/pow(hi,DIM) * (0.25*(2-q)*(2-q)*(2-q));
    else 
        return 0;
}

void GradKernelW(const real ri[],const real rj[], const real hi, real* ans)
{
    real length_ij_2 = 0;
    for(int d = 0;d<DIM;++d)
        length_ij_2 += (ri[d]-rj[d])*(ri[d]-rj[d]);
    const real length_ij =  sqrt(length_ij_2);

    const real q = length_ij / hi;
    real c;

    if(0 <= q && q <= 1)
        c = sigma2[DIM-1]/pow(hi,DIM) *(q - 4.0/3)/(hi*hi);
    else if(1 <= q && q <= 2)
        c = sigma2[DIM-1]/pow(hi,DIM) *(-1/3.0*(2-q)*(2-q)/(q*hi*hi));
    else 
        c = 0;
    
    for(int d = 0;d<DIM;++d)
        ans[d] = c*(ri[d]-rj[d]);
}