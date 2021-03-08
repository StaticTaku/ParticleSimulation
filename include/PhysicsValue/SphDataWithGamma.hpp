#pragma once 
#include <PhysicsValue/SphCoreData.hpp>

struct SphDataWithGamma:public SphCoreDataWithFixedH
{
    real heatCapRatio;
    SphDataWithGamma(std::function<double(double*,double*,double)> _kernelW, std::function<void(double*,double*,double,double*)> _gradKernelW, real _h, real _heatCapRatio, int num):SphCoreDataWithFixedH(_kernelW,_gradKernelW,_h,num),heatCapRatio(_heatCapRatio){}
};