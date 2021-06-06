#pragma once 
#include <PhysicsValue/CoreData.hpp>
#include <functional>


struct SphCoreDataWithFixedH:public CoreData
{
    real density[N];
    real pressure[N];
    real internalEnergy[N];
    real actualInternalEnergy[N];
    real h;

	std::function<double(double*,double*,double)> kernelW; //カーネル関数
	std::function<void(double*,double*,double,double*)> gradKernelW;//カーネル関数の勾配

    SphCoreDataWithFixedH(std::function<double(double*,double*,double)> _kernelW, std::function<void(double*,double*,double,double*)> _gradKernelW, real _h, int num):kernelW(_kernelW), gradKernelW(_gradKernelW), h(_h),CoreData(num) {}
};