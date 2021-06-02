#pragma once
#include <PhysicsValue/CoreGridData.hpp>

struct CpiData:public CoreGridData
{
   real velocity_xDif[N][DIM];
   real density_xDif[N];
   real pressure_xDif[N];
   real gamma;

   CpiData(real _gamma, int number):CoreGridData(number),gamma(_gamma) {}
};