#pragma once
#include <Settings.hpp>

struct CpiCoreResult
{
    real velocity[N][DIM];
    real density[N];
    real pressure[N];

    real velocity_xDif[N][DIM];
    real density_xDif[N];
    real pressure_xDif[N];
};