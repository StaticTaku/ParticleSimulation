#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <utility>

class VelocityVerlet
{
    /*
    initial state:v(t), x(t)
        t += dt
        get x(t+dt) by UpdatePosition
        get a(t+dt) by Calculator
        get v(t+dt) by UpdateVelocity
    then get state:v(t+dt),x(t+dt)
    */
public:
    real _pastAccel[N][DIM];
    real (*pastAccel)[DIM];

    VelocityVerlet():pastAccel(_pastAccel) {}

    void UpdatePosition(CoreData& data, const CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                    data.position[i][j] += (data.velocity[i][j] + 0.5 * data.accel[i][j] * dt )* dt;
            }
        }
    #pragma omp barrier
    }

    void UpdateVelocity(CoreData& data, CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                {
                    data.velocity[i][j] += 0.5 * (data.accel[i][j] + pastAccel[i][j]) * dt;
                }
            }
        }
    #pragma omp barrier
    }
};