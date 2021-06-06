#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <utility>

class LeapFrog
{
    /*
    initial state:v(t), x(t)
        t += dt
        get v(t+dt/2) by HalfUpdataVelocity
        get x(t+dt) by UpdatePosition
        get a(t+dt) by Calculator
        get v(t+dt) by HalfUpdateVelocity
    then get state:v(t+dt),x(t+dt)
    */
public:
    void HalfUpdateVelocity(CoreData& data, const CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                    data.velocity[i][j] += data.accel[i][j]*dt/2.0;
            }
        }
    #pragma omp barrier
    }
    
    void UpdatePosition(CoreData& data, const CoreResult& result, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                    data.position[i][j] += data.velocity[i][j]*dt;
            }
        }
    #pragma omp barrier
    }
};