#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <utility>

class VelocityVerlet_Boundary
{
public:
    void UpdatePosition(CoreData& data, const CoreResult& result, int boundaryCondition, real dt)
    {
    #pragma omp parallel
        {
            int i, j;
    #pragma omp for
            for (i = boundaryCondition; i < data.number; ++i)
            {
                for (j = 0; j < DIM; ++j)
                {
                    data.position[i][j] += (data.velocity[i][j] + 0.5 * data.accel[i][j] * dt )* dt;
                }
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
                    data.velocity[i][j] += 0.5 * (data.accel[i][j] + result.accel[i][j]) * dt;
                }
            }
        }
    #pragma omp barrier
    }
};