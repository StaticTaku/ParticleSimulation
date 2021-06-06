#pragma once
#include <Settings.hpp>

struct CoreData
{
    int number;
    real position[N][DIM];
    real velocity[N][DIM];
    real actualVelocity[N][DIM];
    real _accel[N][DIM];
    real (*accel)[DIM];
    real mass[N];
    real potential[N];

    CoreData(int num):number(num),accel(_accel) {}

    void ClearPotential()
    {
    #pragma omp parallel
        {
    #pragma omp for
            for(int i = 0;i<number;++i)
                potential[i] = 0;
        }
    }

    double GetKineticEnergy() const
    {
        double _KE = 0;
        int i;
    #pragma omp parallel for reduction(+: _KE)
        for (i = 0; i < number; ++i)
        {
            double _velocity[3] = { velocity[i][0], velocity[i][1] ,velocity[i][2] };
            _KE += 0.5 * mass[i] * (_velocity[0] * _velocity[0] + _velocity[1] * _velocity[1] + _velocity[2] * _velocity[2]);
        }

        return _KE;
    }

    double GetPotentialEnergy() const
    {
        double _PE = 0;
        double _PE_Object = 0;
        int i;

    #pragma omp parallel for reduction(+: _PE) reduction(+: _PE_Object)
        for(i = 0;i< number;++i)
        {
            _PE += potential[i];
        }

        return 0.5*_PE + _PE_Object;
    }
};