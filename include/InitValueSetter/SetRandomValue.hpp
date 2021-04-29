#pragma once

#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <omp.h>
#include <random>

inline void SetRandomValue(int seed, real R, real Mass, real maxSpeed, CoreData& data)
{
#pragma omp parallel
	{
		int i, j;
        double _sum;
		double _value[DIM];
		double m = Mass/data.number;
		std::default_random_engine engine(seed + omp_get_thread_num());
		std::uniform_real_distribution<double> dist(0.0, 2 * R);
		std::uniform_real_distribution<double> dist2(0.0, 2 * maxSpeed);

#pragma omp for
		for (i = 0; i < data.number; ++i)
		{
			do
			{  
                _sum = 0;
				for (j = 0; j < DIM; ++j)
                {
					_value[j] = R - dist(engine);
                    _sum += _value[j] * _value[j];
                }
			} while (_sum > R * R);

			for (j = 0; j < DIM; ++j)
				data.position[i][j] = _value[j];
		}


#pragma omp for
		for (i = 0; i < data.number; ++i)
		{
			do
			{
                _sum = 0;
				for (j = 0; j < DIM; ++j)
                {
					_value[j] = maxSpeed - dist2(engine);
                    _sum += _value[j] * _value[j];
                }
			} while (_sum > maxSpeed * maxSpeed);

			for (j = 0; j < DIM; ++j)
				data.velocity[i][j] = _value[j];
		}

#pragma omp for
		for (i = 0; i < data.number; i++)
			data.mass[i] = m;

#pragma omp barrier
	}
}