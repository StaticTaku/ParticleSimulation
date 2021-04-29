#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <cmath>

class CalcSelfGrav
{
public:
	void operator() (real _softParam, CoreData& data, CoreResult& result)
	{
	#pragma omp parallel
		{
			double _accel[DIM], ri[DIM], dr[DIM];
			double ri2, mri5_2,sqrt_ri2;
			double eps2 = _softParam * _softParam;
			int i, j, k;
			double _ri2;
	#pragma omp for
			for (i = 0; i < data.number; ++i)
			{
				data.potential[i] = 0;

				for (j = 0; j < DIM; ++j)
				{
					_accel[j] = 0;
					ri[j] = data.position[i][j];
				}

				for (j = 0; j < data.number; ++j)
				{
					if (i == j)
						continue;

					_ri2 = 0;
					for (k = 0; k < DIM; ++k)
					{
						dr[k] = data.position[j][k] - ri[k];
						_ri2 += dr[k]*dr[k];
					}

					
					ri2 = 1.0 / (_ri2 + eps2);
					sqrt_ri2 = sqrt(ri2);

					mri5_2 = data.mass[j] * ri2 * sqrt_ri2;

					for (k = 0; k < DIM; ++k)
						_accel[k] += mri5_2 * dr[k];

					data.potential[i] += -data.mass[i]*data.mass[j]*sqrt_ri2;
				}

				for (j = 0; j < DIM; ++j)
					result.accel[i][j] += _accel[j];
			}

	#pragma omp barrier
		}
	}
};