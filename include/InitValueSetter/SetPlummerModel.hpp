#pragma once
#include <PhysicsValue/CoreData.hpp>
#include <random>
#include <omp.h>
#include <cmath>

#include <limits>   

namespace Detail
{
    double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
    {
        return curr == prev
            ? curr
            : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
}

/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
*/
double constexpr ct_sqrt(double x)
{
    return x >= 0 && x < std::numeric_limits<double>::infinity()
        ? Detail::sqrtNewtonRaphson(x, x, 0)
        : std::numeric_limits<double>::quiet_NaN();
}

constexpr double M_PI2 = 2 * M_PI;
constexpr double sq2 = ct_sqrt(2.0);

inline void SetPlummerModel(int seed,real Mass,real R,CoreData& data)
{
#pragma omp parallel
	{
		int i;
		double r,rz,v,vz;
		double X1, X2, X3, X4, X5, X6, X7;
		std::default_random_engine engine(seed + omp_get_thread_num());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double _mass = Mass/data.number;
		
#pragma omp for
		for (i = 0; i < data.number; ++i)
		{
			r = R*pow((pow(dist(engine), -2.0 / 3) - 1), -0.5);
			X1 = dist(engine);

			data.position[i][2] = (1 - 2 * dist(engine)) * r;

			rz = data.position[i][2];

			data.position[i][0] = pow((r * r - rz * rz), 0.5) * cos(M_PI2 * X1);
			data.position[i][1] = pow((r * r - rz * rz), 0.5) * sin(M_PI2 * X1);

			do
			{
				X2 = dist(engine);
				X3 = dist(engine);
			} while (X2 * X2 * pow((1 - X2 * X2), 3.5) <= 0.1 * X3);

			v = X2 * (sq2 * pow((1 + r * r), -0.25));
			X4 = dist(engine);

			data.velocity[i][2] = (1 - 2 * dist(engine)) * v;

			vz = data.velocity[i][2];

			data.velocity[i][0] = pow((v * v - vz * vz), 0.5) * cos(M_PI2 * X4);
			data.velocity[i][1] = pow((v * v - vz * vz), 0.5) * sin(M_PI2 * X4);

			data.mass[i] = _mass;
		}
	}
}


