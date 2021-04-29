#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <cmath>

class Bulge
{
private:
	double bulgeMass;
    double bulgeScaleRadius;
    double bulgePosition[3];
	mutable double _potential;

public:
	Bulge(double _mass,double _scaleRadius,double x, double y, double z):bulgeMass(_mass),bulgeScaleRadius(_scaleRadius) 
	{
		bulgePosition[0] = x;
		bulgePosition[1] = y;
		bulgePosition[2] = z;
	}

	double M(double r) const
	{
		double r_ = (1.0/(1.0+bulgeScaleRadius/r));
		return r_*r_*bulgeMass;
	}

	double _GetPotential(double r) const
	{
		return -bulgeMass/(r+bulgeScaleRadius);
	}
    
    void Action(CoreData& data,CoreResult& result) const
    {
        #pragma omp parallel
        {
            double eps2 = softParam * softParam;
            double ri_bulgePos[3];
            double length;
            double sqrt_length;
            int i, j;
        #pragma omp for reduction(+: _potential)
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    ri_bulgePos[j] = data.position[i][j]-bulgePosition[j];
                }
        
                length = ri_bulgePos[0]*ri_bulgePos[0]+ri_bulgePos[1]*ri_bulgePos[1]+ri_bulgePos[2]*ri_bulgePos[2]+eps2;
                sqrt_length = sqrt(length);

                for (j = 0; j < 3; ++j)
                    result.accel[i][j] += -M(sqrt_length)*ri_bulgePos[j]/(length*sqrt_length);


                data.potential[i] += data.mass[i]*_GetPotential(sqrt_length);
            }

        #pragma omp barrier
        }
    }

	double GetPosition(int i) const
	{
		return bulgePosition[i];
	}
};