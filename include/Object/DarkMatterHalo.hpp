#pragma once
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <cmath>


//NFW profile
class DarkMatterHalo
{
private:
	double c;
	double r_s;
	double haloPosition[3];
	mutable double _potential;

public:
	DarkMatterHalo(double _c,double _scaleRadius,double x, double y, double z):c(_c),r_s(_scaleRadius) 
	{
		haloPosition[0] = x;
		haloPosition[1] = y;
		haloPosition[2] = z;
	}

	double M(double r) const
	{
		double r_ = r/r_s;
		return c * ((log(1 + r_) - r_ / (1 + r_)));
	}

	double _GetPotential(double r) const
	{
		return -c/r*log(1+r/r_s);
	}
    
	void DarkMatterHalo::Action(CoreData& data, CoreResult& result) const
    {
        _potential = 0;
        #pragma omp parallel
        {
            double eps2 = softParam * softParam;
            double ri_haloPos[3];
            double length;
            double sqrt_length;
            int i, j;
        #pragma omp for reduction(+: _potential)
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    ri_haloPos[j] = data.position[i][j]-haloPosition[j];
                }
        
                length = ri_haloPos[0]*ri_haloPos[0]+ri_haloPos[1]*ri_haloPos[1]+ri_haloPos[2]*ri_haloPos[2]+eps2;
                sqrt_length = sqrt(length);

                for (j = 0; j < 3; ++j)
                    result.accel[i][j] += -M(sqrt_length)*ri_haloPos[j]/(length*sqrt_length);

                data.potential[i] += data.mass[i]*_GetPotential(sqrt_length);
            }

        #pragma omp barrier
        }
    }

	double GetPosition(int i) const 
	{
		return haloPosition[i];
	}

};