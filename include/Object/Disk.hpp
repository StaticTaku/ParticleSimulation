#pragma once
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <cmath>

class Disk//Kuzumin-Toomre Disk
{
private:
	double diskCentralSurfaceDens;
    double diskScaleRadius;
    double diskPosition[3];
    double _const;
	double _const2;
    mutable double _cache;
	mutable double _potential;

public:
	Disk(double _diskCentralSurfaceDens,double _scaleRadius,double x, double y, double z):diskCentralSurfaceDens(_diskCentralSurfaceDens),diskScaleRadius(_scaleRadius) 
	{
		diskPosition[0] = x;
		diskPosition[1] = y;
		diskPosition[2] = z;
        _const = 2*M_PI*diskScaleRadius*diskScaleRadius*diskCentralSurfaceDens;
	}

    void GetForce_rz(double length,double z,double additionalAccel[2]) const
    {
        _cache = pow(length*length+(diskScaleRadius+abs(z))*(diskScaleRadius+abs(z)),-1.5);
        additionalAccel[0] = -_const*length*_cache;
        additionalAccel[1] = z >= 0 ? -_const*_cache*(abs(z)+diskScaleRadius):_const*_cache*(abs(z)+diskScaleRadius);
    }

	double _GetPotential(double r,double z) const
	{
		_cache = diskScaleRadius+abs(z);
		return -_const*pow(r*r+_cache*_cache,-0.5);
	}

	void Action(CoreData& data, CoreResult& result) const
    {
        #pragma omp parallel
        {
            double eps2 = softParam * softParam;
            double ri_diskPos[3];
            double length;
            double length_xy;
            double sqrt_length;
            double sqrt_length_xy;
            int i, j;
            double addAccel_rz[2];
        #pragma omp for reduction(+: _potential)
            for (i = 0; i < data.number; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    ri_diskPos[j] = data.position[i][j]-diskPosition[j];
                }
        
                length = ri_diskPos[0]*ri_diskPos[0]+ri_diskPos[1]*ri_diskPos[1]+ri_diskPos[2]*ri_diskPos[2]+eps2;
                sqrt_length = sqrt(length);
                length_xy = ri_diskPos[0]*ri_diskPos[0]+ri_diskPos[1]*ri_diskPos[1]+eps2;
                sqrt_length_xy = sqrt(length_xy);
                GetForce_rz(sqrt_length_xy,ri_diskPos[2],addAccel_rz);

                for (j = 0; j < 2; ++j)
                    result.accel[i][j] += addAccel_rz[0]*ri_diskPos[j]/(sqrt_length_xy);

                result.accel[i][2] += addAccel_rz[1];

                data.potential[i] += data.mass[i]*_GetPotential(sqrt_length_xy,ri_diskPos[2]);
            }

        #pragma omp barrier
        }
    }

	double GetPosition(int i)
	{
		return diskPosition[i];
	}
};