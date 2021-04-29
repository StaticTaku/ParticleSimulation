#pragma once
#include <Object/Bulge.hpp>
#include <Object/DarkMatterHalo.hpp>
#include <Object/Disk.hpp>
#include <cmath>

class M31
{
private:
    Bulge* bulge;
    Disk* disk;
    DarkMatterHalo* dmh;
    mutable double _potential;

public:
    M31(Bulge* _bulge, Disk* _disk, DarkMatterHalo* _dmh):bulge(_bulge),disk(_disk),dmh(_dmh) {}

    void Action(CoreData& data, CoreResult& result) const
    {
        #pragma omp parallel
        {
            double eps2 = softParam * softParam;
            double ri_diskPos[3];
            double ri_bulgePos[3];
            double ri_haloPos[3];
            double length_bulge;
            double sqrt_length_bulge;
            double length_dmh;
            double sqrt_length_dmh;
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
                    ri_bulgePos[j] = data.position[i][j]-bulge->GetPosition(j);
                }

                for (j = 0; j < 3; ++j)
                {
                    ri_diskPos[j] = data.position[i][j]-disk->GetPosition(j);
                }

                for (j = 0; j < 3; ++j)
                {
                    ri_haloPos[j] = data.position[i][j]-dmh->GetPosition(j);
                }

                length_bulge = ri_bulgePos[0]*ri_bulgePos[0]+ri_bulgePos[1]*ri_bulgePos[1]+ri_bulgePos[2]*ri_bulgePos[2]+eps2;
                sqrt_length_bulge = sqrt(length_bulge);

                length_dmh = ri_haloPos[0]*ri_haloPos[0]+ri_haloPos[1]*ri_haloPos[1]+ri_haloPos[2]*ri_haloPos[2]+eps2;
                sqrt_length_dmh = sqrt(length_dmh);
        
                length = ri_diskPos[0]*ri_diskPos[0]+ri_diskPos[1]*ri_diskPos[1]+ri_diskPos[2]*ri_diskPos[2]+eps2;
                sqrt_length = sqrt(length);
                length_xy = ri_diskPos[0]*ri_diskPos[0]+ri_diskPos[1]*ri_diskPos[1]+eps2;
                sqrt_length_xy = sqrt(length_xy);
                disk->GetForce_rz(sqrt_length_xy,ri_diskPos[2],addAccel_rz);

                for (j = 0; j < 2; ++j)
                    result.accel[i][j] += (addAccel_rz[0]*ri_diskPos[j]/(sqrt_length_xy)
                                    -dmh->M(sqrt_length_dmh)*ri_haloPos[j]/(length_dmh*sqrt_length_dmh))
                                    -bulge->M(sqrt_length_bulge)*ri_bulgePos[j]/(length_bulge*sqrt_length_bulge);

                result.accel[i][2] += (addAccel_rz[1]
                                -dmh->M(sqrt_length_dmh)*ri_haloPos[2]/(length_dmh*sqrt_length_dmh) 
                                -bulge->M(sqrt_length_bulge)*ri_bulgePos[2]/(length_bulge*sqrt_length_bulge));

                data.potential[i] += data.mass[i]*(bulge->_GetPotential(sqrt_length_bulge) + disk->_GetPotential(sqrt_length_xy,ri_diskPos[2]) + dmh->_GetPotential(sqrt_length_dmh));
            }

        #pragma omp barrier
        }
    }
};