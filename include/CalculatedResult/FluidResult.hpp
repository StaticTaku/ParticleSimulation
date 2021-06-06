#pragma once
#include <Settings.hpp>
#include <CalculatedResult/CoreResult.hpp>


struct FluidResult:public CoreResult
{
    real _internalEnergyDif[N];
    real predictedInternalEnergy[N];
    real _pastInternalEnergyDif[N];
    real max_ui[N];

    real *internalEnergyDif;
    real *pastInternalEnergyDif;

    FluidResult():internalEnergyDif(_internalEnergyDif),pastInternalEnergyDif(_pastInternalEnergyDif) 
    {
        for(int i = 0;i<N;++i)
            max_ui[i] = 0;
    }

    void ClearAccelAndIEDif(int number)
    {
        #pragma omp parallel
        {
        #pragma omp for
            for(int i = 0;i<number;++i)
            {
                for(int d = 0;d<DIM;++d)
                {
                    accel[i][d] = 0;
                }
                std::swap(internalEnergyDif,pastInternalEnergyDif);
                internalEnergyDif[i] = 0;
            }
        }
    }
};