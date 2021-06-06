#pragma once
#include <Settings.hpp>
#include <CalculatedResult/CoreResult.hpp>


struct FluidResult:public CoreResult
{
    real _internalEnergyDif[N];
    real max_ui[N]; //for timestepDicider

    real *internalEnergyDif;

    FluidResult():internalEnergyDif(_internalEnergyDif) 
    {
        for(int i = 0;i<N;++i)
            max_ui[i] = 0;
    }

    void ClearIEDif(int number)
    {
        #pragma omp parallel
        {
        #pragma omp for
            for(int i = 0;i<number;++i)
            {
                internalEnergyDif[i] = 0;
                max_ui[i] = 0;
            }
        }
    }
};