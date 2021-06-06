#pragma once
#include <Settings.hpp>
#include <CalculatedResult/CoreResult.hpp>

inline void AddGravity(int number, CoreResult& result)
{
    for(int i = 0;i < number;++i)
        result.accel[i][1] -= 9.8;
}