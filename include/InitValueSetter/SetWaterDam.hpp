#pragma once
#include <Settings.hpp>
#include <PhysicsValue/SphCoreData.hpp>

constexpr real interval = 0.01;
inline int SetWaterDam(int layer, SphCoreDataWithFixedH& data)
{
    int a = 1/interval;
    int nextNum = 0;
    ++a;//壁を構成するsph粒子の縦の方向に積む数
    
    //左壁
    for(int w = 1;w<=a;++w)
    {
        for(int i = 0;i<layer;++i)
        {
            data.position[layer*(w-1) + i][0] = -1.0 - interval*i;
            data.position[layer*(w-1) + i][1] = 1.0 - interval*(w-1);

            data.velocity[layer*(w-1) + i][0] = 0;
            data.velocity[layer*(w-1) + i][0] = 0;
        }
    }

    nextNum += layer*a;
    

    //右壁
    for(int w = 1;w<=a;++w)
    {
        for(int i = 0;i<layer;++i)
        {
            data.position[nextNum+layer*(w-1) + i][0] = 1.0 + interval*i;
            data.position[nextNum+layer*(w-1) + i][1] = 1.0 - interval*(w-1);

            data.velocity[layer*(w-1) + i][0] = 0;
            data.velocity[layer*(w-1) + i][0] = 0;
        }
    }

    nextNum *= 2;

    a = 2/interval;
    ++a;
    //床
    for(int w = 1;w<=a;++w)
    {
        for(int i = 0;i<layer;++i)
        {
            data.position[nextNum+layer*(w-1) + i][0] = -1.0 + interval*(w-1);
            data.position[nextNum+layer*(w-1) + i][1] = -interval*i;

            data.velocity[nextNum + layer*(w-1) + i][0] = 0;
            data.velocity[nextNum + layer*(w-1) + i][0] = 0;
        }
    }

    nextNum += layer*a;

    a = 0.5/interval;
    ++a;
    //液柱
    for(int w = 1;w<=a;++w)
    {
        for(int i = 0;i<a;++i)
        {
            data.position[nextNum + layer*(w-1) + i][0] = -1.0 + i*interval;
            data.position[nextNum + layer*(w-1) + i][1] = 0.5 - (w-1) * interval;

            data.velocity[nextNum + layer*(w-1) + i][0] = 0;
            data.velocity[nextNum + layer*(w-1) + i][0] = 0;
        }
    }

    
    data.number = nextNum + layer*(a-1) + a;
    int dismiss = nextNum;

    for(int i = 0;i<data.number;++i)
        data.mass[i] = 0.1;

    return dismiss;
}