#pragma once

using real = double;

#if defined(N_BodyPlummerModel)
    using real = double;
    constexpr short DIM = 3;
    constexpr int N = 4000;
    constexpr real ZERO = 1.0e-5;
    constexpr real softParam = 0.01;
    constexpr int seed = 10;
    constexpr real R = 1;
    constexpr real Mass = 1;
    constexpr real max_speed =  0.5;
    constexpr real soft_Param = 0.01;
    constexpr real dt = 0.01;
    
#elif defined(ShockTube)
    using real = double;
    constexpr short DIM = 1;
    constexpr int N = 2000;
    constexpr real ZERO = 1.0e-8;
    constexpr real eta = 1e-5; //0割りの回避のための定数
    constexpr real Ccfl = 0.1;

    constexpr int number = N;
    constexpr real h = 0.05;
    constexpr real heatCapRatio = 1.4;
    constexpr real alpha = 1.0; //人口粘性の強さを決める係数
    constexpr real cbeta = 2*alpha; //人口粘性の強さを決める係数
    real dt = 0;

    constexpr double dens_L = 1.0;
    constexpr double vel_L = 0;
    constexpr double pres_L = 1.0;

    constexpr double dens_R = 0.125;
    constexpr double vel_R = 0.0;
    constexpr double pres_R = 0.1;

    constexpr real length = 12;
#elif defined(WaterDam)
    using real = double;
    constexpr short DIM = 2;
    constexpr int N = 50000;
    constexpr real ZERO = 1.0e-8;
    constexpr real eta = 1e-5; //0割りの回避のための定数
    constexpr real Ccfl = 0.1;
    constexpr int number = N;
    constexpr real h = 0.01;
    constexpr real heatCapRatio = 1.4;
    constexpr real alpha = 1.0; //人口粘性の強さを決める係数
    constexpr real cbeta = 2*alpha; //人口粘性の強さを決める係数
    real dt = 0.001;
#endif
