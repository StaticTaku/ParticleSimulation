#pragma once
#define N_BodyPlummerModel

#if defined(N_BodyPlummerModel)
    using real = double;
    constexpr short DIM = 3;
    constexpr int N = 4000;
    constexpr real ZERO = 1.0e-5;
    constexpr real softParam = 0.01;
#elif defined(ShockTube)
    using real = double;
    constexpr short DIM = 1;
    constexpr int N = 2000;
    constexpr real ZERO = 1.0e-8;
    constexpr real eta = 1e-5; //0割りの回避のための定数
    constexpr real Ccfl = 0.1;
#endif
