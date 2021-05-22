#pragma once

using real = double;

#if defined(N_BodyPlummerModel)
    using real = double;
    constexpr short DIM = 3; //次元
    constexpr int N = 4000; //粒子数
    constexpr real ZERO = 1.0e-5; //絶対値がこれ以下の値は0とみなす
    constexpr real softParam = 0.01; //万有引力のゼロ割を防ぐソフトニングパラメータ
    constexpr int seed = 10; //乱数生成時のシード
    constexpr real R = 1; //プラマーモデルの半径
    constexpr real Mass = 1; //プラマーモデルの総質量
    constexpr real dt = 0.01; //時間の刻み幅(可変にしたければ適宜メインファイルで)
    
#elif defined(ShockTube) //SPH法で実装
    using real = double;
    constexpr short DIM = 1; //次元
    constexpr int N = 2000; //粒子数
    constexpr real ZERO = 1.0e-8; //絶対値がこれ以下の値は0とみなす
    constexpr real eta = 1e-5; //0割りの回避のための定数
    constexpr real Ccfl = 0.1; 

    constexpr int number = N;
    constexpr real h = 0.05; //固定長
    constexpr real heatCapRatio = 1.4; //気体の比熱
    constexpr real alpha = 1.0; //人口粘性の強さを決める係数
    constexpr real cbeta = 2*alpha; //人口粘性の強さを決める係数
    real dt = 0;

    constexpr double dens_L = 1.0; //衝撃波管の左側気体の密度
    constexpr double vel_L = 0; //衝撃波管の左側気体の速度
    constexpr double pres_L = 1.0; //衝撃波管の左側気体の圧力

    constexpr double dens_R = 0.125; //衝撃波管の右側気体の密度
    constexpr double vel_R = 0.0; //衝撃波管の右側気体の速度
    constexpr double pres_R = 0.1; //衝撃波管の右側気体の圧力

    constexpr real length = 12; //衝撃波管の全体の長さ

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
