#pragma once
#include <cmath>

namespace ObservedData
{
    namespace M31
    {
        /*
            XY平面にM31の円盤、Z軸をM31の円盤の回転軸とし、原点にM31の中心を据えた座標系の座標を
            地上からM31を観測するように、x軸を東向き、y軸を北向き、z軸を地上から空に向かう方向に設定した座標系に
            おける座標に変換する回転行列
        */
        double M31ToObs[3][3] = {{-0.1796539,   0.60181502,  0.77816653},
                                { 0.13537892,  0.79863551, -0.58639054},
                                {-0.97437006,  0.        , -0.22495105}};

        /*
            地上からM31を観測するように、x軸を東向き、y軸を北向き、z軸を地上から空に向かう方向に設定した座標系の座標を
            XY平面にM31の円盤、Z軸をM31の円盤の回転軸とし、原点にM31の中心を据えた座標系に
            おける座標に変換する回転行列
        */
        double ObsToM31[3][3] = {{-0.1796539,   0.13537892, -0.97437006},
                                { 0.60181502,  0.79863551,  0.        },
                                { 0.77816653, -0.58639054, -0.22495105}};
        
        //M31によって作られたシェル構造の、地球から観測した際の位置
        inline void set_observed_shell_position(double ShellPos[][4])
        {
        /* edge position of the observed shells */
        /* xi (deg.), eta(deg.), R(deg.), sigma(deg.) */

        ShellPos[ 0][0] =  1.710514e+00;  ShellPos[ 0][1] =  8.964539e-02;  ShellPos[ 0][2] = 1.712861e+00;  ShellPos[ 0][3] = 9.411902e-02;
        ShellPos[ 1][0] =  1.749714e+00;  ShellPos[ 1][1] =  1.530762e-01;  ShellPos[ 1][2] = 1.756398e+00;  ShellPos[ 1][3] = 1.128078e-01;
        ShellPos[ 2][0] =  1.787052e+00;  ShellPos[ 2][1] =  2.194214e-01;  ShellPos[ 2][2] = 1.800472e+00;  ShellPos[ 2][3] = 1.311400e-01;
        ShellPos[ 3][0] =  1.821983e+00;  ShellPos[ 3][1] =  2.885742e-01;  ShellPos[ 3][2] = 1.844694e+00;  ShellPos[ 3][3] = 1.486299e-01;
        ShellPos[ 4][0] =  1.853934e+00;  ShellPos[ 4][1] =  3.603516e-01;  ShellPos[ 4][2] = 1.888630e+00;  ShellPos[ 4][3] = 1.647893e-01;
        ShellPos[ 5][0] =  1.882361e+00;  ShellPos[ 5][1] =  4.345703e-01;  ShellPos[ 5][2] = 1.931873e+00;  ShellPos[ 5][3] = 1.790995e-01;
        ShellPos[ 6][0] =  1.906714e+00;  ShellPos[ 6][1] =  5.109253e-01;  ShellPos[ 6][2] = 1.973981e+00;  ShellPos[ 6][3] = 1.910850e-01;
        ShellPos[ 7][0] =  1.926482e+00;  ShellPos[ 7][1] =  5.889893e-01;  ShellPos[ 7][2] = 2.014508e+00;  ShellPos[ 7][3] = 2.003279e-01;
        ShellPos[ 8][0] =  1.941278e+00;  ShellPos[ 8][1] =  6.684570e-01;  ShellPos[ 8][2] = 2.053143e+00;  ShellPos[ 8][3] = 2.065743e-01;
        ShellPos[ 9][0] =  1.950777e+00;  ShellPos[ 9][1] =  7.488403e-01;  ShellPos[ 9][2] = 2.089568e+00;  ShellPos[ 9][3] = 2.096310e-01;
        ShellPos[10][0] =  1.954785e+00;  ShellPos[10][1] =  8.297729e-01;  ShellPos[10][2] = 2.123607e+00;  ShellPos[10][3] = 2.095312e-01;
        ShellPos[11][0] =  1.953215e+00;  ShellPos[11][1] =  9.107666e-01;  ShellPos[11][2] = 2.155120e+00;  ShellPos[11][3] = 2.063373e-01;
        ShellPos[12][0] =  1.946084e+00;  ShellPos[12][1] =  9.915771e-01;  ShellPos[12][2] = 2.184140e+00;  ShellPos[12][3] = 2.002526e-01;
        ShellPos[13][0] =  1.933512e+00;  ShellPos[13][1] =  1.071777e+00;  ShellPos[13][2] = 2.210696e+00;  ShellPos[13][3] = 1.915678e-01;
        ShellPos[14][0] =  1.915702e+00;  ShellPos[14][1] =  1.151062e+00;  ShellPos[14][2] = 2.234918e+00;  ShellPos[14][3] = 1.805377e-01;
        ShellPos[15][0] =  1.892900e+00;  ShellPos[15][1] =  1.229248e+00;  ShellPos[15][2] = 2.257016e+00;  ShellPos[15][3] = 1.675611e-01;
        ShellPos[16][0] =  1.865396e+00;  ShellPos[16][1] =  1.306152e+00;  ShellPos[16][2] = 2.277222e+00;  ShellPos[16][3] = 1.528622e-01;
        ShellPos[17][0] =  1.833526e+00;  ShellPos[17][1] =  1.381592e+00;  ShellPos[17][2] = 2.295782e+00;  ShellPos[17][3] = 1.368008e-01;
        ShellPos[18][0] =  1.797571e+00;  ShellPos[18][1] =  1.455566e+00;  ShellPos[18][2] = 2.312993e+00;  ShellPos[18][3] = 1.196085e-01;
        ShellPos[19][0] = -1.239271e+00;  ShellPos[19][1] =  1.530396e+00;  ShellPos[19][2] = 1.969239e+00;  ShellPos[19][3] = 8.656696e-02;
        ShellPos[20][0] = -1.276441e+00;  ShellPos[20][1] =  1.468384e+00;  ShellPos[20][2] = 1.945624e+00;  ShellPos[20][3] = 8.914600e-02;
        ShellPos[21][0] = -1.311935e+00;  ShellPos[21][1] =  1.406982e+00;  ShellPos[21][2] = 1.923740e+00;  ShellPos[21][3] = 9.178571e-02;
        ShellPos[22][0] = -1.345926e+00;  ShellPos[22][1] =  1.345947e+00;  ShellPos[22][2] = 1.903442e+00;  ShellPos[22][3] = 9.447306e-02;
        ShellPos[23][0] = -1.378372e+00;  ShellPos[23][1] =  1.285400e+00;  ShellPos[23][2] = 1.884718e+00;  ShellPos[23][3] = 9.718903e-02;
        ShellPos[24][0] = -1.409365e+00;  ShellPos[24][1] =  1.225220e+00;  ShellPos[24][2] = 1.867478e+00;  ShellPos[24][3] = 9.997197e-02;
        ShellPos[25][0] = -1.438988e+00;  ShellPos[25][1] =  1.165283e+00;  ShellPos[25][2] = 1.851640e+00;  ShellPos[25][3] = 1.028142e-01;
        ShellPos[26][0] = -1.467256e+00;  ShellPos[26][1] =  1.105591e+00;  ShellPos[26][2] = 1.837164e+00;  ShellPos[26][3] = 1.057104e-01;
        ShellPos[27][0] = -1.494130e+00;  ShellPos[27][1] =  1.046265e+00;  ShellPos[27][2] = 1.824032e+00;  ShellPos[27][3] = 1.086781e-01;
        ShellPos[28][0] = -1.519786e+00;  ShellPos[28][1] =  9.869385e-01;  ShellPos[28][2] = 1.812125e+00;  ShellPos[28][3] = 1.116862e-01;
        ShellPos[29][0] = -1.544128e+00;  ShellPos[29][1] =  9.278564e-01;  ShellPos[29][2] = 1.801457e+00;  ShellPos[29][3] = 1.147752e-01;
        ShellPos[30][0] = -1.567262e+00;  ShellPos[30][1] =  8.687744e-01;  ShellPos[30][2] = 1.791948e+00;  ShellPos[30][3] = 1.179140e-01;
        ShellPos[31][0] = -1.589190e+00;  ShellPos[31][1] =  8.096924e-01;  ShellPos[31][2] = 1.783571e+00;  ShellPos[31][3] = 1.211250e-01;
        ShellPos[32][0] = -1.609868e+00;  ShellPos[32][1] =  7.507324e-01;  ShellPos[32][2] = 1.776309e+00;  ShellPos[32][3] = 1.243936e-01;
        ShellPos[33][0] = -1.629384e+00;  ShellPos[33][1] =  6.916504e-01;  ShellPos[33][2] = 1.770105e+00;  ShellPos[33][3] = 1.277300e-01;
        ShellPos[34][0] = -1.647712e+00;  ShellPos[34][1] =  6.325073e-01;  ShellPos[34][2] = 1.764942e+00;  ShellPos[34][3] = 1.311226e-01;
        ShellPos[35][0] = -1.664863e+00;  ShellPos[35][1] =  5.732422e-01;  ShellPos[35][2] = 1.760788e+00;  ShellPos[35][3] = 1.345793e-01;
        ShellPos[36][0] = -1.680832e+00;  ShellPos[36][1] =  5.138550e-01;  ShellPos[36][2] = 1.757624e+00;  ShellPos[36][3] = 1.380930e-01;
        ShellPos[37][0] = -1.695610e+00;  ShellPos[37][1] =  4.543457e-01;  ShellPos[37][2] = 1.755427e+00;  ShellPos[37][3] = 1.416683e-01;
        ShellPos[38][0] = -1.709217e+00;  ShellPos[38][1] =  3.945923e-01;  ShellPos[38][2] = 1.754174e+00;  ShellPos[38][3] = 1.453000e-01;
        ShellPos[39][0] = -1.721626e+00;  ShellPos[39][1] =  3.346558e-01;  ShellPos[39][2] = 1.753850e+00;  ShellPos[39][3] = 1.489860e-01;
        ShellPos[40][0] = -1.732840e+00;  ShellPos[40][1] =  2.744446e-01;  ShellPos[40][2] = 1.754439e+00;  ShellPos[40][3] = 1.527233e-01;
        ShellPos[41][0] = -1.742839e+00;  ShellPos[41][1] =  2.139893e-01;  ShellPos[41][2] = 1.755927e+00;  ShellPos[41][3] = 1.565089e-01;
        ShellPos[42][0] = -1.751613e+00;  ShellPos[42][1] =  1.532440e-01;  ShellPos[42][2] = 1.758303e+00;  ShellPos[42][3] = 1.603414e-01;
        ShellPos[43][0] = -1.759145e+00;  ShellPos[43][1] =  9.219360e-02;  ShellPos[43][2] = 1.761559e+00;  ShellPos[43][3] = 1.642155e-01;
        ShellPos[44][0] = -1.765419e+00;  ShellPos[44][1] =  3.081512e-02;  ShellPos[44][2] = 1.765687e+00;  ShellPos[44][3] = 1.681279e-01;
        ShellPos[45][0] = -1.770413e+00;  ShellPos[45][1] = -3.090286e-02;  ShellPos[45][2] = 1.770683e+00;  ShellPos[45][3] = 1.720740e-01;
        ShellPos[46][0] = -1.774109e+00;  ShellPos[46][1] = -9.297943e-02;  ShellPos[46][2] = 1.776543e+00;  ShellPos[46][3] = 1.760502e-01;
        ShellPos[47][0] = -1.776481e+00;  ShellPos[47][1] = -1.554260e-01;  ShellPos[47][2] = 1.783267e+00;  ShellPos[47][3] = 1.800517e-01;
        }
    }
}

namespace TheoreticalData
{
    /*
        At t = 0
        0<x                    x<0
        pressure = 1.0         pressure = 0.1
        density = 1.0          density = 0.125
        velocity = 0           velocity = 0

        heatCapacityRatio = 1.4

        時刻　t　=　0.14154の時の流体の速度、圧力、密度を返す
    */
   namespace SODsShockTube1D
   {
       constexpr double heatCapRatio = 1.4;

       constexpr double t = 0.14154;

       constexpr double pres_L = 1.0, dens_L = 1.0, vel_L = 0;
       const double c_L = sqrt(heatCapRatio*(pres_L/dens_L));

       constexpr double pres_R = 0.1, dens_R = 0.125, vel_R = 0;
       const double c_R = sqrt(heatCapRatio*(pres_R/dens_R));

       inline double GetVelocity(const double x)
       {
            if(x <= -0.1674724)
                return vel_L;
            if(-0.1674724 < x && x <= -0.0099465)
            {
                return 2.0/(heatCapRatio+1) * (x/t + c_L + (heatCapRatio-1)/2.0*vel_L);
            }
            if(-0.0099465 < x && x <= 0.2480001)
            {
                return 0.9274521;
            }
            if(0.2480001 < x)
            {
                return vel_R;
            }
        }

        inline double GetPressure(const double x)
        {
            if(x <= -0.1674724)
                return pres_L;
            if(-0.1674724 < x && x <= -0.0099465)
            {
                return pres_L*pow(abs(x/t - GetVelocity(x))/c_L,2*heatCapRatio/(heatCapRatio-1));
            }
            if(-0.0099465 < x && x <= 0.2480001)
            {
                return 0.30313;
            }
            if(0.2480001 < x)
            {
                return pres_R;
            }
        }

        inline double GetDensity(const double x)
        {
            if(x <= -0.1674724)
                return dens_L;
            if(-0.1674724 < x && x <= -0.0099465)
            {
                return pow(GetPressure(x)*dens_L/pres_L,1/heatCapRatio);
            }
            if(-0.0099465 < x && x <= 0.1312716)
            {
                return 0.4263192;
            }
            if(0.1312716 < x && x <= 0.2480001)
            {
                return 0.2655736;
            }
            if(0.2480001 < x)
            {
                return vel_R;
            }
        }

        inline double GetInternalEnergy(const double x)
        {
            return 0;
        }
    }
}
