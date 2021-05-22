#define N_BodyPlummerModel
#include <Settings.hpp>

#include <iostream>
#include <fstream>
#include <Calculator/CalcSelfGrav.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <InitValueSetter/SetPlummerModel.hpp>
#include <omp.h>

using namespace std;

int main()
{
    cout << DIM;
    CoreData data(N); //現時点でのposition,velocity,accel,massを持った構造体
    VelocityVerlet coreAdvancer; //現時点でのposition,velocityを求め、dataに格納するクラス
    CoreResult result; //accelの計算結果を格納する構造体
    CalcSelfGrav calcSelfGrav; //現時点でのaccelfを求め、resultに格納するクラス

    SetPlummerModel(seed,R,Mass,data); //初期条件を設定する。position,velocity,massを決定する

    data.ClearPotential();
    result.ClearAccel(data.number); 

    calcSelfGrav(softParam,data,result); //現時点でのaccelを求める
    std::swap(data.accel,result.accel);//初期時刻におけるsph粒子のaccelをdata構造体に設定

    real t = 0;

    for(int i = 0;i<4000;++i)
    {
        cout << data.GetPotentialEnergy() + data.GetKineticEnergy() << "\n";
        
        cout << i << "\n";

        t += dt;//時間をdt進める

        //この時点でdata構造体に格納されているデータはすべてdtだけ前のデータ

        coreAdvancer.UpdatePosition(data,result,dt); //現在のpositionを求める

        data.ClearPotential();
        result.ClearAccel(data.number); 
        calcSelfGrav(softParam,data,result); //現在の加速度(result.accelに格納)を求める

        coreAdvancer.UpdateVelocity(data,result,dt); //現在のvelocityを求める
        std::swap(data.accel,result.accel); //現時点での加速度(result.accel)をdata.accelに移動

        //この時点でdata構造体に格納されているデータはすべて最新の状態に
        
    }  
}