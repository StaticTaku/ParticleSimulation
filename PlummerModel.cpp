#include <iostream>
#include <fstream>
#include <Settings.hpp>
#include <Calculator/CalcSelfGrav.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <Advancer/VelocityVerlet.hpp>
#include <InitValueSetter/SetPlummerModel.hpp>
#include <omp.h>

using namespace std;
constexpr int seed = 10;
constexpr real R = 1;
constexpr real Mass = 1;
constexpr real max_speed =  0.5;
constexpr real soft_Param = 0.01;
constexpr real dt = 0.01;

int main()
{
    cout << DIM;
    CoreData data(N);
    VelocityVerlet coreAdvancer;
    CoreResult result;
    CalcSelfGrav calcSelfGrav;

    SetPlummerModel(seed,R,Mass,data);

    data.ClearPotential();
    result.ClearAccel(data.number);

    calcSelfGrav(softParam,data,result);
    std::swap(data.accel,result.accel);

    for(int i = 0;i<4000;++i)
    {
        cout << data.GetPotentialEnergy() + data.GetKineticEnergy() << "\n";
        
        cout << i << "\n";
        coreAdvancer.UpdatePosition(data,result,dt);

        data.ClearPotential();
        result.ClearAccel(data.number);
        calcSelfGrav(softParam,data,result);

        coreAdvancer.UpdateVelocity(data,result,dt);
        std::swap(data.accel,result.accel);
    }  
}