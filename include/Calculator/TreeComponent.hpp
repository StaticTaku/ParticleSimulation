#pragma once
#include <Settings.hpp>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>

enum class NodeType
{
    BODY = 1,
    CELL = 2,
};

namespace BodysData
{
    extern double (*position)[3];
    extern double (*velocity)[3];
    extern double* mass;
    extern double (*accel)[3];
    extern double (*result)[3];
}

void SetBodysPointer(CoreData& data, CoreResult& result)
{
    BodysData::position = data.position;
    BodysData::velocity = data.velocity;
    BodysData::mass = data.mass;
    BodysData::accel = data.accel;
    BodysData::result = result.accel;
}

class CommonNode
{
private:
    NodeType type;
    CommonNode* next = nullptr;

public:
    CommonNode(NodeType _type):type(_type) {}

    void SetNext(CommonNode* _next) { next = _next; }
    CommonNode* GetNext() { return next; }
    NodeType GetType() { return type; }
    virtual double GetPosition(int i) = 0;
    virtual void SetPosition(double* p) = 0;
    virtual void SetPosition(int i,double x) = 0;
    virtual double GetMass() = 0;
};

class Body:public CommonNode
{
private:
    int bodyID;
public:
    Body():CommonNode(NodeType::BODY) {}

    void SetBodyID(int _bodyID)
    {
        bodyID = _bodyID;
    }

    int GetBodyID()
    {
        return bodyID;
    }

    double GetPosition(int i) override
    {
        return BodysData::position[bodyID][i];
    }

    void SetPosition(double* p) override
    {
        for(int i = 0;i<3;++i)
            BodysData::position[bodyID][i] = p[i];
    }

    void SetPosition(int i,double x) override
    {
        BodysData::position[bodyID][i] = x; 
    }

    double GetMass() override
    {
        return BodysData::mass[bodyID];
    }

    double GetVelocity(int i)
    {
        return BodysData::velocity[bodyID][i];
    }

    double GetAccel(int i)
    {
        return BodysData::accel[bodyID][i];
    }
};

class Cell:public CommonNode
{
private:
    double rcrit2;
    CommonNode* more = nullptr;
    CommonNode* subp[8];
    double mass;
    double position[3];

public:
    Cell():CommonNode(NodeType::CELL) {}

    double GetPosition(int i) override
    {
        return position[i];
    }

    void SetPosition(double* p) override
    {
        for(int i = 0;i<3;++i)
            position[i] = p[i];
    }

    void SetPosition(int i,double x) override
    {
        position[i] = x; 
    }

    double GetMass() override
    {
        return mass;
    }

    void SetMass(double _mass) 
    {
        mass = _mass;
    }

    void AddMass(double _mass)
    {
        mass += _mass;
    }

    void SetSubp(int index,CommonNode* p)
    {
        subp[index] = p;
    }

    double GetCriticalRadius2()
    {
        return rcrit2;
    }

    void SetCriticalRadius2(double _rcrit2)
    {
        rcrit2 = _rcrit2;
    }

    void SetMore(CommonNode* p)
    {
        more = p;
    }

    CommonNode* GetSubp(int index)
    {
        return subp[index];
    }

    void ResetSubp()
    {
        for(int i = 0;i<8;++i)
            subp[i] = nullptr;
    }

    CommonNode* GetMore()
    {
        return more;
    }

};
