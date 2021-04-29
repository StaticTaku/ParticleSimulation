#pragma once
#include <iostream>
#include <cmath>
#include <PhysicsValue/CoreData.hpp>
#include <CalculatedResult/CoreResult.hpp>
#include <Calculator/TreeComponent.hpp>

namespace BodysData
{
    double (*position)[3];
    double (*velocity)[3];
    double* mass;
    double (*accel)[3];
    double (*result)[3];
}


class CalcTreeSelfGrav
{
private:
    Cell* root = nullptr;
    double rsize = 2.0;
    Cell* freeCell = nullptr;
    Body* bodyArray;
    double theta;
    int usedCell; int maxLevel; int nfcalc; int n2bcalc; int n2bterm; int nbccalc; int nbcterm;
    double zero[3];
    bool firstCall = true;

    void MakeTree(int number)
    {
        Body* bodyP;
        NewTree();
        root = MakeCell();
        root->SetPosition(zero);
        ExpandBox(number);
        maxLevel = 0;
        for(int i = 0;i<number;++i)
            LoadBody(i);
        SetCenterOfMass(root,rsize);
        ThreadTree(root,nullptr);
    }

    void InitialMakeTree(int number);

    void NewTree()
    {
        CommonNode* p;
        Cell* cellP;

        static bool firstCall = true;
        
        if(!firstCall)
        {
            p = root;
            while(p != nullptr)
            {
                if(p->GetType() == NodeType::CELL)
                {
                    cellP = static_cast<Cell*>(p);
                    p->SetNext(freeCell);
                    freeCell = cellP;
                    p = cellP->GetMore();
                }else
                {
                    p = p->GetNext();
                }
            }
        }else
        {
            firstCall = false;
        }
        
        root = nullptr;
        usedCell = 0;
    }

    void InitialNewTree();

    void ExpandBox(int number)
    {
        double xyzmax = 0;
        for(int i = 0;i<number;++i)
            for(int j = 0;j<3;++j)
                xyzmax = std::max(xyzmax,std::abs(BodysData::position[i][j]-root->GetPosition(j)));
        
        while(rsize < 2*xyzmax)
            rsize *= 2;
    }

    void LoadBody(int i)
    {
        Cell* q; Cell* c;
        int qind,lev,k;
        double qsize;

        q = root;
        qind = GetSubIndex(i,q);
        qsize = rsize;
        lev = 0;
        while(q->GetSubp(qind) != nullptr)
        {
            if((q->GetSubp(qind))->GetType() == NodeType::BODY)
            {
                c = MakeCell();
                for(k = 0;k<3;++k)
                    c->SetPosition(k,q->GetPosition(k)+(BodysData::position[i][k] < q->GetPosition(k) ? -qsize:qsize)/4);

                c->SetSubp(GetSubIndex(static_cast<Body*>(q->GetSubp(qind))->GetBodyID(),c),q->GetSubp(qind));
                q->SetSubp(qind,c);
            }

            q = static_cast<Cell*>(q->GetSubp(qind));
            qind = GetSubIndex(i,q);
            qsize /= 2;
            ++lev;
        }

        q->SetSubp(qind,&bodyArray[i]);
        maxLevel = std::max(maxLevel,lev);
    }

    int GetSubIndex(int bodyID,Cell* cellP)
    {
        int ind;

        ind = 0;
        for(int k = 0;k<3;++k)
            if(cellP->GetPosition(k) <= BodysData::position[bodyID][k])
                ind += 8 >> (k+1);
        
        return ind;
    }

    void SetCenterOfMass(Cell* cellP,double pSize)
    {
        double cmpos[3] = {0};
        double tmpv[3] = {0};
        double mass;
        CommonNode* q;

        cellP->SetMass(0);
        
        for(int i=0;i<8;++i)
        {
            q = cellP->GetSubp(i);
            if(q != nullptr)
            {
                if(q->GetType() == NodeType::CELL)
                    SetCenterOfMass(static_cast<Cell*>(q),pSize/2);
                cellP->AddMass(q->GetMass());
                mass = q->GetMass();
                for(int k = 0;k<3;++k)
                {
                    tmpv[k] = q->GetPosition(k)*mass;
                    cmpos[k] += tmpv[k];
                }

            }
        }

        mass = cellP->GetMass();
        for(int i = 0;i<3;++i)
            cmpos[i] /= mass;
        
        SetCriticalRadius(cellP,cmpos,pSize);
        cellP->SetPosition(cmpos);
    }

    void SetCriticalRadius(Cell* cellP,double* cmpos,double psize)
    {
        double rc;
        if(theta == 0.0)
        {
            rc = 2*rsize;
        }else
        {
            rc = psize/theta;
        }

        cellP->SetCriticalRadius2(rc*rc);
    }

    void ThreadTree(CommonNode* p,CommonNode* n)
    {
        int ndesc,i;
        CommonNode* desc[9],*_p;
        Cell* _cellP;

        p->SetNext(n);
        if(p->GetType() == NodeType::CELL)
        {
            _cellP = static_cast<Cell*>(p);
            ndesc = 0;
            for(int i = 0;i<8;++i)
            {
                _p = _cellP->GetSubp(i);
                if(_p != nullptr)
                    desc[ndesc++] = _p;
            }
            _cellP->SetMore(desc[0]);
            desc[ndesc] = n;
            for(int i = 0;i<ndesc;++i)
                ThreadTree(desc[i],desc[i+1]);
        }
    }

    void SetAccelOf(int bodyID);

    void TreeScan(CommonNode* p);

    bool SubDivP(Cell* p, double* ri)
    {
        double dr[3];
        double length2 = 0;
        for(int i = 0;i<3;++i)
        {
            dr[i] = p->GetPosition(i) - ri[i];
            length2 += dr[i]*dr[i];
        }

        //pmem = p;
        return (length2 < p->GetCriticalRadius2());
    }

    Cell* MakeCell()
    {
        Cell* cellP;

        if(freeCell == nullptr)
        {

            cellP = new Cell;
            
        }
        else
        {
            cellP = freeCell;
            freeCell = static_cast<Cell*>(cellP->GetNext());
        }

        cellP -> ResetSubp();
        ++usedCell;
        return cellP;
    }
public:
    CalcTreeSelfGrav(double _theta):theta(_theta) {}

    void operator()(real _softParam, CoreData& data, CoreResult& result)
    {
        if(firstCall)
        {
            bodyArray = new Body[data.number];

            for(int i = 0;i<data.number;++i)
                bodyArray[i].SetBodyID(i);
            SetBodysPointer(data,result);
        }else
        {
            firstCall = false;
        }
        

        MakeTree(data.number);
        nfcalc = n2bcalc = nbccalc = 0;	

    #pragma omp parallel
        {
            int i,j,k;
            Body* bodyP;
            Cell* pmem;
            CommonNode* q;
            Cell* cellQ;
            double _accel[3], ri[3], dr[3];
            double ri2, mri5_2,lengthi_q;
            double eps2 = _softParam * _softParam;
            double phi0;

    #pragma omp for reduction(+: usedCell) reduction(+: maxLevel) reduction(+: nfcalc) reduction(+: n2bcalc) reduction(+: nbccalc) reduction(+: n2bterm) reduction(+: nbcterm)
            for(i = 0;i<data.number;++i)
            {
                bodyP = &bodyArray[i];
                
                data.potential[i] = 0;

                for(j = 0;j<3;++j)
                {   
                    _accel[j] = 0;
                    ri[j] = data.position[i][j];
                }

                n2bterm = 0; nbcterm = 0;
                
                q = root;
                while(q != nullptr)
                {
                    if(q->GetType() == NodeType::CELL &&
                        SubDivP((cellQ = static_cast<Cell*>(q)),ri))
                        q = (cellQ->GetMore());
                    else
                    {
                        if(q != bodyP)
                        {
                            for (k = 0; k < 3; ++k)
                                dr[k] = q->GetPosition(k) - ri[k];

                            lengthi_q = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] + eps2);
                            ri2 = 1.0 / lengthi_q;

                            mri5_2 = q->GetMass() * ri2 * sqrt(ri2);
                            lengthi_q = sqrt(lengthi_q);

                            for (k = 0; k < 3; ++k)
                                _accel[k] += mri5_2 * dr[k];

                            data.potential[i] += -data.mass[i]*q->GetMass()/lengthi_q;
                            
                            if(q->GetType() == NodeType::BODY)
                                ++n2bterm;
                            else
                                ++nbcterm;
                        }

                        q = q->GetNext();
                    }
                }

                for(j = 0;j<3;++j)
                {
                    result.accel[i][j] += _accel[j];
                    ++nfcalc;
                    n2bcalc += n2bterm;
                    nbccalc += nbcterm;
                }

            }   
    #pragma omp barrier
        }
    }   

    void SetTheta(double _theta)
    {
        theta = _theta;
    }
    
	std::string WhoAmI() const 
	{
		return "IForceCalculator::TreeCalculator";
	}

    ~CalcTreeSelfGrav() 
    {
        CommonNode* p,*q;

        p = root;
        while(true)
        {
            if(p != nullptr)
            {
                if(p->GetType() == NodeType::CELL)
                {
                    q = static_cast<Cell*>(p)->GetMore();
                    delete p;
                    p = q;
                }else
                {
                    p = p->GetNext();
                }
            }else
            {
                break;
            }
        }

        if(freeCell != nullptr)
            delete freeCell;

        if(bodyArray != nullptr)
            delete[] bodyArray;
    }
};