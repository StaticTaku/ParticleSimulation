#pragma once
#include <Settings.hpp>
#include <cmath>

/*
x<0:dens_l,vel_l,pres_l
x>0:dens_r,vel_r,pres_r
を初期条件とする衝撃波管問題の、x=0での解
(x=0での解は時間によらす定数)
ans[0] = density(x=0);
ans[1] = velocity(x=0);
ans[2] = pressure(x=0);

*/

#if defined(GodunovMethod)
class RiemanSolver
{
public:
    inline void GetShockTubeAnswerAtZero(real dens_l, real vel_l, real pres_l, real dens_r, real vel_r, real pres_r, real* ans)
    {
        real c_l = std::sqrt(heatCapRatio*pres_l/dens_l);
        if(vel_l-c_l>0)
        {
            ans[0] = dens_l;
            ans[1] = vel_l;
            ans[2] = pres_l;
            return;
        }

        real c_r = std::sqrt(heatCapRatio*pres_r/dens_r);
        if(vel_r + c_r < 0)
        {
            ans[0] = dens_r;
            ans[1] = vel_r;
            ans[2] = pres_r;
            return;
        }


        real vel,pres,dens1,dens2,Ml,Mr;

        pres = 1;//適当に初期値を決める。
        real past_pres;
        do
        {
            past_pres = pres;
            Ml = Ml_(dens_l,pres_l,pres);
            Mr = Mr_(dens_r,pres_r,pres);

            pres = omega*pres + (1-omega)*P(Ml,Mr,pres_l,pres_r,vel_l,vel_r);
        } while (std::abs(pres-past_pres) >= epsilon);
        
        vel = u(pres_l,pres_r,vel_l,vel_r,Ml,Mr);
        dens1 = dens_1(pres_l,dens_l,vel_l,Ml,pres,vel);

        real c1 = std::sqrt(heatCapRatio*pres/dens1);

        if(vel_l-c_l <= 0 && vel - c1 > 0)
        { 
            ans[1] = vel_phi1(0,vel_l,c_l);
            real cphi1 = c_phi1(0,vel_l,c_l);
            ans[0] = dens_phi(dens_l,pres_l,cphi1);
            ans[2] = pres_phi(ans[0],cphi1);
            return;
        }

        if(vel-c1 <= 0 && vel > 0)
        {
            ans[0] = dens1;
            ans[1] = vel;
            ans[2] = pres;
            return;
        }

        dens2 = dens_2(pres_r,dens_r,vel_r,Ml,pres,vel);
        real c2 = std::sqrt(heatCapRatio*pres/dens2);

        if(vel <= 0 && vel + c2 > 0)
        {
            ans[0] = dens2;
            ans[1] = vel;
            ans[2] = pres;
            return;
        }

        if(vel + c2 <= 0 && vel_r + c_r >= 0)
        {
            ans[1] = vel_phi2(0,vel_r,c_r);
            real cphi2 = c_phi2(0,vel_r,c_r);
            ans[0] = dens_phi(dens_r,pres_r,cphi2);
            ans[2] = pres_phi(ans[0],cphi2);
            return;
        }

        return;

    }

private:
    real Ml_(real dens_l, real pres_l, real pres)
    {
        if(pres_l > pres)
        {
            return std::sqrt(dens_l*pres_l)*(heatCapRatio-1)/(2*std::sqrt(heatCapRatio))*(1-pres/pres_l)/(1-std::pow(pres/pres_l,(heatCapRatio-1)/(2*heatCapRatio)));
        }else
        {
            return std::sqrt(dens_l*pres_l)*std::pow((heatCapRatio-1)/2 + (heatCapRatio+1)/2 * pres/pres_l,0.5);
        }
    }

    real Mr_(real dens_r, real pres_r, real pres)
    {
        if(pres_r > pres)
        {
            return std::sqrt(dens_r*pres_r)*(heatCapRatio-1)/(2*std::sqrt(heatCapRatio))*(1-pres/pres_r)/(1-std::pow(pres/pres_r,(heatCapRatio-1)/(2*heatCapRatio)));
        }else
        {
            return std::sqrt(dens_r*pres_r)*std::pow((heatCapRatio-1)/2 + (heatCapRatio+1)/2 * pres/pres_r,0.5);
        }
    }

    real P(real Ml, real Mr, real pres_l, real pres_r, real vel_l, real vel_r)
    {
        return ((vel_l-vel_r) + (pres_l/Ml + pres_r/Mr))/(1/Ml + 1/Mr);
    }

    real u(real pres_l, real pres_r, real vel_l, real vel_r, real Ml, real Mr)
    {
        return (Ml*vel_l+Mr*vel_r+pres_l-pres_r)/(Ml+Mr);
    }

    real dens_1(real pres_l, real dens_l, real vel_l, real Ml, real pres, real vel)
    {
        if(pres_l>pres)
        {
            return dens_l*std::pow(pres/pres_l,1/heatCapRatio);
        }else
        {
            return dens_l*Ml/(dens_l*(vel - vel_l) + Ml);
        }
    }

    real dens_2(real pres_r, real dens_r, real vel_r, real Mr, real pres, real vel)
    {
        if(pres_r>pres)
        {
            return dens_r*std::pow(pres/pres_r,1/heatCapRatio);
        }else
        {
            return dens_r*Mr/(dens_r*(vel - vel_r) + Mr);
        }
    }

    real vel_phi1(real phi, real vel_l, real c_l)
    {
        return (heatCapRatio-1)/(heatCapRatio+1) * (2/(heatCapRatio-1) * phi + vel_l + 2/(heatCapRatio-1)*c_l);
    }

    real vel_phi2(real phi, real vel_r, real c_r)
    {
        return (heatCapRatio-1)/(heatCapRatio+1) * (2/(heatCapRatio-1) * phi + vel_r - 2/(heatCapRatio-1)*c_r);
    }

    real c_phi1(real phi, real vel_l, real c_l)
    {
        return (heatCapRatio-1)/(heatCapRatio+1) * (-phi + vel_l + 2/(heatCapRatio-1)*c_l);
    }

    real c_phi2(real phi, real vel_r, real c_r)
    {
        return (heatCapRatio-1)/(heatCapRatio+1) * (phi - vel_r + 2/(heatCapRatio-1)*c_r);
    }

    real dens_phi(real dens_, real pres_, real c_phi)
    {
        return std::pow((std::pow(dens_,heatCapRatio)*c_phi*c_phi)/(heatCapRatio*pres_),1/(heatCapRatio-1));
    }

    real pres_phi(real dens_phi, real c_phi)
    {
        return dens_phi*c_phi*c_phi/heatCapRatio;
    }
};

#endif