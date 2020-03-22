#include "sphere.h"

CurrentObject Sphere::intersection(Ray R)
{
    CurrentObject co;
    co.t = -1;
    co.couleur = color();
    co.isSpeculaire = isSpeculaire();
    co.isTransparent = isTransparent();
    co.n = n();
    
    double a = 1;
    double b = 2*(R.u().dot((R.C()-_O)));
    double c = (R.C()-_O).normeCarre() - _r*_r;

    double d = b*b - 4*a*c;

    if (d == 0)
    {
         double t = -b/(2*a);
         if (t > 0)
             co.t = t;
    }
    else if (d > 0)
    {
        double t2 = (-b+sqrt(d))/(2*a);
        double t1 = (-b-sqrt(d))/(2*a);

        if (t2>0)
        {
            if (t1<0)
                co.t = t2;
            else
                co.t = t1;
        }
    }

    co.normal = getNormal(R.C()+R.u()*co.t);
    return co;
}
