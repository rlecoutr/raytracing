#include "ray.h"

Ray Ray::reflechir(Vecteur n, Vecteur P)
{
    double eps = 0.000001;

    Vecteur i = _u;
    i.normalisation();

    return Ray(P+n*eps, i - n*2*i.dot(n));
}

Ray Ray::refracter(Vecteur n, Vecteur P, double nSphere)
{
    double nAir = 1.00;
    double eps = 0.000001;

    Vecteur i = _u;
    i.normalisation();

    if (i.dot(n) < 0) // On rentre dans la sphère
    {
        double racine = 1-(nAir/nSphere)*(nAir/nSphere)*(1-i.dot(n)*i.dot(n));
        if (racine < 0)
            racine = 0.0;

        return Ray(P-n*eps, i*(nAir/nSphere) - n*(nAir/nSphere*i.dot(n) + sqrt(racine)));
    }
    else // On quitte la sphère
    {
        n = -n;
        double racine = 1-(nSphere/nAir)*(nSphere/nAir)*(1-i.dot(n)*i.dot(n));
        if (racine < 0)
            racine = 0.0;

        return Ray(P-n*eps, i*(nSphere/nAir) - n*(nSphere/nAir*i.dot(n) + sqrt(racine)));
    }
}

