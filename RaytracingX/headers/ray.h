#ifndef RAY_H
#define RAY_H

#include "vecteur.h"

class Ray
{
    Vecteur _C, _u;


public:
    Ray():_C(), _u() {}
    Ray(Vecteur C_, Vecteur u_):_C(C_), _u(u_) {}

    Vecteur C() const { return _C; }
    Vecteur u() const { return _u; }

    Ray reflechir(Vecteur n, Vecteur P);
    Ray refracter(Vecteur n, Vecteur P, double nSphere);
};

#endif // RAY_H
