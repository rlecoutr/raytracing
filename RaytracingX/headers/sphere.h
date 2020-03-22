#ifndef SPHERE_H
#define SPHERE_H

#include "vecteur.h"
#include "ray.h"
#include <math.h>
#include "object.h"

class Sphere : public Object
{
    Vecteur _O;
    double _r;
    Vecteur _color;
    bool _isSpeculaire = false;
    bool _isTransparent = false;
    double _n = 1.5;

public:
    Sphere():_O(),_r(),_color() {}
    Sphere(Vecteur O_, double r_, Vecteur color_):_O(O_),_r(r_),_color(color_) {}
    Sphere(Vecteur O_, double r_, Vecteur color_, bool isSpeculaire_):_O(O_),_r(r_),_color(color_) ,_isSpeculaire(isSpeculaire_) {}
    Sphere(Vecteur O_, double r_, Vecteur color_, bool isSpeculaire_, bool isTransparent_):
        _O(O_),_r(r_),_color(color_) ,_isSpeculaire(isSpeculaire_),_isTransparent(isTransparent_) {}

    Vecteur O() const { return _O; }
    double r() const { return _r; }
    Vecteur color() { return _color; }
    bool isSpeculaire() { return _isSpeculaire; }
    bool isTransparent() { return _isTransparent; }
    double n() { return _n; }

    void rotate(Vecteur axe, double angle) {}
    CurrentObject intersection (Ray R);
    Vecteur getNormal(Vecteur P) { Vecteur n = P-_O; n.normalisation(); return n; }
};

#endif // SPHERE_H
