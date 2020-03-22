#ifndef OBJECT_H
#define OBJECT_H

#include "ray.h"
#include "vecteur.h"

typedef struct
{
    int index;
    double t;
    Vecteur normal;
    Vecteur couleur;
    bool isSpeculaire;
    bool isTransparent;
    double n;
    
} CurrentObject;

class Object
{
public:
    Object();
    virtual CurrentObject intersection(Ray R) = 0;
    virtual Vecteur getNormal(Vecteur P) = 0;
    virtual bool isSpeculaire() = 0;
    virtual bool isTransparent() = 0;
    virtual Vecteur color() = 0;
    virtual double n() = 0;

    virtual void rotate(Vecteur axe, double angle) = 0;
};

#endif // OBJECT_H
