#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "sphere.h"
#include "ray.h"
#include <algorithm>
#include <iostream>
#include <random>

class Scene
{
    std::vector<Object*> objectTab;
    Sphere lumiere;

public:
    std::default_random_engine engine;

    Scene() {}

    void addObject(Object *O) { objectTab.push_back(O); }
    void addLumiere(Sphere L) { lumiere = L; }

    CurrentObject intersection(Ray R);
    std::vector<Object*> objects() const { return objectTab; }
    Vecteur getColor(Ray R, int nRebond);
    
    Vecteur randomCos(Vecteur n);
};

#endif // SCENE_H
