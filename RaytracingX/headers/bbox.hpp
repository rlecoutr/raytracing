//
//  bbox.hpp
//  RaytracingX
//
//  Created by Raphael Lecoutre on 30/01/2020.
//  Copyright © 2020 Raphael Lecoutre. All rights reserved.
//

#ifndef bbox_hpp
#define bbox_hpp

#include "vecteur.h"
#include "ray.h"
#include <algorithm>

class Bbox
{
public:
    Bbox() {}
    Bbox(Vecteur bmin_, Vecteur bmax_) : _bmin(bmin_), _bmax(bmax_) {}
    Vecteur _bmin, _bmax;
    bool intersection(Ray r);
};

class BVH
{
public:
    int i0, i1;
    Bbox bb;
    BVH *fg, *fd;
};

#endif /* bbox_hpp */
