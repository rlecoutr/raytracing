//
//  bbox.cpp
//  RaytracingX
//
//  Created by Raphael Lecoutre on 30/01/2020.
//  Copyright © 2020 Raphael Lecoutre. All rights reserved.
//

#include "bbox.hpp"

bool Bbox::intersection(Ray r)
{
    // X
    double tx1 = (_bmin[0] - r.C()[0])/r.u()[0];
    double tx2 = (_bmax[0] - r.C()[0])/r.u()[0];
    double txmin = std::min(tx1, tx2);
    double txmax = std::max(tx1, tx2);
    //Y
    double ty1 = (_bmin[1] - r.C()[1])/r.u()[1];
    double ty2 = (_bmax[1] - r.C()[1])/r.u()[1];
    double tymin = std::min(ty1, ty2);
    double tymax = std::max(ty1, ty2);
    // Z
    double tz1 = (_bmin[2] - r.C()[2])/r.u()[2];
    double tz2 = (_bmax[2] - r.C()[2])/r.u()[2];
    double tzmin = std::min(tz1, tz2);
    double tzmax = std::max(tz1, tz2);
    
    double tmin = std::max(std::max(txmin, tymin), tzmin);
    double tmax = std::min(std::min(txmax, tymax), tzmax);
    
    if (tmax > tmin && tmax > 0)
        return true;

    return false;
}
