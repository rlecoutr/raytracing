#ifndef MESH_H
#define MESH_H

#include "vecteur.h"
#include "ray.h"
#include "object.h"
#include "bbox.hpp"

#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <list>

#include "stb_image.h"

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
    }
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup;
    
    int& operator[](int i) { if (i==1) return vtxi; if (i==2) return vtxj; return vtxk; }
};



class Mesh : public Object
{
    Vecteur _color = Vecteur(1, 1, 1);
    bool _isSpeculaire = false;
    bool _isTransparent = false;
    double _n = 1.5;
    
    std::string textureFolder;
    bool hasTexture = true;

    int currentTriangle = -1;

    std::vector<TriangleIndices> indices;
    std::vector<Vecteur> vertices;
    std::vector<Vecteur> normals;
    std::vector<Vecteur> uvs;
    std::vector<Vecteur> vertexcolors;

    std::vector<std::vector<unsigned char>> textures;
    std::vector<int> w, h;
    
    BVH bvh;

public:
    Mesh() {}
    Mesh(const char* obj, std::string folder, double scaling, const Vecteur& offset);
    Mesh(const char* obj, std::string folder, double scaling, const Vecteur& offset, bool isSpeculaire_);
    Mesh(const char* obj, std::string folder, double scaling, const Vecteur& offset, bool isSpeculaire_, bool isTransparent_);

    void readOBJ(const char* obj);
    void add_texture(const char* filename);
    void rotate(Vecteur axe, double angle);
    
    Bbox setBox(int i0, int i1);
    void setBVH(BVH* node, int i0, int i1);

    Vecteur color() { return _color; }
    bool isSpeculaire() { return _isSpeculaire; }
    bool isTransparent() { return _isTransparent; }
    double n() { return _n; }

    CurrentObject intersection(Ray R);
    Vecteur getNormal(Vecteur P);

    Vecteur normalTriangle(TriangleIndices t, Vecteur P);
    Vecteur uvTriangle(TriangleIndices t, Vecteur P);
    
    double intersectPlan(Ray R, Vecteur P0, Vecteur n);
    bool pointInTriangle(TriangleIndices t, Vecteur P);
};

#endif // MESH_H
