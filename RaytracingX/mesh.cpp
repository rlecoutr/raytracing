#include "mesh.h"

Mesh::Mesh(const char *obj, std::string folder, double scaling, const Vecteur &offset)
{
    textureFolder = folder;
    readOBJ(obj);
    
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] * scaling + offset;
    }
    
    setBVH(&bvh, 0, indices.size());
}

Mesh::Mesh(const char *obj, std::string folder, double scaling, const Vecteur &offset, bool isSpeculaire_) : _isSpeculaire(isSpeculaire_)
{
    textureFolder = folder;
    readOBJ(obj);
    
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] * scaling + offset;
    }
    
    setBVH(&bvh, 0, indices.size());
}

Mesh::Mesh(const char *obj, std::string folder, double scaling, const Vecteur &offset, bool isSpeculaire_, bool isTransparent_) : _isSpeculaire(isSpeculaire_), _isTransparent(isTransparent_)
{
    textureFolder = folder;
    readOBJ(obj);
    
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] * scaling + offset;
    }
    
    setBVH(&bvh, 0, indices.size());
}

void Mesh::rotate(Vecteur axe, double angle)
{
    axe.normalisation();

    for (int i=0; i<vertices.size(); i++)
        vertices[i] = vertices[i].rotation(axe, angle);
    for (int i=0; i<normals.size(); i++)
        normals[i] = normals[i].rotation(axe, angle);
    
    setBVH(&bvh, 0, indices.size());
}

void Mesh::readOBJ(const char *obj)
{
    char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj, "r");

    std::map<std::string, int> groupNames;
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\r]\n", grp);
            if (groupNames.find(std::string(grp)) != groupNames.end()) {
                curGroup = groupNames[std::string(grp)];
            }
            else {
                curGroup = groupNames.size();
                groupNames[std::string(grp)] = curGroup;
            }
        }
        if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
            sscanf(line, "mtllib %[^\r]\n", matfile);
        }
        if (line[0] == 'v' && line[1] == ' ') {
            Vecteur vec;
            Vecteur col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
                vertices.push_back(vec);
                vertexcolors.push_back(col);
            }
            else {
                if (textureFolder != "girl")
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);  // cube / dracaufeu
                else
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // girl
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vecteur vec;
            if (textureFolder != "girl")
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]); //cube / dracaufeu
            else
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vecteur vec(0., 0., 0.);
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;

            char* consumedline = line + 1;
            int offset;
            t.faceGroup = curGroup;
            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;

                indices.push_back(t);
            }
            else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }


            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.faceGroup = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }

        }
    }
    fclose(f);
    
    const char* name = (textureFolder+"/"+std::string(matfile)).c_str();
    FILE* fi;
    fi = fopen(name, "r");
        
    if (fi)
    {
        while (!feof(fi))
        {
            char line[255];
            fgets(line, 255, fi);
            if (line[0] == 'm' && line[4] == 'K' && line[5] == 'd')
            {
                char textureFile[255];
                sscanf(line, "map_Kd %100s\n", textureFile);
                std::cout << (textureFolder+"/"+std::string(textureFile)).c_str() << std::endl;
                add_texture((textureFolder+"/"+std::string(textureFile)).c_str());
            }
        }
        fclose(fi);
    }
    else
    {
        std::cout << "Couldn't open " << name << std::endl;
        hasTexture = false;
    }
}

void Mesh::add_texture(const char* filename)
{
    textures.resize(textures.size() + 1);
    w.resize(w.size() + 1);
    h.resize(h.size() + 1);
    
    int x, y, n;
    unsigned char *data = stbi_load(filename, &x, &y, &n, 3);
    
    if (data)
    {
        w[w.size() - 1] = x;
        h[h.size() - 1] = y;
        
        int size = 3 * w[w.size() - 1] * h[h.size() - 1];
        textures[textures.size() - 1].resize(0);
        textures[textures.size() - 1].insert(textures[textures.size() - 1].end(), data, data+size);
        std::cout << "Material loaded" << std::endl;
    }
    else
        std::cout << "Material failed loading" << std::endl;
}

double Mesh::intersectPlan(Ray R, Vecteur P0, Vecteur n)
{
    double dot1 = R.u().dot(n);
    double t = -1;

    if (dot1 != 0)
        t = -(R.C()-P0).dot(n)/dot1;

    return t;
}

bool Mesh::pointInTriangle(TriangleIndices t, Vecteur P)
{
    Vecteur v0 = vertices[t.vtxj] - vertices[t.vtxi];
    Vecteur v1 = vertices[t.vtxk] - vertices[t.vtxi];
    Vecteur v2 = P - vertices[t.vtxi];

    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v < 1);
}

Vecteur Mesh::normalTriangle(TriangleIndices t, Vecteur P)
{
    Vecteur v0 = vertices[t.vtxj] - vertices[t.vtxi];
    Vecteur v1 = vertices[t.vtxk] - vertices[t.vtxi];
    Vecteur v2 = P - vertices[t.vtxi];
    
    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);
    
    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    return normals[t.ni]*(1-u-v) + normals[t.nj]*u + normals[t.nk]*v;
}

Vecteur Mesh::uvTriangle(TriangleIndices t, Vecteur P)
{
    Vecteur v0 = vertices[t.vtxj] - vertices[t.vtxi];
    Vecteur v1 = vertices[t.vtxk] - vertices[t.vtxi];
    Vecteur v2 = P - vertices[t.vtxi];
    
    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);
    
    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    Vecteur result = uvs[t.uvi]*(1-u-v) + uvs[t.uvj]*u + uvs[t.uvk]*v;
    
    return result;
}

CurrentObject Mesh::intersection(Ray R)
{
    CurrentObject co;
    co.t = -1;
    
    double tmin = 1e99;
    double max = tmin;
    
    if (!bvh.bb.intersection(R)) return co;
    
    std::list<BVH*> l;
    l.push_front(&bvh);
    
    while (!l.empty())
    {
        BVH* current = l.front();
        l.pop_front();
        
        if (current->fg && current->fg->bb.intersection(R))
            l.push_back(current->fg);
        if (current->fd && current->fd->bb.intersection(R))
            l.push_back(current->fd);
        
        if (!current->fg)
        {
            for (int i=current->i0; i<current->i1; i++)
            {
                Vecteur nt = (vertices[indices[i].vtxj] - vertices[indices[i].vtxi])^(vertices[indices[i].vtxk] - vertices[indices[i].vtxi]);
                
                double t = intersectPlan(R, vertices[indices[i].vtxi], nt);
                
                if (t > 0) // Si le plan du triangle est dans l'axe du rayon
                {
                    Vecteur P = R.C() + R.u()*t; // Point d'intersection dans le plan du triangle
                    if (pointInTriangle(indices[i], P) && t<tmin)
                    {
                        tmin = t;
                        co.t = t;
                        co.normal = normalTriangle(indices[i], P);
                        co.normal.normalisation();
                        co.isTransparent = isTransparent();
                        co.isSpeculaire = isSpeculaire();
                        co.n = n();
                        
                        if (hasTexture)
                        {
                            Vecteur uv = uvTriangle(indices[i], P);
                            
                            int x = uv[0]*(w[indices[i].faceGroup]);
                            int y = (1-uv[1])*(h[indices[i].faceGroup]);
                            
                            double cr = textures[indices[i].faceGroup][(y*w[indices[i].faceGroup] + x) * 3] / 255.0;
                            double cg = textures[indices[i].faceGroup][(y*w[indices[i].faceGroup] + x) * 3 + 1] / 255.0;
                            double cb = textures[indices[i].faceGroup][(y*w[indices[i].faceGroup] + x) * 3 + 2] / 255.0;
                            
                            co.couleur = Vecteur(cr, cg, cb);
                        }
                        else
                            co.couleur = color();
                    }
                }
            }
        }
    }

    if (tmin == max)
        co.t = -1;

    return co;
}

Vecteur Mesh::getNormal(Vecteur P)
{
    Vecteur n = normalTriangle(indices[currentTriangle], P);
    n.normalisation();
    return n;
}

Bbox Mesh::setBox(int i0, int i1)
{
    Bbox box;
    
    box._bmin = vertices[indices[i0].vtxi];
    box._bmax = vertices[indices[i0].vtxi];
    
    for (int i = i0; i < i1; i++) { // Triangles
        for (int j=0; j<3; j++) { // Sommets du triangle
            for (int k=0; k < 3; k++) // Coordonnées du sommet
            {
                box._bmax[k] = std::max(box._bmax[k], vertices[indices[i][j]][k]);
                box._bmin[k] = std::min(box._bmin[k], vertices[indices[i][j]][k]);
            }
        }
    }
    
    return box;
}

void Mesh::setBVH(BVH* node, int i0, int i1)
{
    if (i1 <= i0) return;
        
    node->bb = setBox(i0, i1);
    node->i0 = i0;
    node->i1 = i1;
    node->fg = NULL;
    node->fd = NULL;
    
    Vecteur diag = node->bb._bmax - node->bb._bmin;
    int dim;
    
    if (diag[0] > diag[1] && diag[0] > diag[2])
        dim = 0;
    else if (diag[1] > diag[0] && diag[1] > diag[2])
        dim = 1;
    else
        dim = 2;
    
    int pivot = i0-1;
    double split = node->bb._bmin[dim] + diag[dim]*0.5;
    
    for (int i=i0; i<i1; i++)
    {
        double center = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim])/3.0;
        
        if (center < split)
        {
            pivot++;
            std::swap(indices[i], indices[pivot]);
        }
    }
    
    if (pivot < i0 || pivot >= i1-1 || i1 == i0+1) return;
    
    node->fg = new BVH();
    setBVH(node->fg, i0, pivot+1);
    
    node->fd = new BVH();
    setBVH(node->fd, pivot+1, i1);
}
