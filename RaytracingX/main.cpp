#include <iostream>
#include <algorithm>
#include "vecteur.h"
#include "sphere.h"
#include "mesh.h"
#include "ray.h"
#include "scene.h"
#include <math.h>
#include "omp.h"

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


Ray generateRay(int i, int j, double depth, int W, int H, Vecteur Position, double distance_focale, double ouverture_focale) {

    double x = drand48(), y = drand48(), R=sqrt(-2*log(x));
    double u =  R*cos(2*3.1416*y)*0.5;
    double v =  R*sin(2*3.1416*y)*0.5;

    Vecteur V(j-W/2-0.5, H/2-i+0.5, -depth);
    V.normalisation();
    Vecteur directionTemp(j-W/2+u-0.5, H/2-i-v+0.5, -depth);
    directionTemp.normalisation();

    double dx = (drand48()-0.5)*ouverture_focale;
    double dy = (drand48()-0.5)*ouverture_focale;

    Vecteur desti = Position + directionTemp*distance_focale;
    Vecteur origine = Position + Vecteur(dx, dy, 0);
    Vecteur directionFinale = desti - origine;
    directionFinale.normalisation();

    return Ray(origine, directionFinale); // Avec anti-aliasing
    //return Ray(Position, V); // Sans anti-aliasing
}


double setRight(double n)
{
    return floor(std::min(255.0, pow(std::max(0.0, n), 0.45)));
}

int main() {
    int W = 512;
    int H = 512;

    double fov = M_PI/3;
    double depth = H/(2*tan(fov*0.5));

    int nb_rebonds = 5;
    double nb_alea = 50.0;

    double ouverture_focale = 0.0;
    double distance_focale = 55.0;

    Object *sphere1 = new Sphere(Vecteur(0, 0, 0), 4, Vecteur(1,1,1), true, false);
    Object *sphere2 = new Sphere(Vecteur(10, 0, 20), 10, Vecteur(1,1,1), false, false);
    Object *mur1 = new Sphere(Vecteur(0,0,-1000),940, Vecteur(1,0,0));
    Object *mur2 = new Sphere(Vecteur(0,-1000,0),970, Vecteur(0,0,1));
    Object *mur3 = new Sphere(Vecteur(0,0,1000),940, Vecteur(0,1,0));
    Object *mur4 = new Sphere(Vecteur(0,1000,0),940, Vecteur(0,1,1));
    Object *mur5 = new Sphere(Vecteur(1000,0,0),940, Vecteur(1,1,0));
    Object *mur6 = new Sphere(Vecteur(-1000,0,0),940, Vecteur(1,0,1));
    
    Sphere L(Vecteur(-10, 20, 40), 5, Vecteur(1,1,1));

    //Object *mercedes = new Mesh("mercedes/sls_amg.obj", "mercedes", 10.0, Vecteur(0, -20, 0));
    //mercedes->rotate(Vecteur(0, 1, 0), 90.0);
    
    Object *dracaufeu = new Mesh("dracaufeu/Charizard.obj", "dracaufeu", 2., Vecteur(0, -20, 0));
    dracaufeu->rotate(Vecteur(0, 1, 0), 10.0);
    
    Object *girl = new Mesh("girl/girl.obj", "girl", 20.0, Vecteur(-20, -20, 10));
    
    //Object *cube = new Mesh("cube/cube.obj", "cube", 1, Vecteur(0, 0, 0));
    //cube->rotate(Vecteur(0, 1, 0), 10.0);

    Vecteur centreCamera(0, 0, 55);

    Scene scene;
    //scene.addObject(sphere1);
    //scene.addObject(sphere2);
    //scene.addObject(girl);
    //scene.addObject(mercedes);
    scene.addObject(dracaufeu);
    //scene.addObject(cube);
    scene.addObject(mur1);
    scene.addObject(mur2);
    scene.addObject(mur3);
    scene.addObject(mur4);
    scene.addObject(mur5);
    scene.addObject(mur6);
    scene.addLumiere(L);

    std::vector<unsigned char> image(W*H * 3, 0);

    #pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vecteur couleur(0.0, 0.0, 0.0);

            for (int k=0; k<nb_alea; k++)
            {
                Ray R = generateRay(i, j, depth, W, H, centreCamera, distance_focale, ouverture_focale);
                couleur = couleur + scene.getColor(R, nb_rebonds);
            }

            couleur = couleur*(1/nb_alea);

            image[(i*W + j) * 3 + 0] = setRight(couleur.x());
            image[(i*W + j) * 3 + 1] = setRight(couleur.y());
            image[(i*W + j) * 3 + 2] = setRight(couleur.z());
        }
        if (i % 5 == 0)
            std::cout << i << std::endl;
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
