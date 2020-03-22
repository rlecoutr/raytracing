#include "scene.h"

CurrentObject Scene::intersection(Ray R)
{
    CurrentObject co;
    CurrentObject coTemp;

    co.index = -1;

    double tmin = 1e99;
    double max = tmin;

    for (int i=0; i<objectTab.size(); i++)
    {
        coTemp = objectTab[i]->intersection(R);

        if (coTemp.t != -1 && coTemp.t < tmin)
        {
            tmin = coTemp.t;
            co = coTemp;
            co.index = i;
        }
    }

    if (tmin == max)
    {
        co.index = -1;
    }

    return co;
}

Vecteur Scene::randomCos(Vecteur n)
{
    std::uniform_real_distribution<double> distrib(0,1);
    double r1 = distrib(engine);
    double r2 = distrib(engine);
    double x = cos(2*M_PI*r1)*sqrt(1-r2);
    double y = sin(2*M_PI*r1)*sqrt(1-r2);
    double z = sqrt(r2);
    
    Vecteur rando(distrib(engine) - 0.5, distrib(engine) - 0.5, distrib(engine) - 0.5);
    Vecteur u2 = n^rando;
    u2.normalisation();
    Vecteur u1 = u2^n;
    
    return u1*x+u2*y+n*z;
}

Vecteur Scene::getColor(Ray R, int nRebond)
{
    Vecteur couleur(0.0, 0.0, 0.0);
    double eps = 0.000001;

    CurrentObject co = intersection(R);

    if (co.index != -1)
    {
        Vecteur P = R.C() + R.u()*co.t; // Point d'intersection P
        Vecteur n = co.normal;

        if (co.isSpeculaire && nRebond > 0) // Si on a affaire à une surface miroir
        {
            couleur = getColor(R.reflechir(n, P), nRebond-1);
        }
        else if (co.isTransparent && nRebond > 0) // Si on a affaire à une surface transparente
        {
            couleur = getColor(R.refracter(n, P, co.n), nRebond-1);
        }
        else // Surface diffuse
        {
            double coef_diffu = 0.25;
            Vecteur couleurDiffu(0.0, 0.0, 0.0);
            
            
            // Gestion de la diffusion
            Vecteur d = randomCos(n);
            Ray RDiffu(P+n*eps, d);
            if (nRebond > 0)
            {
                Vecteur couleurTemp = getColor(RDiffu, nRebond-1)*coef_diffu;
                couleurDiffu = couleurTemp*co.couleur;
            }
            
            // Gestion de la lumière (et ombre portée)
            
            Vecteur l = lumiere.O() - P; // Vecteur de P vers le centre de la lumière
            l.normalisation();
            Vecteur OX = randomCos(-l);
            Vecteur x = lumiere.O()+OX*(lumiere.r()+eps);
            Vecteur w = x-P;
            double d2 = w.normeCarre();
            w.normalisation();
            
            Ray RLumiere(P+n*eps, w);
            double t2 = intersection(RLumiere).t;
            
            if (t2 > 0 && t2*t2 < d2) // S'il existe une sphère qui bloque la lumière, on voit du noir
                couleur = couleurDiffu;
            else
            {
                Vecteur np = x-lumiere.O();
                np.normalisation();
                
                double costheta = std::max(0.0, n.dot(w));
                double costhetap = std::max(0.0, np.dot(-w));
                double costhetapp = std::max(0.0, (-l).dot(np));
                
                double I = 1000000000;
                
                Vecteur value = lumiere.color()*(I*costheta*costhetap/(M_PI*costhetapp*d2));
                
                couleur = couleurDiffu + co.couleur*value;
            }
        }
    }

    return couleur;
}
