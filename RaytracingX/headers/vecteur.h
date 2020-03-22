#ifndef VECTEUR_H
#define VECTEUR_H

#include <math.h>

class Vecteur
{
    double _x;
    double _y;
    double _z;

public:
    Vecteur():_x(),_y(),_z() {}
    Vecteur(double x_, double y_, double z_):_x(x_),_y(y_),_z(z_) {}

    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }

    Vecteur operator-() const { return Vecteur(-_x, -_y, -_z); }
    double dot(Vecteur y);

    double normeCarre() const { return (_x*_x + _y*_y + _z*_z); }
    double norme() const { return sqrt(normeCarre()); }
    void normalisation();
    double& operator[](int i) { if (i==0) return _x; else if (i==1) return _y; else return _z; }

    Vecteur rotation(Vecteur axe, double angle);
};

Vecteur operator*(Vecteur const& x, double a);
Vecteur operator*(Vecteur const& x, Vecteur const& y);
Vecteur operator+(Vecteur const& x, Vecteur const& y);
Vecteur operator-(Vecteur const& x, Vecteur const& y);
Vecteur operator^(Vecteur const& x, Vecteur const& y);

#endif // VECTEUR_H
