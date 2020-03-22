#include "vecteur.h"

void Vecteur::normalisation()
{
    double n = norme();

    _x /= n;
    _y /= n;
    _z /= n;
}

Vecteur operator*(Vecteur const& x, double a)
{
    Vecteur result(x.x()*a, x.y()*a, x.z()*a);

    return result;
}

Vecteur operator*(Vecteur const& x, Vecteur const& y)
{
    Vecteur result(x.x()*y.x(), x.y()*y.y(), x.z()*y.z());
    
    return result;
}

Vecteur operator+(Vecteur const& x, Vecteur const& y)
{
    Vecteur result(x.x()+y.x(), x.y()+y.y(), x.z()+y.z());

    return result;
}

Vecteur operator-(Vecteur const& x, Vecteur const& y)
{
    Vecteur result(x.x()-y.x(), x.y()-y.y(), x.z()-y.z());

    return result;
}

double Vecteur::dot(Vecteur y)
{
    double result = _x*y.x() + _y*y.y() + _z*y.z();

    return result;
}

Vecteur operator^(Vecteur const& x, Vecteur const& y)
{
    double x_, y_, z_;

    x_ = x.y()*y.z() - x.z()*y.y();
    y_ = x.z()*y.x() - x.x()*y.z();
    z_ = x.x()*y.y() - x.y()*y.x();

    Vecteur result(x_,y_,z_);

    return result;
}

Vecteur Vecteur::rotation(Vecteur axe, double angle)
{
    double c = cos(angle*M_PI/180);
    double s = sin(angle*M_PI/180);

    double R11 = axe.x()*axe.x()*(1-c) + c;
    double R12 = axe.x()*axe.y()*(1-c) - axe.z()*s;
    double R13 = axe.x()*axe.z()*(1-c) + axe.y()*s;
    double R21 = axe.x()*axe.y()*(1-c) + axe.z()*s;
    double R22 = axe.y()*axe.y()*(1-c) + c;
    double R23 = axe.y()*axe.z()*(1-c) - axe.x()*s;
    double R31 = axe.x()*axe.z()*(1-c) - axe.y()*s;
    double R32 = axe.y()*axe.z()*(1-c) + axe.x()*s;
    double R33 = axe.z()*axe.z()*(1-c) + c;

    return Vecteur(R11*_x+R12*_y+R13*_z, R21*_x+R22*_y+R23*_z, R31*_x+R32*_y+R33*_z);
}
