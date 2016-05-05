#include "Vecteur.hpp"


Vecteur::Vecteur()
{
    x = 0;
    y = 0;
    z = 0;
}

Vecteur::Vecteur(const Vecteur& v)
{
    x = v.getx();
    y = v.gety();
    z = v.getz();
}

Vecteur::Vecteur(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}
void Vecteur::operator+=(Vecteur& v)
{
    x+=v.x;
    y+=v.y;
    z+=v.z;
}

void Vecteur::operator-=(Vecteur& v)
{
    x-=v.x;
    y-=v.y;
    z-=v.z;
}

void Vecteur::operator*=(float f)
{
    x*=f;
    y*=f;
    z*=f;
}

void Vecteur::operator*(float f)
{
    (*this)*=f;
    /*x=x*f;
    y=y*f;
    z=z*f;
    return *this;*/
}

Vecteur& Vecteur::operator/(float l)
{
    x=x/l;
    y=y/l;
    z=z/l;
    return *this;
}

void Vecteur::operator+(Vecteur& v)
{
    (*this)+=v;
    /*x+=v.x;
    y+=v.y;
    z+=v.z;
    return *this;*/
}

std::ostream& operator<<(std::ostream& out, Vecteur& v)
{
    out << "(" << v.getx() << ";" << v.gety() << ";" << v.getz() <<")" << std::endl;
    return out;
}
