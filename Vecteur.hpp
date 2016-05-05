#ifndef VECTEUR
#define VECTEUR

#include <iostream>

class Vecteur{
private:
    float x;
    float y;
    float z;
public:
    Vecteur();
    Vecteur(const Vecteur& v);
    Vecteur(float x, float y, float z);
    void operator+=(Vecteur& v);
    void operator-=(Vecteur& v);
    friend std::ostream& operator<<(std::ostream& out, Vecteur& v);
    Vecteur& operator/(float l);
    void operator*(float f);
    void operator*=(float f);
    void operator+(Vecteur& v);
    float getx() const {return x;};
    float gety() const {return y;};
    float getz() const {return z;};
    void setx(float x){this->x = x;};
    void sety(float y){this->y = y;};
    void setz(float z){this->z = z;};
};
#endif // VECTEUR
