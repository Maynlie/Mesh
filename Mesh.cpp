# define M_PI           3.14159265358979323846  /* pi */
#include "Mesh.hpp"
#include <gl/gl.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

using namespace std;

void Mesh::translate(Vecteur& v)
{
    for(unsigned int i = 0; i < geom.size(); i++)
    {
        geom[i]+=v;
    }
}

void Mesh::rotation(Vecteur& r)
{
    for(int i = 0; i < geom.size(); i++)
    {
        //Rotation en X
        geom[i].sety(geom[i].gety()*cos(r.getx())-geom[i].getz()*sin(r.getx()));
        geom[i].setz(geom[i].gety()*sin(r.getx())+geom[i].getz()*cos(r.getx()));

        //Rotation en Y
        geom[i].setx(geom[i].getx()*cos(r.gety())+geom[i].getz()*sin(r.gety()));
        geom[i].setz(-geom[i].getx()*sin(r.gety())+geom[i].getz()*cos(r.gety()));

        //Rotation en Z
        geom[i].setx(geom[i].getx()*cos(r.getz())-geom[i].gety()*sin(r.getz()));
        geom[i].sety(geom[i].getx()*sin(r.getz())+geom[i].gety()*cos(r.getz()));
    }
}

void Mesh::draw()
{
    for(unsigned int i = 0; i < topo.size(); i+=3)
    {
        glBegin(GL_TRIANGLES);

        glColor3f(1.0f, 0.0f, 0.0f);   glVertex3f(geom[topo[i]].getx(), geom[topo[i]].gety(), geom[topo[i]].getz());
        glColor3f(0.0f, 1.0f, 0.0f);   glVertex3f(geom[topo[i+1]].getx(), geom[topo[i+1]].gety(), geom[topo[i+1]].getz());
        glColor3f(0.0f, 0.0f, 1.0f);   glVertex3f(geom[topo[i+2]].getx(), geom[topo[i+2]].gety(), geom[topo[i+2]].getz());

        glEnd();
    }
}

void Mesh::save(string titre)
{
    ofstream flux(string(titre+".obj").c_str());
    if(flux)
    {
        flux << "o " << titre << endl;

        flux << endl;

        for(unsigned int i = 0; i < geom.size(); i++)
        {
            flux << "v ";
            flux << geom[i].getx() << " ";
            flux << geom[i].gety() << " ";
            flux << geom[i].getz() << " " << endl;
        }

        flux << endl;
        for(unsigned int i = 0; i < topo.size(); i+=3)
        {
            flux << "f " << topo[i]+1 << " " << topo[i+1]+1 << " " << topo[i+2]+1 << endl;
        }
    }
}

void Mesh::cube(Vecteur& a, Vecteur& b)
{

    geom.push_back(a);
    geom.push_back(Vecteur(b.getx(), a.gety(), a.getz()));
    geom.push_back(Vecteur(a.getx(), b.gety(), a.getz()));
    geom.push_back(Vecteur(b.getx(), b.gety(), a.getz()));

    geom.push_back(Vecteur(a.getx(), a.gety(), b.getz()));
    geom.push_back(Vecteur(b.getx(), a.gety(), b.getz()));
    geom.push_back(Vecteur(a.getx(), b.gety(), b.getz()));
    geom.push_back(b);

    //Face 1
    topo.push_back(1);
    topo.push_back(3);
    topo.push_back(0);

    topo.push_back(2);
    topo.push_back(0);
    topo.push_back(3);
    //Face 2
    topo.push_back(5);
    topo.push_back(7);
    topo.push_back(1);

    topo.push_back(3);
    topo.push_back(1);
    topo.push_back(7);
    //Face 3
    topo.push_back(4);
    topo.push_back(6);
    topo.push_back(5);

    topo.push_back(7);
    topo.push_back(5);
    topo.push_back(6);
    //Face 4
    topo.push_back(0);
    topo.push_back(2);
    topo.push_back(4);

    topo.push_back(6);
    topo.push_back(4);
    topo.push_back(2);

    //Face 5
    topo.push_back(3);
    topo.push_back(7);
    topo.push_back(2);

    topo.push_back(6);
    topo.push_back(2);
    topo.push_back(7);
    //Face 6
    topo.push_back(1);
    topo.push_back(5);
    topo.push_back(0);

    topo.push_back(4);
    topo.push_back(0);
    topo.push_back(5);
}

void Mesh::TestSphere()
{
  int latitudeBands = 30;
  int longitudeBands = 30;
  int radius = 2;

  //var vertexPositionData = [];
    vector<double> normalData;
    vector<double> textureCoordData;
    for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
      double theta = latNumber * M_PI / latitudeBands;
      double sinTheta = sin(theta);
      double cosTheta = cos(theta);

      for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
        double phi = longNumber * 2 * M_PI / longitudeBands;
        double sinPhi = sin(phi);
        double cosPhi = cos(phi);

        double x = cosPhi * sinTheta;
        double y = cosTheta;
        double z = sinPhi * sinTheta;
        double u = 1 - (longNumber / longitudeBands);
        double v = 1 - (latNumber / latitudeBands);

        normalData.push_back(x);
        normalData.push_back(y);
        normalData.push_back(z);
        textureCoordData.push_back(u);
        textureCoordData.push_back(v);
        geom.push_back(Vecteur(radius * x, radius * y, radius * z));
      }
    }

    cout << geom.size() << endl;

    for (int latNumber = 0; latNumber < latitudeBands; latNumber++) {
      for (int longNumber = 0; longNumber < longitudeBands; longNumber++) {
        int first = (latNumber * (longitudeBands + 1)) + longNumber;
        int second = first + longitudeBands + 1;
        topo.push_back(first);
        topo.push_back(second);
        topo.push_back(first + 1);

        topo.push_back(second);
        topo.push_back(second + 1);
        topo.push_back(first + 1);
      }
    }
    cout << topo.size() << endl;
}

void Mesh::Cercle()
{
    int R = 5;

    geom.push_back(Vecteur(0, 0, 0));
    //sphere_vertex.push(1, 1, 0);
    double Ecart = 2*M_PI/128;

    for(int angle = 0; angle < 128; angle++)
    {

        geom.push_back(Vecteur(R*cos(Ecart*angle), R*sin(Ecart*angle), 0));
        //sphere_vertex.push(1, 1, 0);
        if(angle!=0)
        {
            topo.push_back(0);
            topo.push_back(angle);
            topo.push_back(angle+1);
        }
    }
    topo.push_back(0);
    topo.push_back(127);
    topo.push_back(1);
}

void Mesh::cylindre(Vecteur& a, Vecteur& b, float radius)
{
    int diviseur = 128;
    geom.push_back(a);
    geom.push_back(b);
    int compteur = 0;
    double angle;
    Vecteur i, j, k;
    Vecteur ab(b.getx()-a.getx(), b.gety()-a.gety(), b.getz()-a.getz());
    float l = sqrt(pow(ab.getx(), 2) + pow(ab.gety(), 2) + pow(ab.getz(), 2));
    k = ab/l;
    if(k.getx() == 0 && k.gety() == 0)
    {
        i = Vecteur(0, -k.getz(), 0);
    }
    else
    {
        i = Vecteur(-k.gety(), k.getx(), 0);
    }
    i = i/sqrt(pow(i.getx(),2) + pow(i.gety(),2) + pow(i.getz(),2));
    j = Vecteur((i.gety()*k.getz()) - (k.gety()*i.getz()), (i.getz()*k.getx()) - (k.getz()*i.getx()), (i.getx()*k.gety()) - (k.getx()*i.gety()));
    j = j/sqrt(pow(j.getx(),2) + pow(j.gety(),2) + pow(j.getz(),2));

    i*=radius;
    j*=radius;

    //cout << "a" << a << "b" << b << "ab" << ab << "i" << i << "j" << j << "k" << k << endl;
    float ecartAngle = 2*M_PI/diviseur;

    for(angle=0; angle<=diviseur; angle++)
    {
        if(angle!=diviseur)
        {
            Vecteur c = a;
            Vecteur i2 = i;
            Vecteur j2 = j;
            i2*=cos(angle*ecartAngle);
            j2*=sin(angle*ecartAngle);
            c+= i2;
            c+= j2;

            Vecteur d = b;
            d+= i2;
            d+= j2;
            //Vecteur d = b + i*cos(angle) + j*sin(angle);

            //Vecteur c = a + ((i*cos(angle) + j*sin(angle))*radius);
            //Vecteur d = b + ((i*cos(angle) + j*sin(angle))*radius);

            geom.push_back(c);
            geom.push_back(d);
            //geom.push_back(d);
        }

        if(compteur!=0)
        {
            if(angle==diviseur)
            {
                topo.push_back(0);
                topo.push_back(2);
                topo.push_back(2*compteur);

                topo.push_back(1);
                topo.push_back(2*compteur+1);
                topo.push_back(3);

                topo.push_back(2*compteur);
                topo.push_back(2*compteur+1);
                topo.push_back(3);

                topo.push_back(2);
                topo.push_back(2*compteur);
                topo.push_back(3);
            }
            else
            {
                topo.push_back(0);
                topo.push_back(2*compteur+2);
                topo.push_back(2*compteur);

                topo.push_back(1);
                topo.push_back(2*compteur+1);
                topo.push_back(2*compteur+3);

                topo.push_back(2*compteur);
                topo.push_back(2*compteur+1);
                topo.push_back(2*compteur+3);

                topo.push_back(2*compteur+2);
                topo.push_back(2*compteur);
                topo.push_back(2*compteur+3);
            }
        }

        compteur++;
    }
}

void Mesh::triangle(Vecteur& a, float base, float hauteur)
{
    Vecteur b = Vecteur(a.getx()+base, a.gety(), a.getz());
    Vecteur c = Vecteur(a.getx(), a.gety()+base, a.getz());
    Vecteur d = Vecteur(b.getx(), c.gety(), a.getz());

    geom.push_back(a);
    geom.push_back(b);
    geom.push_back(d);
    geom.push_back(c);
    geom.push_back(Vecteur((b.getx()-a.getx())/2+a.getx(), (c.gety()-a.gety())/2+a.gety(), a.getz()+hauteur));

    topo.push_back(0);
    topo.push_back(1);
    topo.push_back(2);

    topo.push_back(0);
    topo.push_back(2);
    topo.push_back(3);

    topo.push_back(0);
    topo.push_back(1);
    topo.push_back(4);

    topo.push_back(1);
    topo.push_back(2);
    topo.push_back(4);

    topo.push_back(2);
    topo.push_back(3);
    topo.push_back(4);

    topo.push_back(3);
    topo.push_back(0);
    topo.push_back(4);

}


void Mesh::sphere(Vecteur& a, float radius)
{
    std::vector<Vecteur> listPoint;

    int diviseur = 128;
    geom.push_back(a);
    int compteur = 0;
    double angle;
    Vecteur i, j, k;
    float l = sqrt(pow(a.getx(), 2) + pow(a.gety(), 2) + pow(a.getz(), 2));
    k = Vecteur(1,1,1);
    i = Vecteur(-k.gety(), k.getx(), 0);
    i = i/sqrt(pow(i.getx(),2) + pow(i.gety(),2) + pow(i.getz(),2));
    j = Vecteur((i.gety()*k.getz()) - (k.gety()*i.getz()), (i.getz()*k.getx()) - (k.getz()*i.getx()), (i.getx()*k.gety()) - (k.getx()*i.gety()));
    j = j/sqrt(pow(j.getx(),2) + pow(j.gety(),2) + pow(j.getz(),2));

    i*=radius;
    j*=radius;

    float ecartAngle = 2*M_PI/diviseur;

    for(angle=0; angle<=diviseur/2-1; angle++)
    {
        Vecteur c = a;
        Vecteur i2 = i;
        Vecteur j2 = j;
        i2*=cos(angle*ecartAngle);
        j2*=sin(angle*ecartAngle);
        c+= i2;
        c+= j2;

        listPoint.push_back(c);

    }

    for(int x=0; x<listPoint.size(); x++)
    {

        Vecteur b = listPoint[x];

        compteur = 0;
        Vecteur ab(b.getx()-a.getx(), b.gety()-a.gety(), b.getz()-a.getz());
        float l = sqrt(pow(ab.getx(), 2) + pow(ab.gety(), 2) + pow(ab.getz(), 2));
        k = ab/l;
        i = Vecteur(-k.gety(), k.getx(), 0);
        i = i/sqrt(pow(i.getx(),2) + pow(i.gety(),2) + pow(i.getz(),2));
        j = Vecteur((i.gety()*k.getz()) - (k.gety()*i.getz()), (i.getz()*k.getx()) - (k.getz()*i.getx()), (i.getx()*k.gety()) - (k.getx()*i.gety()));
        j = j/sqrt(pow(j.getx(),2) + pow(j.gety(),2) + pow(j.getz(),2));

        i*=radius;
        j*=radius;

        for(angle=0; angle<=diviseur; angle++)
        {
            if(angle!=diviseur)
            {
                Vecteur c = a;
                Vecteur i2 = i;
                Vecteur j2 = j;
                i2*=cos(angle*ecartAngle);
                j2*=sin(angle*ecartAngle);
                c+= i2;
                c+= j2;

                geom.push_back(c);
            }

            if(x>0 && angle!=diviseur && compteur)
            {
                topo.push_back((x*(diviseur))+compteur);
                topo.push_back(((x-1)*(diviseur))+compteur);
                topo.push_back(((x-1)*(diviseur))+compteur+1);

                topo.push_back((x*(diviseur))+compteur-1);
                topo.push_back(((x-1)*(diviseur))+compteur);
                topo.push_back((x*(diviseur))+compteur);
            }

            if(x==listPoint.size()-1)
            {
                topo.push_back(((x)*(diviseur))+compteur);
                topo.push_back(((x)*(diviseur))+compteur-1);
                topo.push_back((diviseur/2+diviseur-compteur+2)%diviseur);

                topo.push_back(((x)*(diviseur))+compteur);
                topo.push_back((diviseur/2+diviseur-compteur+2)%diviseur);
                topo.push_back((diviseur/2+diviseur-compteur+1)%diviseur);
            }

            compteur++;
        }
    }
}

void noise(float h, float p, Vecteur v)
{

}

void Mesh::toit(Vecteur& a, Vecteur& b)
{
    cube( a,  b);
    Vecteur c = Vecteur((b.getx()-a.getx())/2+a.getx(), a.gety(), b.getz());
    Vecteur temp = Vecteur(0,0,(b.getz() - a.getz()));
    c += temp;
    geom.push_back(c);
    Vecteur d = Vecteur((b.getx()-a.getx())/2+a.getx(), b.gety(), b.getz());
    d += temp;
    geom.push_back(d);

    topo.push_back(4);
    topo.push_back(5);
    topo.push_back(8);

    topo.push_back(7);
    topo.push_back(6);
    topo.push_back(9);

    topo.push_back(4);
    topo.push_back(8);
    topo.push_back(6);

    topo.push_back(6);
    topo.push_back(8);
    topo.push_back(9);

    topo.push_back(5);
    topo.push_back(8);
    topo.push_back(7);

    topo.push_back(7);
    topo.push_back(8);
    topo.push_back(9);
}

void Mesh::countTriangles()
{
    cout << topo.size()/3 << endl;
}

void Mesh::temple(Vecteur& a, Vecteur& b, float h)
{
    Mesh base;
    base.cube(a, b);
    fusionMesh(base);
    base.Clear();
    float longueur = a.gety() - b.gety();
    if(longueur<0) longueur = b.gety() - a.gety();
    float largeur = a.getx() - b.getx();
    if(largeur<0) largeur = b.getx() - a.getx();
    largeur-=2.0f;
    Vecteur c, d;
    c = a;
    c.setz(b.getz());
    c.setx(a.getx()+1);
    c.sety(a.gety()+1);
    d = c;
    Vecteur hauteur(0, 0, h);
    Vecteur ecart(0, 3, 0);
    Vecteur large(largeur, 0, 0);
    d += hauteur;
    for(float i = 0.0f; i <= longueur-0.25f; i+=3.0f)
    {
        if(i!=0)
        {
            c+=ecart;
            d+=ecart;
        }

        base.collone(c, d, 0.5f);
        fusionMesh(base);
        base.Clear();
        c+=large;
        d+=large;
        base.collone(c, d, 0.5f);
        fusionMesh(base);
        base.Clear();
        c-=large;
        d-=large;
    }
    hauteur.setz(h+(b.getz()-a.getz()));
    a+=hauteur;
    b+=hauteur;
    base.toit(a, b);
    fusionMesh(base);
    base.Clear();
}

void Mesh::collone(Vecteur& a, Vecteur& b, float radius)
{
    Mesh addMesh;
    cylindre(a, b, radius);
    Vecteur ab(b.getx()-a.getx(), b.gety()-a.gety(), b.getz()-a.getz());
    ab = ab/sqrt(pow(ab.getx(), 2) + pow(ab.gety(), 2) + pow(ab.getz(), 2));
    addMesh.tor(a, ab, radius);
    fusionMesh(addMesh);
    addMesh.Clear();
    addMesh.tor(b, ab, radius);
    fusionMesh(addMesh);
}

void Mesh::tor(Vecteur& a, Vecteur& normaleT, float radius)
{

    std::vector<Vecteur> listPoint;

    int diviseur = 128;
    geom.push_back(a);
    int compteur = 0;
    double angle;
    Vecteur i, j, k;
    float l = sqrt(pow(a.getx(), 2) + pow(a.gety(), 2) + pow(a.getz(), 2));
    k = normaleT;
    if(k.getx() == 0 && k.gety() == 0)
    {
        i = Vecteur(0, -k.getz(), 0);
    }
    else
    {
        i = Vecteur(-k.gety(), k.getx(), 0);
    }
    normaleT = normaleT/sqrt(pow(normaleT.getx(),2) + pow(normaleT.gety(),2) + pow(normaleT.getz(),2));
    i = i/sqrt(pow(i.getx(),2) + pow(i.gety(),2) + pow(i.getz(),2));
    j = Vecteur((i.gety()*k.getz()) - (k.gety()*i.getz()), (i.getz()*k.getx()) - (k.getz()*i.getx()), (i.getx()*k.gety()) - (k.getx()*i.gety()));
    j = j/sqrt(pow(j.getx(),2) + pow(j.gety(),2) + pow(j.getz(),2));

    i*=radius;
    j*=radius;
    normaleT*=radius/2;

    float ecartAngle = 2*M_PI/diviseur;

    for(angle=0; angle<=diviseur; angle++)
    {
        Vecteur c = a;
        Vecteur i2 = i;
        Vecteur j2 = j;
        i2*=cos(angle*ecartAngle);
        j2*=sin(angle*ecartAngle);
        c+= i2;
        c+= j2;

        listPoint.push_back(c);

    }

    for(int x=0; x<listPoint.size(); x++)
    {
        Vecteur b = listPoint[x];

        compteur = 0;
        Vecteur ab((b.getx()-a.getx()), (b.gety()-a.gety()), (b.getz()-a.getz()));
        float l = sqrt(pow(ab.getx(), 2) + pow(ab.gety(), 2) + pow(ab.getz(), 2));
        k = ab/l;

        ab*=radius/2;

        for(angle=0; angle<=diviseur; angle++)
        {
            if(angle!=diviseur)
            {
                Vecteur c = b;
                Vecteur ab2 = ab;
                Vecteur normaleT2 = normaleT;
                ab2*=cos(angle*ecartAngle);
                normaleT2*=sin(angle*ecartAngle);
                c+= ab2;
                c+= normaleT2;

                geom.push_back(c);
            }

            if(x>0)
            {
                if(angle==diviseur)
                {
                    topo.push_back(((x-1)*(diviseur))+compteur);
                    topo.push_back(((x-1)*(diviseur))+1);
                    topo.push_back((x*(diviseur))+1);

                    topo.push_back(((x-1)*(diviseur))+compteur);
                    topo.push_back((x*(diviseur))+1);
                    topo.push_back((x*(diviseur))+compteur);
                }
                else if(x!=1 && compteur!=0)
                {
                    topo.push_back(((x-1)*(diviseur))+compteur);
                    topo.push_back(((x-1)*(diviseur))+compteur+1);
                    topo.push_back((x*(diviseur))+compteur+1);

                    topo.push_back(((x-1)*(diviseur))+compteur);
                    topo.push_back((x*(diviseur))+compteur+1);
                    topo.push_back((x*(diviseur))+compteur);
                }

            }

            if(x==0 && compteur!=0)
            {
                topo.push_back(compteur);
                topo.push_back(compteur+1);
                topo.push_back(diviseur+compteur+1);

                topo.push_back(compteur);
                topo.push_back(diviseur+compteur+1);
                topo.push_back(diviseur+compteur);
            }

            compteur++;
        }
    }
}

void Mesh::fusionMesh(Mesh& mesh)
{
    int taille = geom.size();
    for(int i=0; i<mesh.geom.size(); i++){
        geom.push_back(mesh.geom[i]);
    }

    for(int i=0; i<mesh.topo.size(); i++){
        topo.push_back(mesh.topo[i]+taille);
    }
}

void Mesh::Clear()
{
    geom.clear();
    topo.clear();
}
