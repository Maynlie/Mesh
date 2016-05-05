#ifndef MESH
#define MESH

#include <vector>
#include <string>
#include "Vecteur.hpp"

class Mesh{

private:
    std::vector<Vecteur> geom;
    std::vector<int> topo;
public:
    void translate(Vecteur& v);
    void rotation(Vecteur& r);
    void cube(Vecteur& a, Vecteur& b);
    void Cercle();
    void TestSphere();
    void cylindre(Vecteur& a, Vecteur& b, float radius);
    void toit(Vecteur& a, Vecteur& b);
    void sphere(Vecteur& a, float radius);
    void Sphere2(Vecteur& a, float radius);
    void triangle(Vecteur& a, float base, float hauteur);
    void collone(Vecteur& a, Vecteur& b, float radius);
    void tor(Vecteur& a, Vecteur& normaleT, float radius);
    void temple(Vecteur& a, Vecteur& b, float h);
    void fusionMesh(Mesh& mesh);
    void Noise(float h, float p);
    void countTriangles();
    void Clear();
    void draw();
    void save(std::string titre = "Mesh");
};

#endif // MESH
