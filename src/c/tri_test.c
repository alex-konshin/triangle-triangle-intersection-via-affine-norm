// src/c/tri_test.c
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct { double x,y,z; } Vec3;
typedef struct { Vec3 v[3]; } Tri;

// forward decl from tri_tri_affine.c
bool tri_tri_intersect_affine(const Tri* A, const Tri* B, double eps);

static Vec3 V(double x,double y,double z){ Vec3 r={x,y,z}; return r; }
static Tri T(Vec3 a,Vec3 b,Vec3 c){ Tri t={{a,b,c}}; return t; }

int main(int argc, char** argv){
    double eps = (argc>1)? atof(argv[1]) : 1e-9;

    // A: unit right triangle in xy
    Tri A = T(V(0,0,0), V(1,0,0), V(0,1,0));

    // B1: intersects (straddles z)
    Tri B1 = T(V(0.2,0.2,-1), V(0.8,0.2, 1), V(0.2,0.8,-1));
    printf("B1 intersect: %s\n", tri_tri_intersect_affine(&A,&B1,eps) ? "YES":"NO");

    // B2: no intersection (all below)
    Tri B2 = T(V(2,2,-1), V(3,2,-1), V(2,3,-1));
    printf("B2 intersect: %s\n", tri_tri_intersect_affine(&A,&B2,eps) ? "YES":"NO");

    // B3: touch at a point
    Tri B3 = T(V(0.5,0.5,-1), V(0.5,0.5,1), V(0.6,0.6,1));
    printf("B3 intersect: %s\n", tri_tri_intersect_affine(&A,&B3,eps) ? "YES":"NO");

    return 0;
}
