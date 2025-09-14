// src/c/tri_tri_affine.c
// CC0 1.0 — Public Domain

// This C implementation is provided for demonstration purposes only.
// It illustrates the core idea of the algorithm and is not intended
// to cover all possible edge or degenerate cases.
// It is also not optimized for high performance.

#include <math.h>
#include <stdbool.h>

typedef struct { double x,y,z; } Vec3;
typedef struct { double x,y;   } Vec2;
typedef struct { Vec3 v[3];    } Tri;

static inline Vec3 v3(double x,double y,double z){ return (Vec3){x,y,z}; }
static inline Vec3 add(Vec3 a, Vec3 b){ return v3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 sub(Vec3 a, Vec3 b){ return v3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline Vec3 muls(Vec3 a, double s){ return v3(a.x*s, a.y*s, a.z*s); }
static inline Vec3 lerp3(Vec3 a, Vec3 b, double t){ return add(muls(a,1.0-t), muls(b,t)); }
static inline Vec3 cross(Vec3 a, Vec3 b){
    return v3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static inline double det3(Vec3 c1, Vec3 c2, Vec3 c3){
    return c1.x*(c2.y*c3.z - c2.z*c3.y)
         - c2.x*(c1.y*c3.z - c1.z*c3.y)
         + c3.x*(c1.y*c2.z - c1.z*c2.y);
}

// Inverse of M=[c1 c2 c3] using adjugate/determinant. Returns false if near singular.
static bool invert3(Vec3 c1, Vec3 c2, Vec3 c3, Vec3 outR[3], double eps){
    double D = det3(c1,c2,c3);
    if (fabs(D) < eps) return false;
    // Adjugate rows: (cofactors transposed)
    Vec3 r0 = v3(  (c2.y*c3.z - c2.z*c3.y),
                  -(c1.y*c3.z - c1.z*c3.y),
                   (c1.y*c2.z - c1.z*c2.y));
    Vec3 r1 = v3( -(c2.x*c3.z - c2.z*c3.x),
                   (c1.x*c3.z - c1.z*c3.x),
                  -(c1.x*c2.z - c1.z*c2.x));
    Vec3 r2 = v3(  (c2.x*c3.y - c2.y*c3.x),
                  -(c1.x*c3.y - c1.y*c3.x),
                   (c1.x*c2.y - c1.y*c2.x));
    double invD = 1.0 / D;
    outR[0] = muls(r0, invD);
    outR[1] = muls(r1, invD);
    outR[2] = muls(r2, invD);
    return true;
}

// Multiply row-vector form (rows of Minv) by column vector p: Minv * p
static inline Vec3 mulMinv(const Vec3 R[3], Vec3 p){
    return v3(R[0].x*p.x + R[0].y*p.y + R[0].z*p.z,
              R[1].x*p.x + R[1].y*p.y + R[1].z*p.z,
              R[2].x*p.x + R[2].y*p.y + R[2].z*p.z);
}

// Clip segment against half-plane a*x + b*y + c >= 0 (z=0 plane already enforced)
static bool clip_seg_halfplane(Vec2* p0, Vec2* p1, double a, double b, double c, double eps){
    double d0 = a*p0->x + b*p0->y + c;
    double d1 = a*p1->x + b*p1->y + c;
    bool in0 = (d0 >= -eps), in1 = (d1 >= -eps);
    if (in0 && in1) return true;          // unchanged
    if (!in0 && !in1) return false;       // fully out
    double t = d0 / (d0 - d1);            // intersection param
    Vec2 q = (Vec2){ p0->x + (p1->x - p0->x)*t, p0->y + (p1->y - p0->y)*t };
    if (!in0) *p0 = q; else *p1 = q;
    return true;
}

// Segment vs unit right triangle: x>=0, y>=0, x+y<=1
static bool seg_vs_unit_triangle(Vec2 p0, Vec2 p1, double eps){
    if (!clip_seg_halfplane(&p0, &p1,  1, 0, 0, eps)) return false; // x>=0
    if (!clip_seg_halfplane(&p0, &p1,  0, 1, 0, eps)) return false; // y>=0
    if (!clip_seg_halfplane(&p0, &p1, -1,-1, 1, eps)) return false; // x+y<=1
    return true;
}

typedef struct {
    Vec3 MinvRows[3]; // rows of Minv
    Vec3 t0;          // -Minv * a0
    bool valid;
} AffineCanon;

// Precompute canonical mapping for A
static AffineCanon canon_from_triangle(const Tri* A, double eps){
    Vec3 e1 = sub(A->v[1], A->v[0]);
    Vec3 e2 = sub(A->v[2], A->v[0]);
    Vec3 n  = cross(e1, e2);
    AffineCanon C = {0};
    if (!invert3(e1, e2, n, C.MinvRows, eps)) { C.valid = false; return C; }
    // t0 = -Minv * a0
    Vec3 a0p = mulMinv(C.MinvRows, A->v[0]);
    C.t0 = muls(a0p, -1.0);
    C.valid = true;
    return C;
}

// Transform point: b' = Minv * b + t0
static inline Vec3 to_canon(const AffineCanon* C, Vec3 p){
    return add(mulMinv(C->MinvRows, p), C->t0);
}

bool tri_tri_intersect_affine(const Tri* A, const Tri* B, double eps){
    AffineCanon C = canon_from_triangle(A, eps);
    if (!C.valid) return false; // A degenerate

    Vec3 bp[3];
    int pos=0, neg=0, cop=0;
    for (int i=0;i<3;i++){
        bp[i] = to_canon(&C, B->v[i]);  // (x,y,z)
        if      (bp[i].z >  eps) pos++;
        else if (bp[i].z < -eps) neg++;
        else                     cop++;
    }
    if (pos==3 || neg==3) return false;

    // Coplanar: 2D triangle overlap in (x,y) — implement as needed
    if (pos+neg==0){
        // Minimal fallback: check edges vs edges + vertex inclusion (left as exercise)
        // For production: 2D SAT or segment-segment tests.
        // Returning true here if projections overlap can be overly permissive.
        // Implement tri2d_overlap_xy(bp, eps) for robustness.
        // Placeholder conservative check (not sufficient for all coplanar cases):
        return true;
    }

    // Collect intersections of B edges with z=0
    Vec2 seg[2]; int k=0;
    for (int i=0;i<3;i++){
        int j=(i+1)%3;
        double zi = bp[i].z, zj = bp[j].z;
        if ((zi > eps && zj < -eps) || (zi < -eps && zj > eps)) {
            double t = zi / (zi - zj);
            Vec3 P = lerp3(bp[i], bp[j], t);
            if (k < 2) seg[k++] = (Vec2){P.x, P.y};
        } else if (fabs(zi) <= eps && k < 2) {
            seg[k++] = (Vec2){bp[i].x, bp[i].y};
        }
        if (k==2) break;
    }
    if (k==0) return false;
    if (k==1) seg[1] = seg[0];

    return seg_vs_unit_triangle(seg[0], seg[1], eps);
}

// Optional: map canonical z=0 endpoints back to world space
static inline Vec3 from_canon_point(Vec2 pxy, Vec3 c1, Vec3 c2, Vec3 n, Vec3 a0){
    // p = M * (x,y,0) + a0 = x*c1 + y*c2 + a0
    return add(add(muls(c1, pxy.x), muls(c2, pxy.y)), a0);
}
