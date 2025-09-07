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

static bool invert3(Vec3 c1, Vec3 c2, Vec3 c3, Vec3 outR[3], double eps){
    double D = det3(c1,c2,c3);
    if (fabs(D) < eps) return false;
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

static inline Vec3 mulMinv(const Vec3 R[3], Vec3 p){
    return v3(R[0].x*p.x + R[0].y*p.y + R[0].z*p.z,
              R[1].x*p.x + R[1].y*p.y + R[1].z*p.z,
              R[2].x*p.x + R[2].y*p.y + R[2].z*p.z);
}

static bool clip_seg_halfplane(Vec2* p0, Vec2* p1, double a, double b, double c, double eps){
    double d0 = a*p0->x + b*p0->y + c;
    double d1 = a*p1->x + b*p1->y + c;
    bool in0 = (d0 >= -eps), in1 = (d1 >= -eps);
    if (in0 && in1) return true;
    if (!in0 && !in1) return false;
    double t = d0 / (d0 - d1);
    Vec2 q = (Vec2){ p0->x + (p1->x - p0->x)*t, p0->y + (p1->y - p0->y)*t };
    if (!in0) *p0 = q; else *p1 = q;
    return true;
}

static bool seg_vs_unit_triangle(Vec2 p0, Vec2 p1, double eps){
    if (!clip_seg_halfplane(&p0, &p1,  1, 0, 0, eps)) return false;
    if (!clip_seg_halfplane(&p0, &p1,  0, 1, 0, eps)) return false;
    if (!clip_seg_halfplane(&p0, &p1, -1,-1, 1, eps)) return false;
    return true;
}

typedef struct {
    Vec3 MinvRows[3];
    Vec3 t0;
    bool valid;
} AffineCanon;

static AffineCanon canon_from_triangle(const Tri* A, double eps){
    Vec3 e1 = sub(A->v[1], A->v[0]);
    Vec3 e2 = sub(A->v[2], A->v[0]);
    Vec3 n  = cross(e1, e2);
    AffineCanon C = {0};
    if (!invert3(e1, e2
