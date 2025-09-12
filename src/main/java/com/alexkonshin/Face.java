package com.alexkonshin;

/**
 * Reference implementation of triangle–triangle intersection via affine normal basis.
 * Core idea:
 * - Build local coords for triangle A: q1 at origin; X along (q2-q1); Y along (q3-q1); Z along normal.
 * - In this basis, A lies in plane z=0 within the simplex: x>=0, y>=0, x+y<=1.
 * - Transform triangle B into this basis.
 * - If all B vertices are strictly on one side of z=0 => no intersection.
 * - Compute line of plane intersection in z=0 (rnx*x + rny*y = d), clip it by A's simplex to get [p0,p1] in s=x+y.
 * - Find intersection of B's edges with z=0, get two points and their s=x+y as [s0,s1].
 * - Intervals overlap iff triangles intersect.
 */
public final class Face {

    // Geometry
    private final Vertex[] v = new Vertex[3];

    public Face(Vertex a, Vertex b, Vertex c) {
        v[0] = a; v[1] = b; v[2] = c;
    }

    // Cached AABB [xmin, xmax, ymin, ymax, zmin, zmax]
    private float[] minmax; // computed lazily

    // Normal of this face in world space
    private float nx = Float.NaN, ny = Float.NaN, nz = Float.NaN;

    // World->local transform matrix (columns: e1, e2, n), we cache its inverse rows r1*, r2*, r3*
    private float r11, r12, r13;
    private float r21, r22, r23;
    private float r31, r32, r33 = Float.NaN; // NaN = not computed yet

    // Numerical tolerances
    private static final float EPS = 1e-6f;
    private static final float PLANE_EPS = 1e-6f;

    // --- Public API

    public boolean intersects(Face other) {
        // 1) Fast AABB reject
        float[] a = getMinMax();
        float[] b = other.getMinMax();
        if (b[1] < a[0] || b[0] > a[1] ||
            b[3] < a[2] || b[2] > a[3] ||
            b[5] < a[4] || b[4] > a[5]) {
            return false;
        }

        // 2) Ensure transform of this face is ready (this face becomes local XY; its normal is +Z)
        ensureTransform();

        // 3) Transform other’s vertices into this local space
        Vec3 q1 = v[0].toVec3();
        Vec3 b1 = worldToLocal(other.v[0].toVec3().sub(q1));
        Vec3 b2 = worldToLocal(other.v[1].toVec3().sub(q1));
        Vec3 b3 = worldToLocal(other.v[2].toVec3().sub(q1));

        // 4) Quick plane-side test for z=0
        int pos = (b1.z > PLANE_EPS ? 1 : 0) + (b2.z > PLANE_EPS ? 1 : 0) + (b3.z > PLANE_EPS ? 1 : 0);
        int neg = (b1.z < -PLANE_EPS ? 1 : 0) + (b2.z < -PLANE_EPS ? 1 : 0) + (b3.z < -PLANE_EPS ? 1 : 0);
        if (pos == 3 || neg == 3) return false; // entirely on one side of plane z=0

        // 5) Compute normal of other triangle in local space (only x,y matter in z=0)
        Vec3 e12 = b2.sub(b1);
        Vec3 e13 = b3.sub(b1);
        Vec3 nB = e12.cross(e13); // components (rnx, rny, rnz) in local coordinates

        // If nB is near-zero, other face is degenerate (collinear or tiny)
        if (nB.length() < EPS) {
            // Fallback: intersect triangle A with segment/point of B projected onto z=0
            return degenerateTriangleIntersectInLocal(b1, b2, b3);
        }

        // 6) Line of intersection with z=0 plane: rnx*x + rny*y + rnz*0 = d, where d = nB dot b1
        float rnx = nB.x, rny = nB.y, rnz = nB.z;
        float d = rnx * b1.x + rny * b1.y + rnz * b1.z;

        // Guard divides; if rnx or rny ~ 0, choose the other intercept robustly
        boolean hasX0 = Math.abs(rnx) > EPS;
        boolean hasY0 = Math.abs(rny) > EPS;

        if (!hasX0 && !hasY0) {
            // Line is parallel to both x and y axes in z=0 → either no intersection or entire plane; fall back
            return coplanarIntersectInLocal(b1, b2, b3);
        }

        float x0 = hasX0 ? (d / rnx) : Float.POSITIVE_INFINITY; // where y=0
        float y0 = hasY0 ? (d / rny) : Float.POSITIVE_INFINITY; // where x=0

        // 7) Clip that line by triangle A’s simplex: x>=0, y>=0, x+y<=1
        // Represent chord within A by the range s = x + y in [p0, p1]
        float p0, p1;
        if (x0 >= 0f && x0 <= 1f) { // TODO compare with EPS
            if (y0 <= 0f || y0 > 1f) { // Side effect for the special case x0 = y0 = 0 => p0=0 & p1=1
                p0 = x0;
                p1 = 1f;
            } else if (x0 > y0) {
                p0 = y0;
                p1 = x0;
            } else {
                p0 = x0;
                p1 = y0;
            }
        } else if (y0 >= 0f && y0 <= 1f) {
            p0 = y0;
            p1 = 1f;
        } else {
            // The chord with A might hit the diagonal x+y=1 twice or miss; use robust clipping
            float[] chord = clipLineAgainstSimplex(rnx, rny, d);
            if (chord == null) return false;
            p0 = chord[0];
            p1 = chord[1];
        }

        // 8) Intersect other triangle edges with z=0, collect two points, compute s=x+y for them
        float[] sB = intersectTriangleEdgesWithPlaneZ0(b1, b2, b3);
        if (sB == null) {
            // Either coplanar with z=0 or no intersections within edges
            return coplanarIntersectInLocal(b1, b2, b3);
        }
        // sB is already sorted (s0 <= s1)
        float s0 = sB[0];
        float s1 = sB[1];

        // 9) Overlap test of intervals [s0, s1] vs [p0, p1]
        if (s1 < p0 - EPS || s0 > p1 + EPS) return false;
        return true;
    }

    // --- Internal helpers

    private float[] getMinMax() {
        if (minmax != null) return minmax;
        float xmin = Math.min(v[0].x, Math.min(v[1].x, v[2].x));
        float xmax = Math.max(v[0].x, Math.max(v[1].x, v[2].x));
        float ymin = Math.min(v[0].y, Math.min(v[1].y, v[2].y));
        float ymax = Math.max(v[0].y, Math.max(v[1].y, v[2].y));
        float zmin = Math.min(v[0].z, Math.min(v[1].z, v[2].z));
        float zmax = Math.max(v[0].z, Math.max(v[1].z, v[2].z));
        minmax = new float[]{xmin, xmax, ymin, ymax, zmin, zmax};
        return minmax;
    }

    private void ensureTransform() {
        if (!Float.isNaN(r33)) return; // already computed
        if (Float.isNaN(nx)) computeNormal();

        Vertex q1 = v[0];
        Vertex q2 = v[1];
        Vertex q3 = v[2];

        float a11 = q2.x - q1.x;
        float a21 = q2.y - q1.y;
        float a31 = q2.z - q1.z;
        float a12 = q3.x - q1.x;
        float a22 = q3.y - q1.y;
        float a32 = q3.z - q1.z;

        // Matrix M = [e1 e2 n], want R = M^{-1}, but we store rows r1*, r2*, r3* directly.
        float det = a11 * a22 * nz + a12 * ny * a31 + nx * a21 * a32
                  - nx * a22 * a31 - a12 * a21 * nz - a11 * ny * a32;

        if (Math.abs(det) < EPS)
            throw new IllegalStateException("Degenerate basis for transform; triangle is nearly collinear.");

        r11 = (a22 * nz - ny * a32) / det;
        r21 = (ny * a31 - a21 * nz) / det;
        r31 = (a21 * a32 - a22 * a31) / det;

        r12 = (a32 * nx - nz * a12) / det;
        r22 = (nz * a11 - a31 * nx) / det;
        r32 = (a12 * a31 - a11 * a32) / det;

        r13 = (a12 * ny - nx * a22) / det;
        r23 = (nx * a21 - a11 * ny) / det;
        r33 = (a11 * a22 - a12 * a21) / det;
    }

    private void computeNormal() {
        Vec3 e1 = v[1].toVec3().sub(v[0].toVec3());
        Vec3 e2 = v[2].toVec3().sub(v[0].toVec3());
        Vec3 n = e1.cross(e2);
        float len = n.length();
        if (len < EPS) {
            throw new IllegalStateException("Degenerate triangle: normal is undefined.");
        }
        n = n.scale(1f / len);
        nx = n.x; ny = n.y; nz = n.z;
    }

    private Vec3 worldToLocal(Vec3 d) {
        // Apply rows of inverse matrix to vector (point relative to q1)
        return new Vec3(
            r11 * d.x + r12 * d.y + r13 * d.z,
            r21 * d.x + r22 * d.y + r23 * d.z,
            r31 * d.x + r32 * d.y + r33 * d.z
        );
    }

    // Robust clipping of the line rnx*x + rny*y = d against simplex x>=0, y>=0, x+y<=1, returned as s=x+y interval [p0,p1]
    private float[] clipLineAgainstSimplex(float rnx, float rny, float d) {
        // Intersections with edges:
        // y=0 => x = d/rnx
        // x=0 => y = d/rny
        // x+y=1 => substitute y=1-x => (rnx - rny)*x + rny = d => x = (d - rny) / (rnx - rny)
      
        float s0 = Float.NaN;
        float s1 = Float.NaN;
        
        int count = 0;
        
        // candidate at (x0,0)
        if (Math.abs(rnx) > EPS) {
            float x0 = d / rnx;
            if (x0 >= -EPS && x0 <= 1f + EPS) {
              count = 1;
              s0 = x0;
            }
        }
        // candidate at (0,y0)
        if (Math.abs(rny) > EPS) {
            float y0 = d / rny;
            if (y0 >= -EPS && y0 <= 1f + EPS) {
              if (count == 0) 
                s0 = y0;
              else
                s1 = y0;
            }
        }
        if (count == 0) return null; // no intersections. 
        
        if (count == 2 && s0 > s1) {
            // swap s1 and s2
            float t = s0;
            s0 = s1;
            s1 = t;
        }
        
        // candidate at x+y=1
        if (Math.abs(rnx - rny) > EPS) {
            float x = (d - rny) / (rnx - rny);
            if (x > -EPS && x < 1f + EPS) {
                // add x + y = 1 to result
                s1 = 1f;
                count++;
            }
        }
        if (count < 2) s1 = s0; // touching single point        
        return new float[] {s0, s1};
    }

    // Intersect the three edges of triangle B with plane z=0 in local space; return two s=x+y values
    private float[] intersectTriangleEdgesWithPlaneZ0(Vec3 b1, Vec3 b2, Vec3 b3) {
        int count = 0;
        float[] result = new float[3];
        count = maybeAppendIntersectionS(b1, b2, result, count);
        count = maybeAppendIntersectionS(b1, b3, result, count);
        count = maybeAppendIntersectionS(b2, b3, result, count);
        switch(count) {
        
        case 0: return null;
        
        case 1: {
          // single point touch 
          result[1] = result[0];
          break; 
        }
        
        case 2: {
          float s0 = result[0];
          float s1 = result[1];
          if (s0 > s1) { // swap
            result[0] = s1;
            result[1] = s0;
          }
          break;
        }
        
        default: // count > 2  ==> Multiple intersections due to coplanarity or tangency; keep extreme values
          float min = result[0];
          float max = min;
          
          float t = result[1];
          if (t < min) min = t;
          else if (t > max) max = t;
          
          t = result[2];
          if (t < min) min = t;
          else if (t > max) max = t;
          
          result[0] = min;
          result[1] = max;
        }
        return result;
    }
    
    private int maybeAppendIntersectionS(Vec3 a, Vec3 b, float[] result, int count) {
        final float za = a.z, zb = b.z;

        // If both endpoints are very close to plane, treat as coplanar edge (skip here, handle elsewhere)
        if (Math.abs(za) <= PLANE_EPS && Math.abs(zb) <= PLANE_EPS) return count;

        // If signs differ or one is on plane, compute intersection
        if ((za > PLANE_EPS && zb < -PLANE_EPS) || (za < -PLANE_EPS && zb > PLANE_EPS) ||
            Math.abs(za) <= PLANE_EPS || Math.abs(zb) <= PLANE_EPS) {

            float t;
            if (Math.abs(za - zb) > EPS) {
                t = za / (za - zb); // intersection parameter from a to b where z=0
            } else {
                t = 0f; // fallback, nearly parallel in z; treat as touching at a
            }
            if (t >= -EPS && t <= 1f + EPS) {
                float x = a.x + t * (b.x - a.x);
                float y = a.y + t * (b.y - a.y);
                result[count++] = x + y;
            }
        }
        return count;
    }

    // Coplanar case in local coords (z ~ 0 for all vertices or edges) → 2D triangle overlap in (x,y) with simplex constraint of A
    private boolean coplanarIntersectInLocal(Vec3 b1, Vec3 b2, Vec3 b3) {
        // Treat A as triangle with corners (0,0), (1,0), (0,1).
        Vec2[] A = new Vec2[]{new Vec2(0,0), new Vec2(1,0), new Vec2(0,1)};
        Vec2[] B = new Vec2[]{new Vec2(b1.x, b1.y), new Vec2(b2.x, b2.y), new Vec2(b3.x, b3.y)};
        return trianglesOverlap2D(A, B);
    }

    // Degenerate B (nearly collinear) → intersect its segment hull with A’s simplex
    private boolean degenerateTriangleIntersectInLocal(Vec3 b1, Vec3 b2, Vec3 b3) {
        // Project to 2D, take convex hull (segment/point), test against A in 2D
        Vec2[] A = new Vec2[]{new Vec2(0,0), new Vec2(1,0), new Vec2(0,1)};
        Vec2[] B = segmentHull2D(new Vec2(b1.x, b1.y), new Vec2(b2.x, b2.y), new Vec2(b3.x, b3.y));
        return polygonTriangleOverlap2D(B, A);
    }

    // --- Minimal 2D helpers for coplanar handling

    private boolean trianglesOverlap2D(Vec2[] T1, Vec2[] T2) {
        // SAT on edges of both triangles
        if (!polygonsOverlapSAT(T1, T2)) return false;
        return true;
    }

    private boolean polygonTriangleOverlap2D(Vec2[] P, Vec2[] T) {
        return polygonsOverlapSAT(P, T);
    }

    private boolean polygonsOverlapSAT(Vec2[] A, Vec2[] B) {
        return !hasSeparatingAxis(A, B) && !hasSeparatingAxis(B, A);
    }

    private boolean hasSeparatingAxis(Vec2[] A, Vec2[] B) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            Vec2 p = A[i];
            Vec2 q = A[(i + 1) % n];
            Vec2 edge = q.sub(p);
            Vec2 axis = new Vec2(-edge.y, edge.x).normalize(EPS);
            float[] aProj = projectPolygon(A, axis);
            float[] bProj = projectPolygon(B, axis);
            if (aProj[1] < bProj[0] - EPS || bProj[1] < aProj[0] - EPS) return true;
        }
        return false;
    }

    private float[] projectPolygon(Vec2[] P, Vec2 axis) {
        float min = Float.POSITIVE_INFINITY, max = Float.NEGATIVE_INFINITY;
        for (Vec2 p : P) {
            float d = p.dot(axis);
            min = Math.min(min, d);
            max = Math.max(max, d);
        }
        return new float[]{min, max};
    }

    // --- Data structures

    public static final class Vertex {
        public final float x, y, z;
        public Vertex(float x, float y, float z) { this.x = x; this.y = y; this.z = z; }
        public Vec3 toVec3() { return new Vec3(x, y, z); }
        @Override
        public String toString() { return String.format( "(%f, %f, %f)", x, y, z ); }
    }

    private static final class Vec3 {
        final float x, y, z;
        Vec3(float x, float y, float z) { this.x = x; this.y = y; this.z = z; }
        //Vec3 add(Vec3 o) { return new Vec3(x + o.x, y + o.y, z + o.z); }
        Vec3 sub(Vec3 o) { return new Vec3(x - o.x, y - o.y, z - o.z); }
        Vec3 scale(float s) { return new Vec3(x * s, y * s, z * s); }
        //float dot(Vec3 o) { return x * o.x + y * o.y + z * o.z; }
        Vec3 cross(Vec3 o) {
            return new Vec3(
                y * o.z - z * o.y,
                z * o.x - x * o.z,
                x * o.y - y * o.x
            );
        }
        float length() { return (float)Math.sqrt(x * x + y * y + z * z); }
    }

    private static final class Vec2 {
        final float x, y;
        Vec2(float x, float y) { this.x = x; this.y = y; }
        Vec2 sub(Vec2 o) { return new Vec2(x - o.x, y - o.y); }
        float dot(Vec2 o) { return x * o.x + y * o.y; }
        Vec2 normalize(float eps) {
            float len = (float)Math.sqrt(x*x + y*y);
            if (len < eps) return new Vec2(0, 0);
            return new Vec2(x/len, y/len);
        }
    }

    private Vec2[] segmentHull2D(Vec2 p, Vec2 q, Vec2 r) {
        // Return minimal polygon (point or segment) covering nearly collinear triplet
        // Pick farthest pair
        float dPQ = dist2(p,q), dPR = dist2(p,r), dQR = dist2(q,r);
        if (dPQ >= dPR && dPQ >= dQR) return new Vec2[]{p, q};
        if (dPR >= dPQ && dPR >= dQR) return new Vec2[]{p, r};
        return new Vec2[]{q, r};
    }

    private float dist2(Vec2 a, Vec2 b) {
        float dx = a.x - b.x, dy = a.y - b.y;
        return dx*dx + dy*dy;
    }

    @Override
    public String toString() {  return "{ "+v[0].toString()+", "+v[1].toString()+", "+v[2].toString()+" }";  }
    
}
