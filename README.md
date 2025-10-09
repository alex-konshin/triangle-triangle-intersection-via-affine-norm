# Triangle–Triangle Intersection via Affine Normalization

Author: Alex Konshin  
First described: 2010  
License: CC0 1.0 (Public Domain)

---

## Motivation

We only need a yes/no intersection test. That frees us to apply any transform preserving incidence. Map one triangle to the unit right triangle on the plane z=0; transform the other with the same affine map. The 3D problem reduces to “does the second triangle cross z=0, and does that crossing segment intersect the unit triangle x ≥ 0, y ≥ 0, x + y ≤ 1?”

This approach avoids explicit plane–plane intersection lines, axis-selection heuristics, and most degenerate-case branching. It also easily yields the actual intersection geometry by transforming canonical results back.

---

## Core idea

- Fix triangle A = {a0, a1, a2}. Define:
  - e1 = a1 − a0
  - e2 = a2 − a0
  - n = e1 × e2
  - M = [ e1 e2 n ] (3×3 with these columns)
- The inverse M⁻¹ defines the world→canonical mapping:
  - p′ = M⁻¹ (p − a0)
  - In canonical space: A becomes {(0,0,0), (1,0,0), (0,1,0)} on z=0.

Transform triangle B with the same mapping; test it in this canonical frame.

<!-- TODO: Insert diagram of affine normalization here -->

---

## Algorithm

1. Build Minv and t0 for triangle A:
   - Minv = M⁻¹, t0 = −Minv · a0 (precompute and reuse when A is static).
2. Transform B’s vertices: b′i = Minv · bi + t0 = (xi, yi, zi).
3. Quick reject by signs of z:
   - If all zi > ε or all zi < −ε → no intersection.
   - If all |zi| ≤ ε → coplanar case: do a 2D triangle–triangle overlap in (x,y).
4. Non-coplanar:
   - Intersect B’s edges with z=0. For edge (i,j) where zi and zj have opposite signs:
     - t = zi / (zi − zj)
     - p = (1−t) b′i + t b′j (lies on z=0)
   - You get a segment S = [p0, p1] on z=0 (or a single point).
5. Test S against the unit triangle:
   - Clip S by the half-planes: x ≥ 0, y ≥ 0, x + y ≤ 1. If anything remains, triangles intersect.

Optional refinement: use s = x + y as a monotone parameter along S when z=0.

---

## Recovering intersection geometry

If intersecting, map canonical endpoints back: P = M · p + a0.  
For coplanar overlap, compute the 2D intersection polygon in (x,y) and map each vertex back similarly.

---

## Reference implementation

### Java implementation
See the full Java source code [here](/triangle-triangle-intersection-via-affine-norm/blob/main/rc/main/java/com/alexkonshin/Face.java).

### C implementation
See [`src/c/tri_tri_affine.c`](src/c/tri_tri_affine.c) for a C implementation.
> **Note:** The C implementation is provided for demonstration purposes only.  
> It illustrates the core idea of the algorithm and is not intended to cover all possible edge or degenerate cases.  
> It is also not optimized for high performance.

---

### Build & Run

Gradle build script is provided to build and test Java source code.

You can build and run the C test harness included in `src/c` to verify the algorithm.
The test program runs a few sample triangle–triangle intersection cases and prints YES or NO for each.

#### Using Make
```bash
cd src/c
make
./tri_test
```

#### Using CMake
```cd src/c
mkdir build && cd build
cmake ..
cmake --build .
./tri_test
```

---

## License

This work is dedicated to the public domain under CC0 1.0.
