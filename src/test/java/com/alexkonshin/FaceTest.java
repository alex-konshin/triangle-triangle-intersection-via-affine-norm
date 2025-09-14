/*
 * Created on September 8, 2025
 * @author Alex Konshin <akonshin@gmail.com>
 */
package com.alexkonshin;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

import com.alexkonshin.Face.Vertex;


class FaceTest {

  @Test
  void test() {
    testIntersection("1. Tilted intersection",
        t( 1.0f, 2.0f, 3.0f,
           4.0f, 2.0f, 3.0f,
           2.5f, 4.0f, 3.0f
        ),
        t( 2.5f, 3.0f, 2.0f,
           2.5f, 3.0f, 4.0f,
           3.0f, 3.5f, 3.0f
        ),
        true
    );
    
    testIntersection("2. Parallel planes, no intersection",
        t( 0.0f, 0.0f, 0.0f,
           2.0f, 0.0f, 0.0f,
           0.0f, 2.0f, 0.0f
        ),
        t( 0.0f, 0.0f, 1.0f,
           2.0f, 0.0f, 1.0f,
           0.0f, 2.0f, 1.0f
        ),
        false
    );
    
    testIntersection("3. Touching at vertex",
        t( -1.0f, -1.0f, 0.0f,
           1.0f, -1.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.0f, 1.0f, -1.0f,
           0.0f, 1.0f, 1.0f,
           0.5f, 2.0f, 0.0f
        ),
        true
    );
    
    testIntersection("4. Embedded above, no intersection",
        t( 0.0f, 0.0f, 0.0f,
           2.0f, 0.0f, 0.0f,
           0.0f, 2.0f, 0.0f
        ),
        t( 0.5f, 0.5f, 1.0f,
           1.5f, 0.5f, 1.0f,
           0.5f, 1.5f, 1.0f
        ),
        false
    );
        
    testIntersection("5. Edge intersection",
        t( 1.0f, 1.0f, 0.0f,
           3.0f, 1.0f, 0.0f,
           1.0f, 3.0f, 0.0f
        ),
        t( 2.0f, 0.0f, -1.0f,
           2.0f, 0.0f, 1.0f,
           2.0f, 2.0f, 0.0f
        ),
        true
    );
    testIntersection("6.1 Canonical triangle test (no transform)",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.25f, 0.25f, -1.0f,
           0.25f, 0.25f, 1.0f,
           0.5f, 0.5f, 0.0f
        ),
        true
    );

    testIntersection("6.2 Canonical triangle test (no transform)",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.5f, 0.5f, -1.0f,
           0.5f, 0.5f, 1.0f,
           0.75f, 0.75f, 0.0f
        ),
        true
    );

  }
  
  @Test
  void test2() {
    System.err.println("test2 ======== Edge cases =========");
    
    // 1. Intersection along OY axis
    testIntersection("1.1 Touching along OY edge",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.0f, 0.0f, -1.0f,
           0.0f, 1.0f, -1.0f,
           0.0f, 0.5f, 1.0f
        ),
        true
    );
    
    testIntersection("1.2 Touching at vertex on OY",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.0f, 0.5f, 0.0f,
           0.5f, 0.5f, 0.0f,
           0.5f, 1.0f, 0.0f
        ),
        true
    );
    
    // 2. Intersection along diagonal x + y = 1
    testIntersection("2.1 Touching along diagonal x+y=1",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 1.0f, 0.0f, -1.0f,
           0.0f, 1.0f, -1.0f,
           1.2f, -0.2f, 0.0f
        ),
        true
    );
    
    testIntersection("2.2 Touching at vertex on diagonal x+y=1",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.5f, 0.5f, 0.0f,
           0.7f, 0.3f, 0.0f,
           0.3f, 0.7f, 0.0f
        ),
        true
    );
    
    // 3. Vertex of unit triangle inside the other triangle
    testIntersection("3.1 Vertex (0,0,0) inside other triangle",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( -0.1f, -0.1f, 0.0f,
           0.5f, -0.1f, 0.0f,
           -0.1f, 0.5f, 0.0f
        ),
        true
    );
    
    testIntersection("3.2 Vertex (1,0,0) inside other triangle",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.9f, -0.1f, 0.0f,
           1.5f, -0.1f, 0.0f,
           0.9f, 0.5f, 0.0f
        ),
        true
    );
    
    // 4. Coplanar cases
    testIntersection("4.1 Coplanar identical triangles",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        true
    );
    
    testIntersection("4.2 Coplanar partial overlap",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 0.5f, 0.0f, 0.0f,
           1.5f, 0.0f, 0.0f,
           0.5f, 1.0f, 0.0f
        ),
        true
    );
    
    testIntersection("4.3 Coplanar touching at edge",
        t( 0.0f, 0.0f, 0.0f,
           1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f
        ),
        t( 1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f,
           1.0f, 1.0f, 0.0f
        ),
        true
    );
    
  }
  
  
  private void testIntersection( String comment, Face tri1, Face tri2, boolean expected ) {
    boolean actual = tri1.intersects( tri2 );

    String message = comment+": tri1="+tri1.toString()+", tri2="+tri2.toString();
    if ( actual==expected ) {
      System.err.println( message );
    } else {
      fail( message );
    }
  }
  
  public static Face t(
    float ax, float ay, float az,
    float bx, float by, float bz,
    float cx, float cy, float cz
  ) { 
    return new Face(
      new Vertex(ax, ay, az),      
      new Vertex(bx, by, bz),      
      new Vertex(cx, cy, cz)
    );
  }
  
}
