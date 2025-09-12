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
