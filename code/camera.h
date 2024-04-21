#ifndef _CAMERA_H_
#define _CAMERA_H_

#include <assert.h>
#include "vectors.h"

// ====================================================================

class Camera {

public:
  // CONSTRUCTOR & DESTRUCTOR
  // Gives center, direction and up as arguments
  Camera(Vec3f c, Vec3f d, Vec3f u);
 ~Camera() {}

  // GL NAVIGATION
  // Explenation and implementation of functions is to be found in camera.C
  virtual void glInit(int w, int h) = 0;
  void glPlaceCamera(void);
  virtual void dollyCamera(float dist) = 0;
  virtual void truckCamera(float dx, float dy) = 0;
  virtual void rotateCamera(float rx, float ry) = 0;
  virtual void Print() = 0;

protected:
  Camera() { assert(0); } // don't use

  // HELPER FUNCTIONS
  const Vec3f getHorizontal() const {
    Vec3f answer;
    // calculate cross product of direction and up and puts the answer in answer => see vectors.h
    // cross product => a method of multiplying two vectors that produces a vector perpendicular to both vectors involved in the multiplication
    // This function computes a vector representing the horizontal direction relative to the given direction and up vector.
    Vec3f::Cross3(answer, direction, up);
    answer.Normalize();
    return answer; 
  }

  // This function computes a vector representing the direction considered as "up" on the screen.
  const Vec3f getScreenUp() const {
    Vec3f answer;
    Vec3f::Cross3(answer, getHorizontal(), direction);
    return answer; 
  }

  // REPRESENTATION
  Vec3f center;
  Vec3f direction;
  Vec3f up;
};

// ====================================================================

class PerspectiveCamera : public Camera {

public:
  // CONSTRUCTOR & DESTRUCTOR
  // this extension of camera also gets an angle as argument
  PerspectiveCamera(Vec3f c, Vec3f d, Vec3f u, float a);
 ~PerspectiveCamera(void) { }

  // GL NAVIGATION
  void glInit(int w, int h);
  void dollyCamera(float dist);
  void truckCamera(float dx, float dy);
  void rotateCamera(float rx, float ry);
  void Print() {
    printf ("PerspectiveCamera {\n");
    printf ("    center    ");
    // see vectors.h
    center.Write(stdout);
    printf ("    direction ");
    direction.Write(stdout);
    printf ("    up        ");
    up.Write(stdout);
    printf ("    angle      %f\n", angle);
    printf ("}\n");
  }    

private:
  PerspectiveCamera() { assert(0); } // don't use

  // REPRESENTATION
  float angle;
  Vec3f lowerLeft;
  Vec3f xAxis;
  Vec3f yAxis;
};

// ====================================================================

#endif
