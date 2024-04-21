#include "boundingbox.h"
#include <GL/gl.h>
#include <GL/glu.h>

// ====================================================================
// ====================================================================

// rendering a wireframe representation of a bounding box in a 3D scene using OpenGL
// This function is a member of the BoundingBox class and thus not modify the state of the a-object where the method is called on
void BoundingBox::paint() const {
  // draw a wireframe box to represent the boundingbox
  
  // sets current draw color to white
  glColor3f(1,1,1);
  // set linewidth to one pixel
  glLineWidth(1);
  // Disables lighting calculations, ensuring that the box is drawn without being affected by lighting.
  glDisable(GL_LIGHTING);
  // Begins the definition of a series of lines to be drawn.
  glBegin(GL_LINES);

  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(min.x(),max.y(),min.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(min.x(),max.y(),min.z());

  glVertex3f(min.x(),min.y(),min.z());
  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(min.x(),max.y(),min.z());
  glVertex3f(min.x(),max.y(),max.z());
  glVertex3f(max.x(),min.y(),min.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(max.x(),max.y(),min.z());
  glVertex3f(max.x(),max.y(),max.z());

  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(min.x(),min.y(),max.z());
  glVertex3f(min.x(),max.y(),max.z());
  glVertex3f(max.x(),max.y(),max.z());
  glVertex3f(max.x(),min.y(),max.z());
  glVertex3f(max.x(),max.y(),max.z());
  glVertex3f(min.x(),max.y(),max.z());

  // Ends the definition of the lines.
  glEnd();
  // Re-enables lighting calculations, restoring the default OpenGL state.
  glEnable(GL_LIGHTING);	   
}

// ====================================================================
// ====================================================================
