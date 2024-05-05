#ifndef _VERTEX_H
#define _VERTEX_H

#include <eigen3/Eigen/Dense>
#include <stdio.h>
#include <assert.h>

#include "vectors.h"

class Vertex;

// ==========================================================

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  // triangleCount is added by us
  Vertex(int i, const Vec3f &pos) : position(pos), triangleCount(0) { 
    index = i; 
  }

  virtual ~Vertex() { }
  
  // =========
  // ACCESSORS
  int getIndex() const { 
    return index; 
  }

  double x() const { 
    return position.x(); 
  }

  double y() const { 
    return position.y(); 
  }

  double z() const { 
    return position.z(); 
  }

  const Vec3f& get() const { 
    return position; 
  }

  // We added this
  const Eigen::Matrix4d getQem() const {
    return qem;
  }

  const Vec3f getNormal() const {
    return normal;
  }

  // added by us
  int getTriangleCount() const { 
    return triangleCount; 
  }

  // =========
  // MODIFIERS
  void set(Vec3f v) { 
    position = v; 
  }

  void set(double x, double y, double z) { 
    position.Set(x, y, z); 
  }

  // We added this
  void set(Eigen::Matrix4d q) {
    qem = q;
  }

  // We added this
  void clearQEM() {
    qem = Eigen::Matrix4d::Zero();
  }
  void setNormal(Vec3f n) {
    normal = n;
  }
  void addNormal(Vec3f n) {
    normal += n;
  }

  // added by us
  void incrementTriangleCount() { triangleCount++; }
  void decrementTriangleCount() { triangleCount--; }

private:

  // don't use these constructors
  Vertex() { assert(0); }
  Vertex& operator=(const Vertex&) { assert(0); }
  Vertex(const Vertex&) { assert(0); }
  
  // ==============
  // REPRESENTATION
  Vec3f position;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure
  int index;  

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation.

    // Count of triangles associated with this vertex
  int triangleCount;

  //normal of the vertex used in gouraud shading
  Vec3f normal = Vec3f(0,0,0);

  // We added this
  // used for saving QEM of vertex
  // Initialised very large 
  Eigen::Matrix4d qem = Eigen::Matrix4d::Zero();
};

// ==========================================================

#endif

