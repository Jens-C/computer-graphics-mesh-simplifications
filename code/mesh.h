#ifndef MESH_H
#define MESH_H

#include <eigen3/Eigen/Dense>
#include "vectors.h"
#include "array.h"
#include "bag.h"
#include "boundingbox.h"
#include "argparser.h"
#include <vector>

class Vertex;
class Edge;
class Triangle;
class VertexParent;

// ======================================================================
// ======================================================================

class Mesh {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Mesh();
  virtual ~Mesh();
  void Load(const char *input_file);
    
  // ========
  // VERTICES
  int numVertices() const { 
    return vertices->Count(); 
  }

  Vertex* addVertex(const Vec3f &pos);

  // this creates a relationship between 3 vertices (2 parents, 1 child)
  void setParentsChild(Vertex *p1, Vertex *p2, Vertex *child);

  // this accessor will find a child vertex (if it exists) when given
  // two parent vertices
  Vertex* getChildVertex(Vertex *p1, Vertex *p2) const;

  // look up vertex by index from original .obj file
  Vertex* getVertex(int i) const {
    assert (i >= 0 && i < numVertices());
    Vertex *v = (*vertices)[i];
    assert (v != NULL);
    return v; 
  }

  // We implemented this
  void calculateCostOfVerticesAndEdges();

  // We added this
  Eigen::Matrix4d assignQEM(Edge* edge);

  // =====
  // EDGES
  int numEdges() const { 
    return edges->Count(); 
  }

  // this efficiently looks for an edge with the given vertices, using a hash table
  Edge* getEdge(Vertex *a, Vertex *b) const;

  // =========
  // TRIANGLES
  int numTriangles() const { 
    return triangles->Count(); 
  }

  void addTriangle(Vertex *a, Vertex *b, Vertex *c);
  void removeTriangle(Triangle *t);

  // ===============
  // OTHER ACCESSORS
  BoundingBox* getBoundingBox() const { 
    return bbox; 
  }

  // ===============
  // OTHER FUNCTIONS
  void Paint(ArgParser *args);
  void LoopSubdivision();
  void Simplification(int target_tri_count);

  

  // We added this
  bool collapseEdge(Edge* triangle);
  // We added this
  // Edge* selectEdgeToCollapse();
  // We added this
  void removeUnusedVertices();

private:

  // ==============
  // REPRESENTATION
  // All points in the mesh in an array
  Array<Vertex*> *vertices;
  // All adges of the mesh in a hash table
  Bag<Edge*> *edges;
  // All triangles of the mesh in a hash table
  Bag<Triangle*> *triangles;
  BoundingBox *bbox;
  // All the vertex child-parent relations in a hash table
  Bag<VertexParent*> *vertex_parents;

  // We added this
  std::vector<Edge*> edgesSorted;
};

// ======================================================================
// ======================================================================


#endif




