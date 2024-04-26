#ifndef EDGE_H
#define EDGE_H

#include <limits.h>
#include <stdio.h>
#include <assert.h>

class Vertex;
class Triangle;

// ===================================================================
// half-edge data structure => pointer in both directions along the edge

class Edge { 

public:

  // ========================
  // CONSTRUCTORS & DESTRUCTOR
  Edge(Vertex *v, Triangle *t);
  ~Edge();

  // here the values are assigned that are used in the hash function so they can be efficiently accessed within the Bag data structure
  static void extract_func(Edge *e, int &a, int &b, int &c);
 
  // =========
  // ACCESSORS
  Vertex* getVertex() const { 
    assert (vertex != NULL); 
    return vertex; 
  }
  //function written by us to remove vertex and add new one.
  void setVertex(Vertex* v) {
    assert(v != nullptr);
    vertex = v;
  }

  Edge* getNext() const { 
    assert (next != NULL); 
    return next; 
  }

  Triangle* getTriangle() const { 
    assert (triangle != NULL); 
    return triangle; 
  }

  Edge* getOpposite() const {
    // warning!  the opposite edge might be NULL!
    return opposite; 
  }

  float getCrease() const { 
    return crease; 
  }

  Vertex* operator[](int i) const { 
    if (i==0) 
      return getVertex();
    if (i==1) 
      return getNext()->getNext()->getVertex();
    assert(0);
  }

  // =========
  // MODIFIERS
  void setOpposite(Edge *e) {
    assert (opposite == NULL); 
    assert (e != NULL);
    assert (e->opposite == NULL);
    opposite = e; 
    e->opposite = this; 
  }

  void clearOpposite() { 
    if (opposite == NULL) 
      return; 
    assert (opposite->opposite == this); 
    opposite->opposite = NULL;
    opposite = NULL; 
  }

  void setNext(Edge *e) {
    //assert (next == NULL);
    assert (e != NULL);
    assert (triangle == e->triangle);
    next = e;
  }
  void setCrease(float c) { 
    crease = c; 
  }

private:

  Edge(const Edge&) { assert(0); }
  Edge& operator=(const Edge&) { assert(0); }

  // ==============
  // REPRESENTATION
  // in the half edge data adjacency data structure, the edge stores everything!
  // Pointer to the vertex associated with the edge. A vertex has an x, y and z position
  Vertex *vertex;
  // Pointer to the triangle associated with the edge.
  // Triangle contains the 3 edges associated with it
  Triangle *triangle;
  // Pointer to the opposite edge (may be NULL).
  // It is the same edge, but in the other direction
  Edge *opposite;
  // Pointer to the next edge in the same triangle.
  Edge *next;
  // indicates wheter the edge is sharp or not
  float crease;
};

// ===================================================================

#endif
