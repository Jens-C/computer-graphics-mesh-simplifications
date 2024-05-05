#include <eigen3/Eigen/Dense>
#include <stdio.h>
#include <assert.h>
#include <GL/gl.h>
#include <fstream>
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"
#include "vertex_parent.h"
#include "glCanvas.h"
// We included this
#include "iostream"
#include <cstdlib>
#include <list>
#include <chrono>


#define INITIAL_VERTEX 10000
#define INITIAL_EDGE 10000
#define INITIAL_TRIANGLE 10000



// =======================================================================
// CONSTRUCTORS & DESTRUCTORS
// =======================================================================

Mesh::Mesh() {
  vertices = new Array<Vertex*>(INITIAL_VERTEX);
  edges = new Bag<Edge*>(INITIAL_EDGE, Edge::extract_func);
  triangles = new Bag<Triangle*>(INITIAL_TRIANGLE, Triangle::extract_func);
  vertex_parents = new Bag<VertexParent*>(INITIAL_VERTEX, VertexParent::extract_func);
  bbox = NULL;
}

Mesh::~Mesh() {
  delete vertices;  
  vertices = NULL;
  delete edges;
  edges = NULL;
  delete triangles;
  triangles = NULL;
  delete bbox;
  bbox = NULL;
}

// =======================================================================
// MODIFIERS:   ADD & REMOVE
// =======================================================================

Vertex* Mesh::addVertex(const Vec3f &position) {
  int index = numVertices();
  Vertex *v = new Vertex(index, position);
  vertices->Add(v);
  if (bbox == NULL) 
    bbox = new BoundingBox(position, position);
  else 
    bbox->Extend(position);
  return v;
}

void Mesh::addTriangle(Vertex *a, Vertex *b, Vertex *c) {

  // create the triangle
  Triangle *t = new Triangle();

  // create the edges
  Edge *ea = new Edge(a, t);
  Edge *eb = new Edge(b, t);
  Edge *ec = new Edge(c, t);

  // point the triangle to one of its edges
  t->setEdge(ea);

  // connect the edges to each other
  ea->setNext(ec);
  eb->setNext(ea);
  ec->setNext(eb);

  // add them to the master list
  edges->Add(ea);
  edges->Add(eb);
  edges->Add(ec);


  // connect up with opposite edges (if they exist)
  // getEdge returns the edge that exists between two points, these two points are accessed by [0] and [1]
  Edge *ea_op = getEdge((*ea)[1], (*ea)[0]);
  Edge *eb_op = getEdge((*eb)[1], (*eb)[0]);
  Edge *ec_op = getEdge((*ec)[1], (*ec)[0]);  

  // The opposite edge is actually the same as the one that exists, but is points in the different direction
  if (ea_op != NULL) { ea_op->setOpposite(ea); }// else{printf("ea_op is NULL\n");}
  if (eb_op != NULL) { eb_op->setOpposite(eb); }// else{printf("eb_op is NULL\n");}
  if (ec_op != NULL) { ec_op->setOpposite(ec); }// else{printf("ec_op is NULL\n");}

  // add the triangle to the master list
  triangles->Add(t); 
}



void Mesh::removeTriangle(Triangle *t) {
  Edge *ea = t->getEdge();
  Edge *eb = ea->getNext();
  Edge *ec = eb->getNext();
  assert (ec->getNext() == ea);
  

  // remove elements from master lists
  // We added the if statements
  if (edges->Member(ea))
    edges->Remove(ea);
  if (edges->Member(eb))
    edges->Remove(eb);
  if (edges->Member(ec))
    edges->Remove(ec);


  
  // Check if the triangle is already removed before attempting to remove it
  // We added the if statements
  if (triangles->Member(t))
    triangles->Remove(t);

  delete ea;
  delete eb;
  delete ec;
  delete t;
}

Edge* Mesh::getEdge(Vertex *a, Vertex *b) const {
  assert (edges != NULL);
  return edges->Get(a->getIndex(), b->getIndex());
}

// return the vertices in order => from small to big
Vertex* Mesh::getChildVertex(Vertex *p1, Vertex *p2) const {
  VertexParent *vp = vertex_parents->GetReorder(p1->getIndex(), p2->getIndex());
  if (vp == NULL) 
    return NULL;
  return vp->get();
}

void Mesh::setParentsChild(Vertex *p1, Vertex *p2, Vertex *child) {
  vertex_parents->Add(new VertexParent(p1,p2,child));
}

Vec3f ComputeNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
  Vec3f v12 = p2;
  v12 -= p1;
  Vec3f v23 = p3;
  v23 -= p2;
  Vec3f normal;
  Vec3f::Cross3(normal,v12,v23);
  normal.Normalize();
  return normal;
}
//calculate normals of the faces, then add them to the vertices (normal of face is not normalized yet )
void Mesh::CalcNormals(){
  Iterator<Triangle*> *iter = triangles->StartIteration();
    while (Triangle *t = iter->GetNext()) {
      Vec3f normal = ComputeNormal((*t)[0]->get(), (*t)[1]->get(), (*t)[2]->get());
      normal.Normalize();
      (*t)[0]->addNormal(normal);
      (*t)[1]->addNormal(normal);
      (*t)[2]->addNormal(normal);
    }
    triangles->EndIteration(iter);

}

// =======================================================================
// PAINT
// =======================================================================



void InsertNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
  Vec3f normal = ComputeNormal(p1,p2,p3);
  glNormal3f(normal.x(), normal.y(), normal.z());
}

void Mesh::Paint(ArgParser *args) {

  // scale it so it fits in the window
  Vec3f center; 
  bbox->getCenter(center);
  float s = 1/bbox->maxDim();
  glScalef(s, s, s);
  glTranslatef(-center.x(), -center.y(), -center.z());

  // this offset prevents "z-fighting" bewteen the edges and faces
  // the edges will always win.
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.1, 4.0);

  // draw the triangles
  glBegin (GL_TRIANGLES); 
  if(args->gouraud){
    CalcNormals();
    glColor3f(1,1,1); //set collor to white 
    Iterator<Triangle*> *iter = triangles->StartIteration();
    while (Triangle *t = iter->GetNext()) {
      for (int i = 0; i < 3; ++i) {
        Vec3f vertex = (*t)[i]->get();
        Vec3f normal = (*t)[i]->getNormal();
        normal.Normalize();
        //add normal for every vertex
        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(vertex.x(), vertex.y(), vertex.z());
      }
    }
    triangles->EndIteration(iter);
  }else{
    glColor3f(1,1,1); //set collor to white
    Iterator<Triangle*> *iter = triangles->StartIteration();
    while (Triangle *t = iter->GetNext()) {
      Vec3f a = (*t)[0]->get();
      Vec3f b = (*t)[1]->get();
      Vec3f c = (*t)[2]->get();
      InsertNormal(a, b, c); 
      glVertex3f(a.x(), a.y(), a.z());
      glVertex3f(b.x(), b.y(), b.z());
      glVertex3f(c.x(), c.y(), c.z());
    }
    triangles->EndIteration(iter);
  }
  
  glEnd();

  glDisable(GL_POLYGON_OFFSET_FILL); 
 
  if (args->wireframe) {
    glDisable(GL_LIGHTING);

    // draw all the interior, non-crease edges
    glLineWidth(1);
    glColor3f(0, 0, 0);
    glBegin (GL_LINES);
    Iterator<Edge*> *iter = edges->StartIteration();
    while (Edge *e = iter->GetNext()) {
      if (e->getOpposite() == NULL || e->getCrease() > 0) 
        continue;
      Vec3f a = (*e)[0]->get();
      Vec3f b = (*e)[1]->get();
      glVertex3f(a.x(), a.y(), a.z());
      glVertex3f(b.x(), b.y(), b.z());
    }
    edges->EndIteration(iter);
    glEnd();

    // draw all the interior, crease edges
    glLineWidth(3);
    glColor3f(1, 1, 0);
    glBegin (GL_LINES);
    iter = edges->StartIteration();
    while (Edge *e = iter->GetNext()) {
      if (e->getOpposite() == NULL || e->getCrease() == 0) 
        continue;
      Vec3f a = (*e)[0]->get();
      Vec3f b = (*e)[1]->get();
      glVertex3f(a.x(), a.y(), a.z());
      glVertex3f(b.x(), b.y(), b.z());
    }
    edges->EndIteration(iter);
    glEnd();

    // draw all the boundary edges, this is the red line when an edge doesn't have an opposite edge
    glLineWidth(3);
    glColor3f(1, 0, 0);
    glBegin (GL_LINES);
    iter = edges->StartIteration();
    while (Edge *e = iter->GetNext()) {
      if (e->getOpposite() != NULL) 
        continue;
      assert (e->getCrease() == 0);
      Vec3f a = (*e)[0]->get();
      Vec3f b = (*e)[1]->get();
      glVertex3f(a.x(), a.y(), a.z());
      glVertex3f(b.x(), b.y(), b.z());
    }
    edges->EndIteration(iter);
    glEnd();

    glEnable(GL_LIGHTING);
  }
  
  HandleGLError(); 
}

// =======================================================================
// the load function parses very simple .obj files
// the basic format has been extended to allow the specification 
// of crease weights on the edges.
// =======================================================================
bool compareEdges(const Edge* a, const Edge* b) {
    return a->getQem() < b->getQem(); // Order by qem in ascending order
}

void Mesh::Load(const char *input_file) {
  
  // Open a text file for reading. The file must exist.
  FILE *objfile = fopen(input_file, "r");
  if (objfile == NULL) {
    printf ("ERROR! CANNOT OPEN '%s'\n",input_file);
    return;
  }

  char line[200];
  char token[100];
  char atoken[100];
  char btoken[100];
  char ctoken[100];
  char dtoken[100];
  char etoken[100];
  float x,y,z;
  int a,b,c,d,e;
  
  int index = 0;
  int vert_count = 0;
  int vert_index = 1;
  
  // The fgets() function reads characters from the current stream position up to and including the first new-line character (\n), up to the end of the stream, or until the number of characters read is equal to n-1, whichever comes first. The fgets() function stores the result in string and adds a null character (\0) to the end of the string. The string includes the new-line character, if read. If n is equal to 1, the string is empty.
  while (fgets(line, 200, objfile)) {   
    
    // Checks if the second-to-last character of the line array is a backslash
    if (line[strlen(line)-2] == '\\') {
      // Reads the next line from the file objfile and stores it in the token array, with a maximum length of 100 characters
      fgets(token, 100, objfile);	
      // Calculates the index of the character just before the backslash in the line array.
      int tmp = strlen(line) - 2;
      // Copies the contents of the token array into the line array, starting from the index tmp. This effectively concatenates the continuation line to the end of the original line.
      strncpy(&line[tmp], token, 100);
    }

    // reads formatted input from a string (line in this case) and stores it into corresponding variables. If successful, it will be 1; otherwise, it will be -1
    int token_count = sscanf (line, "%s\n",token);
    if (token_count == -1) 
      continue;
    a = b = c = d = e = -1;

    if (!strcmp(token, "usemtl") || !strcmp(token, "g")) {
      vert_index = 1; //vert_count + 1;
      index++;
    } else if (!strcmp(token, "v")) {
      vert_count++;
      sscanf (line, "%s %f %f %f\n", token, &x, &y, &z);
      addVertex(Vec3f(x,y,z));
    } else if (!strcmp(token,"f")) {
      int num = sscanf (line, "%s %s %s %s %s %s\n",token,
			atoken,btoken,ctoken,dtoken,etoken);
      sscanf (atoken, "%d", &a);
      sscanf (btoken, "%d", &b);
      sscanf (ctoken, "%d", &c);
      if (num > 4) 
        sscanf (dtoken, "%d", &d);
      if (num > 5) 
        sscanf (etoken, "%d", &e);
      a -= vert_index;
      b -= vert_index;
      c -= vert_index;
      if (d >= 0) 
        d -= vert_index;
      if (e >= 0) 
        e -= vert_index;
      assert (a >= 0 && a < numVertices());
      assert (b >= 0 && b < numVertices());
      assert (c >= 0 && c < numVertices());

      addTriangle(getVertex(a), getVertex(b), getVertex(c));

      if (d > -1) { 
        assert (d < numVertices()); 
        addTriangle(getVertex(a), getVertex(c), getVertex(d));
      }
      if (e > -1) { 
        assert (e < numVertices()); 
        addTriangle(getVertex(a), getVertex(d), getVertex(e));
      }
    } else if (!strcmp(token,"e")) {
      int num = sscanf (line, "%s %s %s %s\n", token, atoken, btoken, ctoken);
      assert (num == 4);
      sscanf (atoken, "%d", &a);
      sscanf (btoken, "%d", &b);
      if (!strcmp(ctoken,"inf")) 
        x = 1000000; // this is close to infinity...
      else sscanf (ctoken, "%f", &x);
      Vertex *va = getVertex(a);
      Vertex *vb = getVertex(b);
      Edge *ab = getEdge(va, vb);
      Edge *ba = getEdge(vb, va);
      assert (ab != NULL);
      assert (ba != NULL);
      ab->setCrease(x);
      ba->setCrease(x);
    } else if (!strcmp(token,"vt")) {
    } else if (!strcmp(token,"vn")) {
    } else if (token[0] == '#') {
    } else {
      printf ("LINE: '%s'",line);
    }
  }
  // We added this function

  //create the initial edgesSorted list
  Iterator<Edge*> *iter = edges->StartIteration();
  while(Edge* e = iter->GetNext()) {
    Eigen::Matrix4d qem1_matrix = assignQEM(e);
    Eigen::Matrix4d qem2_matrix = assignQEM(e->getOpposite());
    Vec3f newP = (e->getVertex()->get() + e->getOpposite()->getVertex()->get()) * 0.5;
    Eigen::Vector4d v = Eigen::Vector4d(newP.x(), newP.y(), newP.z(), 1.0);
    float qem = std::abs(v.transpose().dot(qem1_matrix * v)) + std::abs(v.transpose().dot(qem2_matrix * v));
    e->setQEM(qem);
    edgesSorted.push_back(e);
  }
  edges->EndIteration(iter);
  std::sort(edgesSorted.begin(), edgesSorted.end(), compareEdges);
}




// We added this function




void Mesh::calculateCostOfEdges(Array<Edge*> *EdgesToUpdate) {
  //recreating intire list quicker then removing edges being removed and adding the new ones.
  edgesSorted.clear();
  Iterator<Edge*> *iter = edges->StartIteration();
  while(Edge* e = iter->GetNext()) {
    //only update cost of edges that are changed
      if(EdgesToUpdate->Member(e)){
        Vec3f newP = (e->getVertex()->get() + e->getOpposite()->getVertex()->get()) * 0.5;
      // moved this to edge colapse function so we only do it for the vertices that changed
      //Eigen::Matrix4d qem1_matrix = assignQEM(e);
      //Eigen::Matrix4d qem2_matrix = assignQEM(e->getOpposite());
      Eigen::Vector4d v = Eigen::Vector4d(newP.x(), newP.y(), newP.z(), 1.0);
      float qem = std::abs(v.transpose().dot(e->getVertex()->getQem() * v)) + std::abs(v.transpose().dot(e->getOpposite()->getVertex()->getQem() * v));
      e->setQEM(qem);
    }
    edgesSorted.push_back(e);
  }
  edges->EndIteration(iter);

  //sorting list
  std::sort(edgesSorted.begin(), edgesSorted.end(), compareEdges);

}

// We added this function
Eigen::Matrix4d Mesh::assignQEM(Edge* e) {
  Vertex* vertex = e->getVertex();
  Edge* edge = e;
  Eigen::Matrix4d quadricMatrixTotal = Eigen::Matrix4d::Zero();
  do {
    Edge* next_edge = edge->getNext()->getOpposite();
    Vertex* v1 = edge->getVertex();
    Vertex* v2 = edge->getNext()->getVertex();
    Vertex* v3 = edge->getNext()->getNext()->getVertex();
    Vec3f normal1 = ComputeNormal(v1->get(), v2->get(), v3->get());
    float d = -normal1.Dot3(v1->get());
    Eigen::Matrix4d quadricMatrix = Eigen::Matrix4d::Zero();
    
    float a = 2.0;
    float b = 1.0;
    float c = 3.0;
    normal1.Get(a, b, c);

    quadricMatrix << a*a, a*b, a*c, a*d,
                    a*b, b*b, b*c, b*d,
                    a*c, b*c, c*c, d*c,
                    a*d, b*d, d*c, d*d;

    quadricMatrixTotal += quadricMatrix;
    edge = next_edge;
  } while(edge != e);

  vertex->set(quadricMatrixTotal);
  return quadricMatrixTotal;
}

// =================================================================
// save mesh to file
// =================================================================

void Mesh::exportToFile() {
  
  std::ofstream outFile("./data/out.obj");
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << "/data/out.obj" << " for writing." << std::endl;
        return;
    }
    int* vertexid = new int[numVertices()]; 
    // Write vertices
    for (int i = 0; i < numVertices(); ++i) {
        Vec3f vertex = getVertex(i)->get();
        vertexid[getVertex(i)->getIndex()] = i+1;
        outFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
    }
    // Write faces
      Iterator<Triangle*> *iter = triangles->StartIteration();
      while (Triangle *t = iter->GetNext()) {
        outFile << "f " << vertexid[(*t)[0]->getIndex()] << " " << vertexid[(*t)[1]->getIndex()]<< " " << vertexid[(*t)[2]->getIndex()] << std::endl;
      }
      
    triangles->EndIteration(iter);

    outFile.close();
    std::cout << "Mesh saved to " << "filename" << std::endl;
    delete[] vertexid;

}

// =================================================================
// SIMPLIFICATION
// =================================================================

// Collapse edge after selecting the edge with a point picker
void Mesh::collapseSelectedEdge(float x, float y, float z) {
  Iterator<Edge*> *iter = edges->StartIteration();
  Edge* shortestEdgeDistance = iter->GetNext();
  Vertex* v1 = (*shortestEdgeDistance)[0];
  Vertex* v2 = (*shortestEdgeDistance)[1];
  Vec3f newPos = (v1->get() + v2->get()) * 0.5;
  float shortestDistance = sqrt(pow((newPos.x()-x), 2) + pow((newPos.y()-y), 2) + pow((newPos.z()-z), 2));

  while(Edge* e = iter->GetNext()) {
    v1 = (*e)[0];
    v2 = (*e)[1];
    newPos = (v1->get() + v2->get()) * 0.5;
    float distance = sqrt(pow((newPos.x()-x), 2) + pow((newPos.y()-y), 2) + pow((newPos.z()-z), 2));
    if(distance < shortestDistance) {
      shortestEdgeDistance = e;
      shortestDistance = distance;
    }
  }
  edges->EndIteration(iter);
  bool gelukt = collapseEdge(shortestEdgeDistance);
  if(!gelukt) {
    printf("Onmogelijk deze edge te verwijderen");
  }
}


// Reduce the number of triangles to the desired count (target_tri_count)
void Mesh::Simplification(int target_tri_count) {
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    printf("Simplify the mesh! %d -> %d\n", numTriangles(), target_tri_count);

    // We added everything from here in this function
    // Perform simplification until the target number of triangles is reached
    while (numTriangles() > target_tri_count) {
        // Collapse the selected edge
        bool isCollapsed = false;
        int i = 0;
         //int s = edgesSorted.size();
        while(!isCollapsed) {
          //Edge* e = edgesSorted[s-i-1];
          Edge* e = edgesSorted[i];
          assert(edgesSorted[i]!= NULL);
          //Edge* e = edges->ChooseRandom();
          isCollapsed = collapseEdge(e);
          i++;
        }
    }

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    std::chrono::duration<double> duration = end - start;
    // Output the duration in seconds
    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
}


// We added this function to collapse an edge which can be selected by the quadratic error matric, random or with a point picker
bool Mesh::collapseEdge(Edge* edge) {
  
  if (edge == NULL) return false;

  // Get vertices and triangles involved in the edge collapse
  // Edge points to v1
  
  Vertex* v1 = (*edge)[0]; // this is the vertex associated with the edge
  Vertex* v2 = (*edge)[1];
  Vec3f newPos = (v1->get() + v2->get()) * 0.5;
  v2->set(newPos);

  //get the triangles that are connected to the edge
  Triangle* t1 = edge->getTriangle();
  Triangle* t2 = edge->getOpposite()->getTriangle();

  // Get the next and previous half-edges of the edge

  int end = 0;
  Edge* current = edge->getNext()->getOpposite();
  // store edge poiting to v2 that wont be changed
  Edge* pointTov2 = edge->getOpposite()->getNext()->getOpposite();

  //check for illegal edge collapses
  Array<Vertex*> *vertTemp = new Array<Vertex*>(INITIAL_VERTEX);
  vertTemp->Add(v1);
  vertTemp->Add(edge->getNext()->getVertex());

  //looping trough triangles connected to v1 and checking for double vertices
  while(current->getNext()->getOpposite()->getTriangle() != t2){
    if(vertTemp->Member(current->getNext()->getVertex())){
      return false;
    }
    vertTemp->Add(current->getNext()->getVertex());
    Edge* temp_next = current->getNext()->getOpposite();
    current = temp_next;
  }
  current = edge->getOpposite()->getNext()->getOpposite();
  vertTemp->Add(edge->getOpposite()->getNext()->getVertex());
  //looping trough triangles cottectd to v2 and checking for double vertices
  while(current->getNext()->getOpposite()->getTriangle() != t1){
    if(vertTemp->Member(current->getNext()->getVertex())){
      return false;
    }
    vertTemp->Add(current->getNext()->getVertex());
    Edge* temp_next = current->getNext()->getOpposite();
    current = temp_next;
  }
  //remove t1 after checking if the collapse is legal
  removeTriangle(t1);

  current = edge->getNext()->getOpposite();
  //looping trough triangles connected to v1 and deleting and readding them with v2
  while(end == 0) {
    if (current->getNext()->getOpposite()->getTriangle() == t2) {
      end = 1;
      removeTriangle(t2);
    }
    Edge* temp_next = current->getNext()->getOpposite();
    // current->getVertex() => this is v1 !!
    Vertex* b = current->getNext()->getVertex();
    Vertex* c = current->getNext()->getNext()->getVertex();
    Triangle* temp = current->getTriangle();
    removeTriangle(temp);
    addTriangle(v2, c, b);
    current = temp_next;
  }; 
  vertices->Remove(v1);   
  delete v1;

 //  update QEM for only for vertices and eges being changed
  Array<Edge*> *EdgesToUpdate = new Array<Edge*>(INITIAL_VERTEX);
  current = pointTov2;
  assignQEM(pointTov2);
  //loop around v2
  do{
    Edge* current2 = current->getNext();
    Edge* pointTovx = current2;
    //loop around every vertex surrounding v2
    do{
    EdgesToUpdate->AddNoDuplicates(current2);
    EdgesToUpdate->AddNoDuplicates(current2->getNext());
    EdgesToUpdate->AddNoDuplicates(current2->getNext()->getNext());
    //printf("f%d\n",current2->getNext()->getVertex()->getIndex());
      Edge* temp_next = current2->getNext()->getOpposite();
      current2 = temp_next;
    }while(current2 != pointTovx);
    EdgesToUpdate->AddNoDuplicates(current);
    EdgesToUpdate->AddNoDuplicates(current->getNext());
    EdgesToUpdate->AddNoDuplicates(current->getNext()->getNext());
    assignQEM(current->getNext());
    //printf("d%d\n",current->getNext()->getVertex()->getIndex());
    Edge* temp_next = current->getNext()->getOpposite();
    current = temp_next;
  }  while(current != pointTov2);
  //calc cust of the edges
  calculateCostOfEdges(EdgesToUpdate);
  //avoid memory leaks
  delete vertTemp;
  delete EdgesToUpdate;
  return true;
}


// =================================================================
