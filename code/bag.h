#ifndef BAG_H
#define BAG_H
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#include "math.h"
#include <random>

template <class ITERATOR_ELEMENT> class Iterator;

// enumeration of statusses of elements in the bag
enum BAG_ELEMENT_MARK { BAG_MARK_NULL, BAG_MARK_DELETE, BAG_MARK_PRESENT };

#define MAX_ITERATORS 1

// NOTE:  These primes are used for hashing, and should
// not be the same as any of the primes used as sizes of 
// the hash tables (see utils.C)
#define LARGE_PRIME_A 10007
#define LARGE_PRIME_B 11003
#define LARGE_PRIME_C 12007

// Gives the closest prime number to x that is bigger than x
int NextLargestPrime(unsigned int x);  // defined in utils.C

// ======================================================================

// A bag is implemented with a hash table to allow efficient access and removal
// Hash tables are a type of data structure in which the address/ index value of the data element is generated from a hash function
template <class BAG_ELEMENT>
class Bag {

public:
  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  // takes an initial size and function pointer as argument
  Bag(int s, void (*e_func)(BAG_ELEMENT, int &a, int &b, int &c)) {
    // here's the hash function to use for edges or vertices so they
    // can be efficiently accessed within the Bag data structure
    extract_func = e_func;
    num_iterators = 0;
    // choosing to set your hash table length to a large prime number will greatly reduce the occurrence of collisions.
    size = NextLargestPrime(s);
    // This line dynamically allocates memory for an array of enum BAG_ELEMENT_MARK elements, which are used to mark the status of elements in the bag (whether they are present, deleted, or null). This array is used internally by the Bag class to manage the elements in the bag.
    marks = new enum BAG_ELEMENT_MARK[size];
    // dynamically allocates memory for an array of type BAG_ELEMENT, the result is assigned to the pointer data
    // When you see BAG_ELEMENT in the context of this code, it's indicating that the Bag class is templated, allowing you to create a bag of elements of any type you choose. In this way, BAG_ELEMENT serves as a placeholder for the type of elements the bag will contain, and you provide the actual type when you instantiate the Bag class
    data = new BAG_ELEMENT[size];
    for (int i = 0; i < size; i++)
      marks[i] = BAG_MARK_NULL;
    // counts the number of elements in the bag
    count = 0;
    // counts the number of elements deleted in the bag
    del_count = 0; 
  }

  // Destructor to deallocate memory using delete ...to avoid memory leaks
  virtual ~Bag() {
    assert (num_iterators == 0);
    delete [] data;
    delete [] marks; 
  }

  // =========
  // ACCESSORS 
  // returns the number of elements in the bag
  int Count() const { 
    return count; 
  }

  // Checks if an element is a member of the bag
  int Member(const BAG_ELEMENT e) const {
    int a, b, c;
    // This condition checks if e is not equal to zero.
    assert (e != (BAG_ELEMENT)0);
    // assigns values to a, b and c => see edge.C or vertex.C
    extract_func(e, a, b, c);
    // hash is gemaakt op basis van de index
    int orig = hash(a, b, c);
    int x = orig;
    while (1) {
      if (marks[x] == BAG_MARK_NULL)
	      return 0;
      if (marks[x] == BAG_MARK_PRESENT && data[x] == e)
	      return 1;
      // What does this function do ?
      x = skip(orig,x);
    }
  }

  // Choose a random element from the bag
  BAG_ELEMENT ChooseRandom() const {
    assert (Count() > 0);
    while (1) {    
      // generate random integer with range 0 to size
      // int random_int = int(floor(drand48() * size));
      // Changed code to generate random edge because it was always the same with this existing code
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> distrib(0, size - 1);
      // Generate a random integer within the range
      int random_int = distrib(gen);
      if (marks[random_int] == BAG_MARK_PRESENT)
        return data[random_int];
    }
  }

  BAG_ELEMENT GetReorder(int a, int b, int c) const {
    assert (a != b && b != c && a != c);
    if (a < b && a < c) 
      return Get(a,b,c);
    else if (b < a && b < c)
      return Get(b,c,a);
    assert (c < a && c < b);
    return Get(c,a,b); 
  }

  BAG_ELEMENT GetReorder(int a, int b) const {
    assert (a != b);
    if (a < b) 
      return Get(a,b);
    return Get(b,a); 
  }

  // Returns an element out of the data array
  BAG_ELEMENT Get(int a, int b) const {
    assert (a != b);
    int orig = hash(a, b, 0);
    int x = orig;
    while (1) {
      assert (x >= 0 && x < size);
      if (marks[x] == BAG_MARK_NULL)
	      return 0;
      if (marks[x] == BAG_MARK_PRESENT) {
	      int _a, _b, _c;
	      assert (data[x] != (BAG_ELEMENT)0);
	      extract_func(data[x], _a, _b, _c);
        if (_a == a && _b == b) {
          return data[x];
        }
      }
      x = skip(orig,x);
    }
  }

  // Start iterations over the elements of the bag
  Iterator<BAG_ELEMENT>* StartIteration() {
    //printf ("start iteration %d   %d\n", num_iterators, (int)this);

    // It prevents more than one iterator being active at the same time when MAX_ITERATIONS is set to 1 => see begin file
    assert (num_iterators < MAX_ITERATORS);
    num_iterators++;
    // Creates a new Iterator with type BAG_ELEMENT, which is a placeholder. The pointer "this" refers to the current object instance of bag
    return new Iterator<BAG_ELEMENT>(this);
  }

  void EndIteration(Iterator<BAG_ELEMENT> *&iter) {
    assert (num_iterators > 0);
    num_iterators--;
    assert (iter != NULL);
    delete iter;
    iter = NULL;    
    //printf ("end iteration %d   %d\n", num_iterators, (int)this);
  }

  // =========
  // MODIFIERS
  // Adds an element to the bag
  void Add(const BAG_ELEMENT e) {
    assert(!Member(e));
    // It first checks if the load factor exceeds 50%, and if so, it resizes the hash table, it does this by checking how many elements are in the hash table by counting the deleted and non deleted ones because the deleted ones stay in the array but get the enum BAG_MARK_DELETE
    if (count + del_count > size / 2)
      Resize(max2(count * 4, size));
    int a, b, c;
    assert (e != (BAG_ELEMENT)0);
    extract_func(e, a, b, c);
    // Computes the hash value of the element and probes for an empty slot using linear probing
    int orig = hash(a,b,c);
    int x = orig;
    while (marks[x] == BAG_MARK_PRESENT)
      x = skip(orig,x);
    if (marks[x] == BAG_MARK_DELETE)
      del_count--;
    marks[x] = BAG_MARK_PRESENT;
    data[x] = e;
    count++;
  }

  void AddNoDuplicates(const BAG_ELEMENT e) { 
    if (!Member(e)) 
      Add(e); 
  }

  // method removes an element from the bag. It computes the hash value of the element and probes for the slot containing the element using linear probing.
  void Remove(const BAG_ELEMENT e) {
    int a, b, c;
    assert (e != (BAG_ELEMENT)0);
    extract_func(e, a, b, c);
    int orig = hash(a, b, c);
    int x = orig;
    while (1) {
      assert (marks[x] != BAG_MARK_NULL);
      if (marks[x] == BAG_MARK_PRESENT && data[x] == e) {
	      marks[x] = BAG_MARK_DELETE;
	      del_count++;
	      count--;
	      break;
      }
      x = skip(orig,x);
    }
  }

  void DeleteAllElements() {
    assert(num_iterators == 0);
    for (int i = 0; i < size; i++) {
      if (marks[i] == BAG_MARK_PRESENT)
	      delete data[i];
      // Should this better be BAG_MARK_DELETE instead of BAG_MARK_NULL
      marks[i] = BAG_MARK_NULL;
    }
    del_count = 0;
    count = 0;
  }
 
  void Clear() {
    assert(num_iterators == 0);
    for (int i = 0; i < size; i++) {
      marks[i] = BAG_MARK_NULL;
    } 
    del_count = 0;
    count = 0;
  }
  
  void Print() {
    printf ("BAG::PRINT %d %d %d\n", size, count, del_count);
    int c=0;
    for (int i = 0; i < size; i++) {      
      //printf ("%3d\n: ",i);
      if (marks[i] == BAG_MARK_PRESENT) {
	    data[i]->Print();
	c++;
      } else if (marks[i] == BAG_MARK_DELETE) {
	//printf ("XXXXXXXXXX\n");
      } else {
	assert (marks[i] == BAG_MARK_NULL);
	//printf ("NULL\n");
      }
    }
    assert(c==count);
   
  }

  private: void Resize(int s) {
    assert (s > 0 && s > count);
    // save old stuff
    int old_size = size;
    BAG_ELEMENT *old_data = data;
    enum BAG_ELEMENT_MARK *old_marks = marks;
    // make new space
    // returns the next largest prime number to use as the size of the new data and marks array
    size = NextLargestPrime(s);
    marks = new enum BAG_ELEMENT_MARK[size];
    data = new BAG_ELEMENT[size];
    count = 0;
    del_count = 0;
    for (int i = 0; i < size; i++) {
      marks[i] = BAG_MARK_NULL;
      // why does every element gets the value 37 ?
      data[i] = (BAG_ELEMENT)37;
    }    
    int tmp_count = 0;
    // copy the stuff
    if (old_data != NULL) {
      for (int i = 0; i < old_size; i++) {
	      if (old_marks[i] == BAG_MARK_PRESENT) {
	        tmp_count++;
	        Add(old_data[i]);
	      }
      }
    }
    // del_count and count are still 0 at this point, is this ok ?
    // cleanup
    delete [] old_data;
    delete [] old_marks;
  }

  public: static void extract_int(int e, int &a, int &b, int &c) {
    a = e; 
    b = c = 0; 
  }
  
  // Create hash value
  protected: int hash (int a, int b, int c) const {
    int _a = (a < 0) ? 1 - 2 * a : 2 * a;
    int _b = (b < 0) ? 1 - 2 * b : 2 * b;
    int _c = (c < 0) ? 1 - 2 * c : 2 * c;
    int tmp = LARGE_PRIME_A * _a + LARGE_PRIME_B * _b + LARGE_PRIME_C * _c;
    // note: abs of the largest negative number is undefined...
    tmp = tmp % size;
    tmp = (tmp < 0) ? -tmp : tmp;
    tmp = tmp % size;
    assert (tmp >= 0);
    assert (tmp < size);
    return tmp;
  }
  
  // ==============
  // REPRESENTATION
  // size of the hash table
  int size;
  // counts the number of elements in the bag
  int count;
  // number of elements deleted
  int del_count;
  // Contains all the elements
  BAG_ELEMENT *data;
  // contains the status of the elements, wheter there is nothing, the element is present or deleted
  enum BAG_ELEMENT_MARK *marks;
  // Number of iterators => must be 0 when calling the function StartIteration
  int num_iterators;

  // extract function (an argument to the constructor)
  // here's the hash function to use for edges or vertices so they
  // can be efficiently accessed within the Bag data structure
  void (*extract_func)(BAG_ELEMENT, int &a, int &b, int &c);

  // skip function
  // Not sure for what this is used, why do we increment a hashed value
  inline int skip(int orig, int current) const { 
    assert (current >= 0 && current < size);
    int tmp = (current + 1) % size; 
    assert (current >= 0 && current < size);
    return tmp;
  }

  friend class Iterator<BAG_ELEMENT>;
  
};



// ======================================================================

template <class ITERATOR_ELEMENT>
class Iterator {

protected:

  // CONSTRUCTOR & DESTRUCTOR
  Iterator(Bag<ITERATOR_ELEMENT> *b) {
    bag = b;
    i = 0; 
  }

  virtual ~Iterator() {}

public:

  // ACCESSOR
  ITERATOR_ELEMENT GetNext() {
    ITERATOR_ELEMENT answer = (ITERATOR_ELEMENT)0;
    while (i < bag->size) {
      if (bag->marks[i] == BAG_MARK_PRESENT) {
        answer = bag->data[i];
        assert (answer != (ITERATOR_ELEMENT)0);
        break;
      }
      i++;      
    }
    i++;
    return answer;
  }
  
protected:
  
  Iterator() {}
  friend class Bag<ITERATOR_ELEMENT>;
  
  // ==============
  // REPRESENTATION
  // index
  int i;
  // the bag where the iteration needs to be performed
  Bag<ITERATOR_ELEMENT> *bag;

};




#endif
