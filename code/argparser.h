// This is a standard way to ensure that the code inside is included only once during compilation, preventing multiple definitions if the header file is included in multiple source files.
#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include "string.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"

// takes in command-line parsed arguments 
class ArgParser {

public:
  // contructor 1, this initializes the class with the default values
  ArgParser() { DefaultValues(); }
  // second constructor takes in command-line arguments and parses them to set the values of member variables accordingly
  ArgParser(int argc, char *argv[]) {
    DefaultValues();

    for (int i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "-input")) {
	      i++;
        // stelling in assert is true (1), doet gewoon verder met het programma
        assert (i < argc); 
	      input_file = argv[i];
      } else if (!strcmp(argv[i], "-size")) {
	      i++; 
        assert (i < argc); 
        // predefined function from the stdlib header file used to convert a string value to an integer value
	      width = height = atoi(argv[i]);
      } else if (!strcmp(argv[i], "-wireframe")) {
        wireframe = true;
      } else if (!strcmp(argv[i], "-gouraud")) {
        gouraud = true;
      } else {
        // if unrecognised command-line argument is encountered
	      printf ("whoops error with command line argument %d: '%s'\n",i,argv[i]);
        // terminates the program
	      assert(0);
      }
    }
  }

  void DefaultValues() {
    input_file = NULL;
    width = 600;
    height = 600;
    wireframe = 0;
  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)

  // pointer to character array (string) representing the input file name
  char *input_file;
  // width of window where object is shown
  int width;
  // height of window where object is shown
  int height;
  // boolean flag wheter or not to show edges
  bool wireframe;
  // boolean to set shading or not
  bool gouraud;

};

#endif
