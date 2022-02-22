#include "Cell.h"

// whenever an object of a derived class is created, a base class constructor is called
// constructor of cell: all members from CustomImage are initialized using the constructor from CustomImage
Cell::Cell(CustomImage& other) :CustomImage{ other } {};

