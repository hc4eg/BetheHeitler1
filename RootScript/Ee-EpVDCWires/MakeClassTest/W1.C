#include "W.h"
#include <iostream>

class W1: public W{
public:
  void Plot();
  void Add(Int_t i , Int_t j) { cerr << i << j << endl;};
  ~W1() { return;};
};

void W1::Plot(){
  cerr << "Error occur." << endl;
}
