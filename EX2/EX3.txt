#include <iostream>
#include "blitz/array.h"

using namespace blitz;

int main() {
    GeneralArrayStorage<3> storage;

    TinyVector<int, 3> base(10, 0, 0);  
    storage.setBase(base);  

    TinyVector<int, 3> shape(5, 20, 20); 

    Array<int, 3> A(shape, storage); 

    A = 7;

    std::cout << "A(11,0,0) = " << A(12,0,0) << std::endl; 
    return 0;
}
