#include<iostream>
#include"blitz/array.h"
using namespace blitz;
int main () {
    Array < int ,3 > A (8 ,8 ,4);
    A = 7;
    for (int k = 0; k < 4; ++k) {
        std::cout << "Slice A(:,:, " << k << "):\n";
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                std::cout << A(i, j, k) << " ";
            }
            std::cout << std::endl;
        }
    }

    Array < int ,3 > B = A ( Range (5 ,7) , Range (5 ,7) , Range (0 ,2));
    B = 4;
    for (int k = 0; k < 4; ++k) {
        std::cout << "Slice A(:,:, " << k << "):\n";
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                std::cout << A(i, j, k) << " ";
            }
            std::cout << std::endl;
        }
    }
    

}