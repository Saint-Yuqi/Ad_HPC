#include <iostream>
#include "blitz/array.h"

int main() {

    blitz :: Array < float ,3 > data (10 ,8 ,6);

    //std::cout << data << std::endl;
    int counter = 0;
    for (int i = 0; i < data.extent(0); ++i) {
        for (int j = 0; j < data.extent(1); ++j) {
            for (int k = 0; k < data.extent(2); ++k) {
                data(i, j, k) = counter++;  // Assigning unique values
            }
        }
    }

    // Print the values manually
    for (int i = 0; i < data.extent(0); ++i) {
        std::cout << "Slice " << i << ":\n";
        for (int j = 0; j < data.extent(1); ++j) {
            for (int k = 0; k < data.extent(2); ++k) {
                std::cout << data(i, j, k) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

}