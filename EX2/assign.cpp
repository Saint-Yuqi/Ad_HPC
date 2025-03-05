// This uses features from C++17, so you may have to turn this on to compile
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include "blitz/array.h"
#include "tipsy.h"
using namespace blitz;

int main(int argc, char *argv[]) {
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100;
    if (argc>2) nGrid = atoi(argv[2]);

    std::ifstream io(argv[1],std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h),sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }
    tipsy::swap(h); // Don't forget to write this function in tipsy.h

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles" << std::endl;
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // Load the particles
    tipsy::dark d;
    for(int i=0; i<N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d),sizeof(d))) {
            std::cerr << "error reading particle" << std::endl;
            return errno;
        }
        tipsy::swap(d); // Don't forget to write this function in tipsy.h
        // Save the position
        r(i,0) = d.pos[0];
        r(i,1) = d.pos[1];
        r(i,2) = d.pos[2];

        // Save the mass
        m(i) = d.mass;
    }

    // Create Mass Assignment Grid
    Array<float,3> grid(nGrid,nGrid,nGrid);

    grid = 0;
    for(int pn=0; pn<N; ++pn) {
        float x = r(pn,0);
        float y = r(pn,1);
        float z = r(pn,2);

        // Convert x, y, z to grid indices (assuming normalized coordinates)
        int i = int((x + 1.0) / 2.0 * nGrid);
        int j = int((y + 1.0) / 2.0 * nGrid);
        int k = int((z + 1.0) / 2.0 * nGrid);

        // Ensure indices are within bounds
        i = std::max(0, std::min(nGrid-1, i));
        j = std::max(0, std::min(nGrid-1, j));
        k = std::max(0, std::min(nGrid-1, k));

        // Deposit the mass onto grid(i,j,k)
        grid(i,j,k) += m(pn);
    }
    // Compute Projected Density
    Array<float,2> projected(nGrid, nGrid);
    projected = 0;

    for(int i=0; i<nGrid; ++i) {
        for(int j=0; j<nGrid; ++j) {
            float max_density = 0;
            for(int k=0; k<nGrid; ++k) {
                max_density = std::max(max_density, grid(i,j,k));
            }
            projected(i,j) = max_density;
        }
    }
    std::cout << "Computing projected density..." << std::endl;

     // Save projected density to file
    std::ofstream outFile("projected_density.dat", std::ios::binary);
    if (!outFile) {
        std::cerr << "Error creating output file!" << std::endl;
        return 1;
    }
    outFile.write(reinterpret_cast<char*>(projected.data()), sizeof(float) * nGrid * nGrid);
    outFile.close();

    std::cout << "Projected density map saved as projected_density.dat" << std::endl;


    // Calculate projected density
    // - create a 2D map and initialize to zero
    // - loop over the 3D grid
    // - if the density of the grid projected onto the map is greater than the current value
    //   - update the current value


    // Write out the 2D map
    // Read into Python and plot


}
