// This uses features from C++17, so you may have to turn this on to compile
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include "blitz/array.h"
#include "tipsy.h"
#include <chrono>
using namespace blitz;
using namespace std::chrono;

inline float W_cic(float r) {
    // r is the distance from the cell center in 1D
    float ar = fabs(r);
    if (ar < 1.0f) {
        return 1.0f - ar;
    }
    return 0.0f;
}

inline float W_tsc(float r) {
    float ar = fabs(r);
    if (ar < 0.5f) {
        // inner region
        return 0.75f - ar*ar;
    } else if (ar < 1.5f) {
        // outer region
        float tmp = 1.5f - ar;
        return 0.5f * tmp * tmp;
    }
    return 0.0f;
}

inline float W_pcs(float r) {
    float ar = fabs(r);
    if (ar < 1.0f) {
        return (4.0f - 6.0f * ar * ar + 3.0f * ar * ar * ar) / 6.0f;
    } else if (ar < 2.0f) {
        float tmp = 2.0f - ar;
        return (tmp * tmp * tmp) / 6.0f;
    } else {
        return 0.0f;
    }
}



int main(int argc, char *argv[]) {
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100;
    if (argc>2) nGrid = atoi(argv[2]);

    std::string scheme = "NGP";
    if (argc > 3)
        scheme = argv[3];

    auto start_total = high_resolution_clock::now(); 
    auto start_read = high_resolution_clock::now();
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
    auto end_read = high_resolution_clock::now();
    std::cout << "Reading file took "
              << duration<double>(end_read - start_read).count()
              << " s" << std::endl;

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles" << std::endl;
    Array<float,2> r(N,3);
    Array<float,1> m(N);
    auto start_mass_assignment = high_resolution_clock::now();
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

    if (scheme == "NGP") {
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
    }
    else if (scheme == "CIC"){
        for(int pn=0; pn<N; ++pn) {
            float X = (r(pn,0) + 1.0f) * 0.5f * nGrid; // Map from [-1,1] to [0,nGrid)
            float Y = (r(pn,1) + 1.0f) * 0.5f * nGrid;
            float Z = (r(pn,2) + 1.0f) * 0.5f * nGrid;
            
            // i0, j0, k0 = "floor" indices for the particle
            int i0 = (int) floor(X);
            int j0 = (int) floor(Y);
            int k0 = (int) floor(Z);
            
            // fractional distance from cell center
            float dx = X - i0;
            float dy = Y - j0;
            float dz = Z - k0;
    
            // For each neighbor in x:
            for (int di = 0; di <= 1; di++) {
                int ix = std::max(0, std::min(nGrid-1, i0+di));
                float wx = W_cic((X - ix)); // distance from the center of cell ix
    
                // Similarly for y
                for (int dj = 0; dj <= 1; dj++) {
                    int jy =  std::max(0, std::min(nGrid-1, j0+dj));
                    float wy = W_cic((Y - jy));
    
                    // Similarly for z
                    for (int dk = 0; dk <= 1; dk++) {
                        int kz = std::max(0, std::min(nGrid-1, k0+dk));
                        float wz = W_cic((Z - kz));
    
                        // Combined 3D weight
                        float w = wx * wy * wz;
    
                        // Deposit mass
                        grid(ix, jy, kz) += m(pn) * w;
                    }
                }
            }
        }
    }
    else if (scheme == "TSC"){
        for(int pn=0; pn<N; ++pn) {
            float X = (r(pn,0) + 1.0f) * 0.5f * nGrid; // Map from [-1,1] to [0,nGrid)
            float Y = (r(pn,1) + 1.0f) * 0.5f * nGrid;
            float Z = (r(pn,2) + 1.0f) * 0.5f * nGrid;
            
            // i0, j0, k0 = "floor" indices for the particle
            int i0 = (int) floor(X);
            int j0 = (int) floor(Y);
            int k0 = (int) floor(Z);
            
            // fractional distance from cell center
            float dx = X - i0;
            float dy = Y - j0;
            float dz = Z - k0;
    
            // For each neighbor in x:
            for (int di = 0; di <= 2; di++) {
                int ix = std::max(0, std::min(nGrid-1, i0+di));
                float wx = W_tsc((X - ix)); // distance from the center of cell ix
    
                // Similarly for y
                for (int dj = 0; dj <= 2; dj++) {
                    int jy =  std::max(0, std::min(nGrid-1, j0+dj));
                    float wy = W_tsc((Y - jy));
    
                    // Similarly for z
                    for (int dk = 0; dk <= 2; dk++) {
                        int kz = std::max(0, std::min(nGrid-1, k0+dk));
                        float wz = W_tsc((Z - kz));
    
                        // Combined 3D weight
                        float w = wx * wy * wz;
    
                        // Deposit mass
                        grid(ix, jy, kz) += m(pn) * w;
                    }
                }
            }
        }    
    }
    else if (scheme == "PCS"){
        for(int pn=0; pn<N; ++pn) {
            float X = (r(pn,0) + 1.0f) * 0.5f * nGrid; // Map from [-1,1] to [0,nGrid)
            float Y = (r(pn,1) + 1.0f) * 0.5f * nGrid;
            float Z = (r(pn,2) + 1.0f) * 0.5f * nGrid;
            
            // i0, j0, k0 = "floor" indices for the particle
            int i0 = (int) floor(X);
            int j0 = (int) floor(Y);
            int k0 = (int) floor(Z);
            
            // fractional distance from cell center
            float dx = X - i0;
            float dy = Y - j0;
            float dz = Z - k0;
    
            // For each neighbor in x:
            for (int di = 0; di <= 3; di++) {
                int ix = std::max(0, std::min(nGrid-1, i0+di));
                float wx = W_pcs((X - ix)); // distance from the center of cell ix
    
                // Similarly for y
                for (int dj = 0; dj <= 3; dj++) {
                    int jy =  std::max(0, std::min(nGrid-1, j0+dj));
                    float wy = W_pcs((Y - jy));
    
                    // Similarly for z
                    for (int dk = 0; dk <= 3; dk++) {
                        int kz = std::max(0, std::min(nGrid-1, k0+dk));
                        float wz = W_pcs((Z - kz));
    
                        // Combined 3D weight
                        float w = wx * wy * wz;
    
                        // Deposit mass
                        grid(ix, jy, kz) += m(pn) * w;
                    }
                }
            }
        }
    }
    else {
        std::cerr << "Unknown scheme: " << scheme << std::endl;
        return 1;
    }
    
    auto end_mass_assignment = high_resolution_clock::now();
    std::cout << "Mass assignment took"
              << duration<double>(end_mass_assignment - start_mass_assignment).count()
              << " s" << std::endl;

    // Compute Projected Density
    Array<float,2> projected(nGrid, nGrid);
    projected = 0;
    auto start_projection = high_resolution_clock::now();
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
    auto end_projection = high_resolution_clock::now();
    std::cout << "Projection took"
              << duration<double>(end_projection - start_projection).count()
              << " s" << std::endl;
     // Save projected density to file
     std::ostringstream outFileName;
     outFileName << "projected_density_"<< scheme << "_" << nGrid << ".dat";
     
     std::ofstream outFile(outFileName.str(), std::ios::binary);
     if (!outFile) {
         std::cerr << "Error creating output file!" << std::endl;
         return 1;
     }
     outFile.write(reinterpret_cast<char*>(projected.data()), sizeof(float) * nGrid * nGrid);
     outFile.close();
     
     std::cout << "Projected density map saved as " << outFileName.str() << std::endl;


    // Calculate projected density
    // - create a 2D map and initialize to zero
    // - loop over the 3D grid
    // - if the density of the grid projected onto the map is greater than the current value
    //   - update the current value


    // Write out the 2D map
    // Read into Python and plot


}
