#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "blitz/array.h"
#include "tipsy.h"
#include <chrono>
#include <algorithm>
#include <omp.h>    // Include OpenMP header

using namespace blitz;
using namespace std::chrono;

// Example weight functions for higher-order schemes:
inline float W_cic(float r) {
    float ar = std::fabs(r);
    return (ar < 1.0f) ? 1.0f - ar : 0.0f;
}

inline float W_tsc(float r) {
    float ar = std::fabs(r);
    if (ar < 0.5f) {
        return 0.75f - ar * ar;
    } else if (ar < 1.5f) {
        float tmp = 1.5f - ar;
        return 0.5f * tmp * tmp;
    } else {
        return 0.0f;
    }
}

inline float W_pcs(float r) {
    float ar = std::fabs(r);
    if (ar < 1.0f) {
        return (4.0f - 6.0f * ar * ar + 3.0f * ar * ar * ar) / 6.0f;
    } else if (ar < 2.0f) {
        float tmp = 2.0f - ar;
        return (tmp * tmp * tmp) / 6.0f;
    } else {
        return 0.0f;
    }
}

// Helper for clamping indices
inline int clampIndex(int idx, int nGrid) {
    return std::max(0, std::min(nGrid - 1, idx));
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size] [scheme]" << std::endl;
        std::cerr << "  scheme options: NGP, CIC, TSC, PCS (default: NGP)" << std::endl;
        return 1;
    }

    // Default grid size and mass assignment scheme
    int nGrid = 100;
    if (argc > 2)
        nGrid = std::atoi(argv[2]);

    std::string scheme = "NGP";
    if (argc > 3)
        scheme = argv[3];

    // Timing: Start reading file
    auto start_total = high_resolution_clock::now();
    auto start_read = high_resolution_clock::now();

    std::ifstream io(argv[1], std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    // Read header and swap if necessary
    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h), sizeof(h))) {
        std::cerr << "Error reading header" << std::endl;
        return errno;
    }
    tipsy::swap(h);
    auto end_read = high_resolution_clock::now();
    std::cout << "Reading file took "
              << duration<double>(end_read - start_read).count()
              << " s" << std::endl;

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles" << std::endl;
    Array<float,2> r(N, 3);
    Array<float,1> m(N);
    auto start_mass_assignment = high_resolution_clock::now();

    tipsy::dark d;
    for (std::uint64_t i = 0; i < N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
            std::cerr << "Error reading particle" << std::endl;
            return errno;
        }
        tipsy::swap(d);
        r(static_cast<int>(i), 0) = d.pos[0];
        r(static_cast<int>(i), 1) = d.pos[1];
        r(static_cast<int>(i), 2) = d.pos[2];
        m(static_cast<int>(i)) = d.mass;
    }

    // Create Mass Assignment Grid and initialize to zero
    Array<float,3> grid(nGrid, nGrid, nGrid);
    grid = 0;

    // Parallelized Mass Assignment Loop using OpenMP
    #pragma omp parallel for schedule(static)
    for (std::uint64_t pn = 0; pn < N; ++pn) {
        // Map coordinates from [-1,1] to [0, nGrid)
        float X = (r(static_cast<int>(pn), 0) + 1.0f) * 0.5f * nGrid;
        float Y = (r(static_cast<int>(pn), 1) + 1.0f) * 0.5f * nGrid;
        float Z = (r(static_cast<int>(pn), 2) + 1.0f) * 0.5f * nGrid;

        if (scheme == "NGP") {
            int i = clampIndex(static_cast<int>(X), nGrid);
            int j = clampIndex(static_cast<int>(Y), nGrid);
            int k = clampIndex(static_cast<int>(Z), nGrid);
            #pragma omp atomic
            grid(i, j, k) += m(static_cast<int>(pn));
        }
        else if (scheme == "CIC") {
            int i0 = static_cast<int>(std::floor(X));
            int j0 = static_cast<int>(std::floor(Y));
            int k0 = static_cast<int>(std::floor(Z));
            for (int di = 0; di <= 1; di++) {
                int ix = clampIndex(i0 + di, nGrid);
                float wx = W_cic(X - (i0 + di));
                for (int dj = 0; dj <= 1; dj++) {
                    int jy = clampIndex(j0 + dj, nGrid);
                    float wy = W_cic(Y - (j0 + dj));
                    for (int dk = 0; dk <= 1; dk++) {
                        int kz = clampIndex(k0 + dk, nGrid);
                        float wz = W_cic(Z - (k0 + dk));
                        #pragma omp atomic
                        grid(ix, jy, kz) += m(pn) * wx * wy * wz;
                    }
                }
            }
        }
        else if (scheme == "TSC") {
            int i0 = static_cast<int>(std::floor(X));
            int j0 = static_cast<int>(std::floor(Y));
            int k0 = static_cast<int>(std::floor(Z));
            for (int di = -1; di <= 1; di++) {
                int ix = clampIndex(i0 + di, nGrid);
                float wx = W_tsc(X - (i0 + di));
                for (int dj = -1; dj <= 1; dj++) {
                    int jy = clampIndex(j0 + dj, nGrid);
                    float wy = W_tsc(Y - (j0 + dj));
                    for (int dk = -1; dk <= 1; dk++) {
                        int kz = clampIndex(k0 + dk, nGrid);
                        float wz = W_tsc(Z - (k0 + dk));
                        #pragma omp atomic
                        grid(ix, jy, kz) += m(pn) * wx * wy * wz;
                    }
                }
            }
        }
        else if (scheme == "PCS") {
            int i0 = static_cast<int>(std::floor(X));
            int j0 = static_cast<int>(std::floor(Y));
            int k0 = static_cast<int>(std::floor(Z));
            for (int di = -2; di <= 2; di++) {
                int ix = clampIndex(i0 + di, nGrid);
                float wx = W_pcs(X - (i0 + di));
                for (int dj = -2; dj <= 2; dj++) {
                    int jy = clampIndex(j0 + dj, nGrid);
                    float wy = W_pcs(Y - (j0 + dj));
                    for (int dk = -2; dk <= 2; dk++) {
                        int kz = clampIndex(k0 + dk, nGrid);
                        float wz = W_pcs(Z - (k0 + dk));
                        #pragma omp atomic
                        grid(ix, jy, kz) += m(pn) * wx * wy * wz;
                    }
                }
            }
        }
        else {
            #pragma omp critical
            {
                std::cerr << "Unknown scheme: " << scheme << std::endl;
            }
            exit(1);
        }
    }
    auto end_mass_assignment = high_resolution_clock::now();
    std::cout << "Mass assignment took "
              << duration<double>(end_mass_assignment - start_mass_assignment).count()
              << " s" << std::endl;

    // Compute Projected Density (example: maximum along the z-axis)
    Array<float,2> projected(nGrid, nGrid);
    projected = 0;
    auto start_projection = high_resolution_clock::now();
    for (int i = 0; i < nGrid; ++i) {
        for (int j = 0; j < nGrid; ++j) {
            float max_density = 0;
            for (int k = 0; k < nGrid; ++k) {
                max_density = std::max(max_density, grid(i, j, k));
            }
            projected(i, j) = max_density;
        }
    }
    auto end_projection = high_resolution_clock::now();
    std::cout << "Projection took "
              << duration<double>(end_projection - start_projection).count()
              << " s" << std::endl;

    // Save projected density to file (filename includes grid size)
    std::ostringstream outFileName;
    outFileName << "projected_density_" << nGrid << ".dat";
    std::ofstream outFile(outFileName.str(), std::ios::binary);
    if (!outFile) {
        std::cerr << "Error creating output file!" << std::endl;
        return 1;
    }
    outFile.write(reinterpret_cast<char*>(projected.data()), sizeof(float) * nGrid * nGrid);
    outFile.close();
    std::cout << "Projected density map saved as " << outFileName.str() << std::endl;

    return 0;
}
