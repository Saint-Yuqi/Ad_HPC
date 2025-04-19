// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <new>
#include <complex>
#include <stdlib.h>
#include "mpi.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "blitz/array.h"
#include "fftw3-mpi.h"
#include "tipsy.h"
#include "aweights.h"
#include <algorithm>
using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// A separate version is created for each different "Order".
// This allows the compiler to optimize the process for each of the four orders
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R,Array<float,1> &M) {
    auto nGrid = grid.rows(); // total number of particles
    // C++ Lambda to apply the periodic wrap of the grid index
    auto wrap = [nGrid](int i) {
        if (i<0) i+=nGrid;
        else if (i>=nGrid) i-=nGrid;
        return i;
    };
    #pragma omp parallel for
    for(int pn=R.lbound(0); pn<=R.ubound(0); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),Hy((y+0.5f)*nGrid),Hz((z+0.5f)*nGrid);
        for(auto i=0; i<Order; ++i) {
            for(auto j=0; j<Order; ++j) {
                for(auto k=0; k<Order; ++k) {
                    #pragma omp atomic
                    grid(wrap(Hx.i+i),wrap(Hy.i+j),wrap(Hz.i+k)) += m * Hx.H[i]*Hy.H[j]*Hz.H[k];
                }
            }
        }
    }
}
//Exercise 2
struct FourFloats {
    float x, y, z, m;
};
//Exercise 2

int main(int argc, char *argv[]) {
    int thread_support;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&thread_support);
    if (thread_support < MPI_THREAD_FUNNELED) {
        cerr << "Insufficient MPI thread support -- Funneled required" << std::endl;
        return 1;
    }
    int irank, nrank; // Index of current MPI rank, and total number of MPI ranks.
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);

    std::locale::global(std::locale("")); // e.g., LC_ALL=en_GB.UTF-8 or de_CH.UTF-8
    std::cerr.imbue(std::locale()); // e.g., LC_ALL=en_GB.UTF-8
    if (argc<=2) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std grid-size[/[L]bin-count] [order]"
                  << std::endl;
        return 1;
    }

#ifdef _OPENMP
    fftwf_init_threads();
#endif

    // The second parameter will have the grid size, and optionally the number of bins to use, e.g.,
    // 100      Use a grid size of 100, and use 50 bins
    // 100/40   Grid size is 100 and use 40 bins
    // 100/L35  Grid siye is 100, and use 35 logarithmic bins
    int nGrid = atoi(argv[2]);
    const auto iNyquist = nGrid/2;
    int nBins = iNyquist;
    bool bLog = false;
    auto p = strchr(argv[2],'/');
    if (p) {
        *p++ = '\0';
        if (*p == 'L') {
            bLog = true;
            ++p;
        }
        nBins = atoi(p);
    }

    int iOrder = 1;
    if (argc>3) iOrder = atoi(argv[3]);

    // Read the position and masses from the file
    auto t0 = hrc::now();
    std::ifstream io(argv[1],std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }
    constexpr std::streamsize buffer_size = 8 * 1024 * 1024;
    char* buffer = new char[buffer_size];
    io.rdbuf()->pubsetbuf(buffer, buffer_size);

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h),sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    if (irank==0) std::cerr << "Loading " << N << " particles" << std::endl;

    // Calculate the start and end particle for this MPI rank
    auto nper = (N + nrank - 1) / nrank;
    auto beg = nper * irank;
    auto end = nper * (irank + 1);
    if (beg>N) beg = N;
    if (end>N) end = N;
    //Exercise2
    Array<float,2> r_m( Range(beg,end-1), Range(0,3) ); // shape: (N_local,4) for local chunk
        // r => first 3 columns (x,y,z)
    Array<float,2> r = r_m(Range(beg,end-1), Range(0,2));

    // m => fourth column
    Array<float,1> m = r_m(Range(beg,end-1), 3);

    // Array<float,2> r(Range(beg,end-1),Range(0,2));
    // Array<float,1> m(Range(beg,end-1));

    io.seekg( sizeof(tipsy::header) + beg*sizeof(tipsy::dark));
    // Load the particles
    tipsy::dark d;

    //FILL them from the Tipsy file
    for (int i = beg; i < end; i++) {
        if (!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
            perror(argv[1]);
            abort();
        }
        r(i,0) = d.pos[0]; // x
        r(i,1) = d.pos[1]; // y
        r(i,2) = d.pos[2]; // z
        m(i) = d.mass;   // mass
    }
    // for(int i=beg; i<end; ++i) {
    //     if (!io.read(reinterpret_cast<char*>(&d),sizeof(d))) {
    //         perror(argv[1]); abort();
    //     }
    //     r(i,0) = d.pos[0];
    //     r(i,1) = d.pos[1];
    //     r(i,2) = d.pos[2];
    //     m(i) = d.mass;
    // }
    delete [] buffer;

    duration dt = hrc::now() - t0;
    auto rate = (sizeof(tipsy::header) + N*sizeof(tipsy::dark)) / dt.count() / 1024 / 1024;
    if (irank==0) std::cerr << "File reading took " << std::setw(9) << dt.count() << " seconds (" << rate <<" MB/s)." << std::endl;

    // Create Mass Assignment Grid
    t0 = hrc::now();

    auto k_nz = nGrid/2 + 1;
    auto n_floats = size_t(1) * nGrid * nGrid * 2*k_nz; // Careful. For odd nGrid 2*k_nz != nGrid + 2
    float *data = new (std:: align_val_t (64)) float [n_floats]; // 512-bit alignment
    Array<float,3> raw_grid(data,shape(nGrid,nGrid,2*k_nz),deleteDataWhenDone);
    Array<float,3> grid(raw_grid(Range(0,nGrid-1),Range(0,nGrid-1),Range(0,nGrid-1)));

    Array <std::complex<float>,3> kgrid(reinterpret_cast<std::complex<float>*>(data),shape(nGrid,nGrid,k_nz),neverDeleteData);

    // Assign the mass to the grid
    raw_grid = 0;

    // This creates four different versions of "assign mass", one for each order
    if (irank==0) std::cerr << "Assigning mass to the grid using order " << iOrder <<std::endl;
    switch(iOrder) {
    case 1:
        assign_mass<1>(grid,r,m);
        break;
    case 2:
        assign_mass<2>(grid,r,m);
        break;
    case 3:
        assign_mass<3>(grid,r,m);
        break;
    case 4:
        assign_mass<4>(grid,r,m);
        break;
    default:
        std::cerr << "Invalid order " << iOrder << " (must be 1, 2, 3 or 4)" << std::endl;
    }
    dt = hrc::now() - t0;
    if (irank==0) std::cerr << "Mass assignment took " << std::setw(9) << dt.count() << " seconds." << std::endl;
    if (irank==0) std::cerr << "Total mass assigned is " << std::setw(9) << blitz::sum(grid) << std::endl;


    ptrdiff_t local_n0, local_0_start, complex_count; //local_n0:number of slabs this rank owns, local_0_start:the starting x-index for this rank
    complex_count = fftwf_mpi_local_size_3d(nGrid,nGrid,k_nz,MPI_COMM_WORLD,&local_n0,&local_0_start);

    // Exchange information
    int start = local_0_start;
    Array<int,1> rank_start(nrank);
    MPI_Allgather(&start,1,MPI_INT,rank_start.data(),1,MPI_INT,MPI_COMM_WORLD);

    int count = local_n0;
    Array<int,1> rank_count(nrank);
    MPI_Allgather(&count,1,MPI_INT,rank_count.data(),1,MPI_INT,MPI_COMM_WORLD);
    //Exercise1
    // Create an array of size nGrid that tells us which rank owns each x-slice
    Array<int,1> slab_to_rank(nGrid);
    slab_to_rank = -1; // Initialize to -1 just in case

    for (int r = 0; r < nrank; r++) {
        int begin = rank_start(r);
        int end   = rank_start(r) + rank_count(r); // this is one past the last
        for (int i = begin; i < end; i++) {
        slab_to_rank(i) = r;
        }
    }
    //end Exercise 1

    //Exercise2
    auto compareByRank = [&](const FourFloats &A, const FourFloats &B){
        auto slabWrap = [&](float x) {
            int s = int( floor((x + 0.5f)*nGrid) );
            if (s<0) s += nGrid;
            else if (s>=nGrid) s -= nGrid;
            return s;
        };
        int sA = slabWrap(A.x);
        int sB = slabWrap(B.x);
        int rankA = slab_to_rank(sA);
        int rankB = slab_to_rank(sB);
        if (rankA < rankB) return true;  // A goes before B
        if (rankA > rankB) return false; // A goes after B
        return false; // tie => stable or arbitrary
    };
    
    float *basePtr = r_m.dataFirst(); // pointer to first element
    // interpret it as an array of 'FourFloats'
    FourFloats *beginPtr = reinterpret_cast<FourFloats*>(basePtr);
    FourFloats *endPtr   = beginPtr + (end - beg);
    
    // now sort in place
    std::sort(beginPtr, endPtr, compareByRank);
    //Exercise2 end

    //EXERCISE 3
    // Here we create an integer array scounts with one entry per rank:
    Array<int,1> scounts(nrank);
    scounts = 0;  // initialize all counts to zero

    // We can reuse the same slabWrap logic from above:
    auto slabWrapForCounting = [&](float x) {
        int s = int( floor((x + 0.5f)*nGrid) );
        if (s < 0)       s += nGrid;
        else if (s>=nGrid) s -= nGrid;
        return s;
    };

    // Now loop over the local particles:
    for (FourFloats* ptr = beginPtr; ptr != endPtr; ++ptr) {
        float x = ptr->x;
        int slabIndex = slabWrapForCounting(x);
        int ownerRank = slab_to_rank(slabIndex);
        scounts(ownerRank)++;
    }

    for (int who=0; who<nrank; who++) {
        if (irank == who) {
            std::cerr << "EXERCISE 3: scounts array on rank " << irank << ":\n";
            for (int r = 0; r < nrank; r++) {
                std::cerr << "  rank " << r << " = " << scounts(r) << " particles\n";
            }
            std::cerr << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //Exercise3 end

    // EXERCISE 4
    Array<int,1> soffset(nrank);
    soffset = 0;  // Not strictly needed, but good to init

    int runningSum = 0;
    for(int r=0; r<nrank; r++){
        soffset(r) = runningSum;  
        runningSum += scounts(r);
    }
    if (irank == 0) {
        std::cerr << "EXERCISE 4: soffset array on rank " << irank << ":\n";
        for (int rr = 0; rr < nrank; rr++){
            std::cerr << "  rank " << rr << " offset = " << soffset(rr) << std::endl;
        }
    }
    //Exercise4 end
    //Exercise5
    for (int r=0; r<nrank; r++) {
        int startIdx = soffset(r);
        int num = scounts(r);
        for (int i = 0; i < num; i++) {
            int idx = startIdx + i;
            FourFloats &p = *(beginPtr + idx);
    
            int s = slabWrapForCounting(p.x);
            int actualRank = slab_to_rank(s);
    
            if (actualRank != r) {
                std::cerr << "ERROR on rank " << irank 
                          << ": local idx " << idx
                          << " => x= " << p.x
                          << " belongs to rank " << actualRank
                          << " but in rank " << r << " segment.\n";
            }
        }
    }
     
    //Exercise5 end

    //Exercise6
    Array<int,1> rcounts(nrank);
    rcounts = 0;

    // 1) Exchange scounts so each rank knows how many particles
    //    it will receive FROM each other rank.
    MPI_Alltoall(
        scounts.data(),     // send buffer
        1,                  // send 1 int to each rank
        MPI_INT,            // type
        rcounts.data(),     // receive buffer
        1,                  // receive 1 int from each rank
        MPI_INT,
        MPI_COMM_WORLD
    );

    // 2) Exclusive scan on rcounts
    Array<int,1> rdispls(nrank);
    rdispls = 0;
    int runningSum2 = 0;
    for(int r=0; r<nrank; r++){
        rdispls(r) = runningSum2;
        runningSum2 += rcounts(r);
    }

    if (irank == 0) {
        std::cerr << "EXERCISE 6: rcounts & rdispls on rank 0:\n";
        for (int rr = 0; rr < nrank; rr++) {
            std::cerr << "  rank " << rr 
                      << " rcounts=" << rcounts(rr) 
                      << " rdispls=" << rdispls(rr) << "\n";
        }
    }
    //Exericse6 end
    //Exercise7  
    // We want a new array r_m2 to receive the incoming data
    // total number of incoming particles = sum(rcounts)
    int totalRecv = runningSum2;  // from the final runningSum2
    // each of those is 1 "particle" = 4 floats (x,y,z,m)

    Array<float,2> r_m2( Range(0, totalRecv-1), Range(0,3) );
    Array<float,2> r2 = r_m2( Range(0,totalRecv-1), Range(0,2) );
    Array<float,1> m2 = r_m2( Range(0,totalRecv-1), 3 );

    float* recvPtr = r_m2.dataFirst(); 

    // Our local data pointer is beginPtr for the old array:
    float* sendPtr = reinterpret_cast<float*>(beginPtr); 

    // 2) We must also compute float-based counts & displacements
    Array<int,1> sendCountsFloat(nrank);
    Array<int,1> sendDisplsFloat(nrank);
    Array<int,1> recvCountsFloat(nrank);
    Array<int,1> recvDisplsFloat(nrank);

    for(int r=0; r<nrank; r++){
        sendCountsFloat(r) = scounts(r)*4;

        sendDisplsFloat(r) = soffset(r)*4;
        recvCountsFloat(r) = rcounts(r)*4;
        recvDisplsFloat(r) = rdispls(r)*4;
    }

    // 3) Now do MPI_Alltoallv with the float-based counts
    MPI_Alltoallv(
        // send
        sendPtr,
        sendCountsFloat.data(),
        sendDisplsFloat.data(),
        MPI_FLOAT,
        // receive
        recvPtr,
        recvCountsFloat.data(),
        recvDisplsFloat.data(),
        MPI_FLOAT,
        MPI_COMM_WORLD
    );
    //Exercise 7 end

    //Exercise 8
    bool anyError = false;
    for (int i = 0; i < totalRecv; i++) {
        float x = r2(i,0); // or r_m2(i,0)
        int s = slabWrapForCounting(x);
        int actualRank = slab_to_rank(s); // from earlier exercises
        if (actualRank != irank) {
            std::cerr << "EX8 ERROR on rank " << irank
                      << ": local idx " << i
                      << " => x=" << x
                      << " belongs to rank " << actualRank
                      << " but is on rank " << irank << "\n";
            anyError = true;
        }
    }
    if (!anyError) {
        std::cerr << "Rank " << irank << " (EX8): all new particles verified, no errors.\n";
    }
    //Exercise8 end
    //Exercise9
    // 1) Create a new array for the mass assignment
    size_t n_floats2 = size_t(nGrid)*nGrid*2*k_nz;

    float *data2 = new (std::align_val_t (64)) float[n_floats2]; 
    Array<float,3> raw_grid2(data2, shape(nGrid,nGrid,2*k_nz), deleteDataWhenDone);
    raw_grid2 = 0.0f; // zero it
    // "grid2" is the subarray [0..nGrid-1, 0..nGrid-1, 0..nGrid-1]
    Array<float,3> grid2(raw_grid2(Range(0,nGrid-1), Range(0,nGrid-1), Range(0,nGrid-1)));

    // 2) Assign mass using r2,m2
    switch(iOrder) {
        case 1:
        assign_mass<1>(grid2, r2, m2);
        break;
        case 2:
        assign_mass<2>(grid2, r2, m2);
        break;
        case 3:
        assign_mass<3>(grid2, r2, m2);
        break;
        case 4:
        assign_mass<4>(grid2, r2, m2);
        break;
        default:
        std::cerr << "Invalid order for second assignment: " << iOrder << "\n";
    }

    // 3)compare local sums vs the old grid
    //For a global sum:
    float local_mass2 = blitz::sum(grid2);
    float total_mass2;
    MPI_Allreduce(&local_mass2, &total_mass2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    if (irank==0) {
        std::cerr << "EXERCISE 9: total mass assigned with r2,m2 is " 
                    << total_mass2 << std::endl
                    << "Compare to old total mass which was printed as " 
                    << "somewhere above in logs." << std::endl;
    }
    //Exercise9 end

    // Create the local slabs
    float *slab_data = new (std:: align_val_t (64)) float [2*complex_count]; // 512-bit alignment
    Array<float,3> raw_slab(slab_data,shape(local_n0,nGrid,2*k_nz),deleteDataWhenDone);
    raw_slab.reindexSelf(TinyVector<int,3>(local_0_start,0,0)); // Change the start index of the first dimension
    Array<float,3> slab(raw_slab(Range(local_0_start,local_0_start+local_n0-1),Range(0,nGrid-1),Range(0,nGrid-1)));

    Array <std::complex<float>,3> kslab(reinterpret_cast<std::complex<float>*>(slab_data),shape(local_n0,nGrid,k_nz),neverDeleteData);
    kslab.reindexSelf(TinyVector<int,3>(local_0_start,0,0)); // Change the start index of the first dimension

    // Create the FFT plan. Generally this should be done BEFORE you set the data as "MEASURE" will use the buffer.
#ifdef _OPENMP
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
    auto plan = fftwf_mpi_plan_dft_r2c_3d(nGrid,nGrid,nGrid,
            slab.dataFirst(),
            reinterpret_cast<fftwf_complex*>(kslab.dataFirst()),
            MPI_COMM_WORLD,FFTW_ESTIMATE);

    // Combine the grids into the distibuted slabs
    int slab_size = nGrid * 2*k_nz;
    for(auto root=0; root<nrank; ++root) {
        MPI_Reduce(&grid(rank_start(root),0,0),slab.data(),slab_size * rank_count(root),
            MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
    }

    // Calculate deleta
    float local_mass = blitz::sum(slab);
    float total_mass;
    MPI_Allreduce(&local_mass, &total_mass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    float diRhoBar = ((1.0f*nGrid*nGrid*nGrid)/total_mass);
    slab = slab * diRhoBar - 1.0f;
    // Correct for FFT normalization
    slab /= (nGrid*nGrid*nGrid);

    // Calculate the FFT
    t0 = hrc::now();
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    dt = hrc::now() - t0;
    if (irank==0) std::cerr << "FFT took " << std::setw(9) << dt.count() << " seconds." << std::endl;

    Array<double,1> ak(nBins);  // Sum of k values in each bin
    Array<double,1> pk(nBins);  // Sum of P(k) values in each bin
    Array<long,1> nk(nBins);  // Number of values in each bin

    ak = 0;
    pk = 0;
    nk = 0;

    if (irank==0) std::cerr << "Using " << nBins << (bLog?" logarithmic":" linear") << " bins." << std::endl;

    AssignmentWindow W(nGrid,iOrder);
    for(auto ii=kslab.begin(); ii!=kslab.end(); ++ii) {
        auto pos = ii.position();
        auto bin = [iNyquist,nGrid](int k) {return k<=iNyquist ? k : k-nGrid;};
        auto kx = bin(pos[0]);
        auto ky = bin(pos[1]);
        auto kz = pos[2];
        float k = sqrt(kx*kx + ky*ky + kz*kz);
        *ii *= W[std::abs(kx)] * W[std::abs(ky)] * W[kz]; // Correction for mass assignment
        int i;
        if (bLog) i = int(log(k) / log(iNyquist) * nBins);
        else  i = int(k / iNyquist * nBins);
        if (i>0 && i < nBins) {
            ak(i) += k;
            pk(i) += std::norm(*ii);
            nk(i) += 1;
        }
    }

    Array<double,1> sum_ak(nBins);
    Array<double,1> sum_pk(nBins);
    Array<long,1> sum_nk(nBins);

    MPI_Reduce(ak.data(), sum_ak.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pk.data(), sum_pk.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(nk.data(), sum_nk.data(), nBins, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (irank==0)  {
        sum_ak = sum_ak / sum_nk;
        sum_pk = sum_pk / sum_nk;
        for(auto i=1; i<nBins; ++i) {
            if (sum_nk(i)) {
                printf("%.10g %.10g %ld\n", sum_ak(i), sum_pk(i),sum_nk(i));
            }
        }
    }

    MPI_Finalize();
}
