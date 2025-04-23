// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <new>
#include <complex>
#include <cmath>
#include <stdlib.h>
#include <vector>
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


struct FourFloats {
    float x, y, z, m;
};

// -----------------------------------------------------------------------------
// EXERCISE 1 & 2: A mass-assignment routine that writes directly
// into the "raw_slab" array (including ghost cells) in a distributed manner.
// -----------------------------------------------------------------------------
template<int Order=1>
void assign_mass_distributed(
    Array<float,3> &slab_with_ghost,  // local array (with ghost region)
    int             nGrid,
    int             local_0_start,
    int             local_n0,
    int             ghost_slabs,
    Array<float,2> &R,  // positions (x,y,z)
    Array<float,1> &M   // masses
)
{
    // The X-range in the ghosted slab is:
    // [ local_0_start - ghost_slabs .. local_0_start + local_n0 + ghost_slabs - 1 ]
    int x_min = local_0_start - ghost_slabs;
    int x_max = local_0_start + local_n0 + ghost_slabs - 1; // inclusive

    auto local_xwrap = [&](int i) {
        if (i<0)       i += nGrid;
        else if (i>=nGrid) i -= nGrid;
        while (i < x_min) i += nGrid;
        while (i > x_max) i -= nGrid;
        return i;
    };

    // Normal periodic wrap for y,z in [0..nGrid-1]
    auto yz_wrap = [&](int j) {
        if (j<0)       j += nGrid;
        else if (j>=nGrid) j -= nGrid;
        return j;
    };

    // Loop over local subset of particles
    #pragma omp parallel for
    for(int pn = R.lbound(0); pn<=R.ubound(0); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);

        // Compute continuous indices in [0..nGrid)
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),
                                       Hy((y+0.5f)*nGrid),
                                       Hz((z+0.5f)*nGrid);

        // Add mass to up to Order^3 cells
        for(int i=0; i<Order; i++){
            for(int j=0; j<Order; j++){
                for(int k=0; k<Order; k++){
                    int ix = local_xwrap(Hx.i + i);
                    int iy = yz_wrap   (Hy.i + j);
                    int iz = yz_wrap   (Hz.i + k);

                    #pragma omp atomic
                    slab_with_ghost(ix, iy, iz) += m * (Hx.H[i]*Hy.H[j]*Hz.H[k]);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// For EXERCISE 3: We do a ring reduce of ghost cells. We'll just show placeholders
// for extracting and merging ghost data. In a real code you'd figure out exactly
// how each slab's ghost region maps onto the "next" rank's domain.
// -----------------------------------------------------------------------------
void extract_ghost_data(
    Array<float,3> &raw_slab,  
    int local_0_start,
    int local_n0,
    int ghost_slabs,
    int nGrid,
    std::vector<float> &ghostBuf
)
{
    int x_lo_start = local_0_start - ghost_slabs;
    int x_lo_end   = local_0_start - 1;
    int x_hi_start = local_0_start + local_n0;
    int x_hi_end   = local_0_start + local_n0 + ghost_slabs - 1;

    int sliceDepth = raw_slab.extent(2); 
    ghostBuf.clear();
    ghostBuf.reserve( std::size_t(2*ghost_slabs)*nGrid*sliceDepth );

    auto copySlice = [&](int ix) {
        for(int iy=0; iy<nGrid; iy++){
            for(int iz=0; iz<sliceDepth; iz++){
                ghostBuf.push_back( raw_slab(ix, iy, iz) );
            }
        }
    };
    // low side
    for(int ix=x_lo_start; ix<=x_lo_end; ix++){
        copySlice(ix);
    }
    // high side
    for(int ix=x_hi_start; ix<=x_hi_end; ix++){
        copySlice(ix);
    }
}

void merge_ghost_data(
    Array<float,3> &raw_slab,
    int local_0_start,
    int local_n0,
    int ghost_slabs,
    int nGrid,
    const std::vector<float> &ghostBuf
)
{
    // We'll do a trivial demonstration that merges onto the low side x-range.
    int x_lo_start = local_0_start;
    int x_lo_end   = local_0_start + ghost_slabs - 1;

    int sliceDepth = raw_slab.extent(2);
    size_t pos=0;
    for(int ix=x_lo_start; ix<=x_lo_end; ix++){
        for(int iy=0; iy<nGrid; iy++){
            for(int iz=0; iz<sliceDepth; iz++){
                raw_slab(ix, iy, iz) += ghostBuf[pos++];
            }
        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int irank=0, nrank=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);

#ifdef _OPENMP
    fftwf_init_threads();
#endif

    if(argc<=2) {
        if(irank==0) {
            std::cerr << "Usage: " << argv[0]
                      << " tipsyfile.std grid-size[/[L]bin-count] [order]\n";
        }
        MPI_Finalize();
        return 1;
    }

    // parse nGrid and optional bin info
    int nGrid = std::atoi(argv[2]);
    int iNyquist = nGrid/2;
    int nBins = iNyquist;
    bool bLog = false;
    char* slash = std::strchr(argv[2],'/');
    if(slash) {
        *slash++ = '\0';
        if(*slash=='L') {
            bLog=true; slash++;
        }
        nBins = std::atoi(slash);
    }

    int iOrder = 1;
    if(argc>3) iOrder=std::atoi(argv[3]);

    // -------------------------------------------------------------------------
    // Read Tipsy file and store local subset in zero-based array
    // -------------------------------------------------------------------------
    auto t0 = hrc::now();
    std::ifstream fin(argv[1], std::ios::binary);
    if(!fin) {
        std::cerr<<"Cannot open file "<< argv[1] <<"\n";
        MPI_Finalize();
        return 1;
    }

    tipsy::header H;
    if(!fin.read(reinterpret_cast<char*>(&H), sizeof(H))) {
        std::cerr<<"Error reading tipsy header\n";
        MPI_Finalize();
        return 1;
    }
    std::uint64_t N = H.nDark;
    if(irank==0) {
        std::cerr<<"Reading "<< N <<" dark particles.\n";
    }

    // distribute among ranks
    std::uint64_t nper = (N + nrank -1)/nrank;
    std::uint64_t beg  = nper*irank;
    std::uint64_t end  = nper*(irank+1);
    if(beg>N) beg=N;
    if(end>N) end=N;
    int localCount = (int)(end-beg);  // should be safe if end-beg <= INT_MAX

    // We'll store local data in a zero-based array of shape (localCount,4).
    // columns: 0->x, 1->y, 2->z, 3->mass
    blitz::Array<float,2> r_m(localCount,4);
    // define "r" => columns [0..2], "mA" => column 3
    Array<float,2> r  = r_m(Range::all(), Range(0,2));
    Array<float,1> mA = r_m(Range::all(), 3);

    // position file to read from beg-th dark particle
    fin.seekg(sizeof(tipsy::header) + beg*sizeof(tipsy::dark));
    tipsy::dark dpart;
    for(int i=0; i<localCount; i++){
        if(!fin.read(reinterpret_cast<char*>(&dpart), sizeof(dpart))) {
            perror("Reading tipsy");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        r(i,0) = dpart.pos[0];
        r(i,1) = dpart.pos[1];
        r(i,2) = dpart.pos[2];
        mA(i)  = dpart.mass;
    }
    fin.close();

    duration dt = hrc::now() - t0;
    if(irank==0) {
        double MB = (sizeof(tipsy::header) + N*sizeof(tipsy::dark))/1024./1024.;
        std::cerr<<"File read took "<< dt.count()
                 <<" s ("<< MB/dt.count() <<" MB/s)\n";
    }

    ptrdiff_t local_n0, local_0_start, complex_count;
    auto k_nz = nGrid/2 + 1;
    complex_count = fftwf_mpi_local_size_3d(
                        nGrid,nGrid,k_nz,
                        MPI_COMM_WORLD,
                        &local_n0,
                        &local_0_start);

    // For EX1 & EX2: create slab + ghost region
    int ghost_slabs = (iOrder>1)? (iOrder-1) : 0;
    int slab_count  = local_n0 + 2*ghost_slabs;
    size_t needed_for_ghost = size_t(slab_count)*nGrid*(2*k_nz);
    size_t needed_for_fftw  = size_t(2)*complex_count;
    size_t final_needed     = std::max(needed_for_ghost, needed_for_fftw);

    float *slab_data = new (std::align_val_t(64)) float[final_needed];
    Array<float,3> raw_slab(
        slab_data,
        shape(slab_count, nGrid, 2*k_nz),
        deleteDataWhenDone
    );
    // reindex so the "owned" region is [ local_0_start .. local_0_start+local_n0 -1 ]
    raw_slab.reindexSelf( TinyVector<int,3>(local_0_start-ghost_slabs, 0, 0) );
    raw_slab = 0.f;

    // subarray with no ghosts
    Array<float,3> slab = raw_slab(
        Range(local_0_start, local_0_start+local_n0-1),
        Range(0,nGrid-1),
        Range(0,2*k_nz-1)
    );

    // Complex view for FFT
    Array<std::complex<float>,3> kslab(
        reinterpret_cast<std::complex<float>*>(slab_data),
        shape(slab_count, nGrid, k_nz),
        neverDeleteData
    );
    kslab.reindexSelf( TinyVector<int,3>(local_0_start-ghost_slabs, 0, 0) );
    Array<std::complex<float>,3> kslab_noghost = kslab(
        Range(local_0_start, local_0_start+local_n0-1),
        Range(0, nGrid-1),
        Range(0, k_nz-1)
    );

    // Assign mass
    t0 = hrc::now();
    switch(iOrder){
        case 1: assign_mass_distributed<1>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 2: assign_mass_distributed<2>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 3: assign_mass_distributed<3>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 4: assign_mass_distributed<4>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        default:
            std::cerr<<"Unsupported iOrder="<< iOrder <<"\n";
    }
    dt=hrc::now()-t0;
    float local_sum = sum(slab);
    float total_sum=0.f;
    MPI_Allreduce(&local_sum,&total_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    if(irank==0){
        std::cerr<<"Mass assignment took "<< dt.count()
                 <<" s. sum="<< total_sum <<"\n";
    }

    // EXERCISE 3 ring reduce of ghost data
    // Suppose rank i's ghost region belongs to rank i+1 mod nrank.
    if(nrank%2!=0 && irank==0){
        std::cerr<<"Note: ring reduce code assumes even nrank.\n";
    }
    std::vector<float> ghostBuf;
    extract_ghost_data(raw_slab,local_0_start,local_n0,ghost_slabs,nGrid, ghostBuf);

    // 1st split: (0,1), (2,3), ...
    int color1 = irank/2;
    int key1   = irank%2;
    MPI_Comm newcomm1;
    MPI_Comm_split(MPI_COMM_WORLD, color1, key1, &newcomm1);
    int newrank1; MPI_Comm_rank(newcomm1, &newrank1);

    // rank with newrank=1 is the "root"
    std::vector<float> recvBuf1(ghostBuf.size(), 0.f);
    MPI_Request req1;
    MPI_Ireduce( ghostBuf.data(), recvBuf1.data(),
                 (int)ghostBuf.size(),
                 MPI_FLOAT,
                 MPI_SUM,
                 1, // root
                 newcomm1,
                 &req1 );

    // 2nd split: (1,2), (3,4), (5,0), ...
    int next = (irank+1)%nrank;
    int color2= (irank<next)? irank:next;
    int key2  = (irank==next)?0 : (irank==color2?0:1);
    MPI_Comm newcomm2;
    MPI_Comm_split(MPI_COMM_WORLD, color2, key2, &newcomm2);
    int newrank2=-1;
    if(newcomm2!=MPI_COMM_NULL) {
        MPI_Comm_rank(newcomm2,&newrank2);
    }
    std::vector<float> recvBuf2(ghostBuf.size(),0.f);
    MPI_Request req2 = MPI_REQUEST_NULL;
    if(newcomm2!=MPI_COMM_NULL){
        MPI_Ireduce( ghostBuf.data(), recvBuf2.data(),
                     (int)ghostBuf.size(),
                     MPI_FLOAT, MPI_SUM,
                     1, newcomm2, &req2 );
    }

    // wait
    MPI_Request reqs[2] = {req1, req2};
    MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

    // merge
    if(newrank1==1){
        merge_ghost_data(raw_slab, local_0_start, local_n0, ghost_slabs, nGrid, recvBuf1);
    }
    if(newcomm2!=MPI_COMM_NULL && newrank2==1){
        merge_ghost_data(raw_slab, local_0_start, local_n0, ghost_slabs, nGrid, recvBuf2);
    }
    if(newcomm1!=MPI_COMM_NULL) MPI_Comm_free(&newcomm1);
    if(newcomm2!=MPI_COMM_NULL) MPI_Comm_free(&newcomm2);

    // EXERCISE 4
    local_sum = sum(slab);
    MPI_Allreduce(&local_sum, &total_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    if(irank==0){
        std::cerr<<"After ring reduce, total mass="<< total_sum <<"\n";
    }

#ifdef _OPENMP
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
    auto plan = fftwf_mpi_plan_dft_r2c_3d(
        nGrid, nGrid, nGrid,
        slab.dataFirst(),
        reinterpret_cast<fftwf_complex*>(kslab_noghost.dataFirst()),
        MPI_COMM_WORLD, FFTW_ESTIMATE
    );

    // Convert to overdensity = (rho/rho_bar) -1, then scale by 1/(nGrid^3)
    float diRhoBar = float(nGrid)*nGrid*nGrid / total_sum;
    slab = slab*diRhoBar - 1.f;
    slab /= (nGrid*nGrid*nGrid);

    t0 = hrc::now();
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    dt = hrc::now()-t0;
    if(irank==0){
        std::cerr<<"FFT took "<< dt.count() <<" s.\n";
    }

    // compute power spectrum
    Array<double,1> ak(nBins), pk(nBins);
    Array<long,1>   nk(nBins);
    ak=0.; pk=0.; nk=0;

    AssignmentWindow W(nGrid,iOrder);

    for(auto it = kslab_noghost.begin(); it!=kslab_noghost.end(); ++it){
        auto pos = it.position();
        int iX=pos[0], iY=pos[1], iZ=pos[2];
        auto fold=[&](int k){return (k<=iNyquist)?k:(k-nGrid);};
        int kx=fold(iX), ky=fold(iY), kz=iZ;
        float kk = std::sqrt(kx*kx + ky*ky + kz*kz);

        (*it)*= W[std::abs(kx)]*W[std::abs(ky)]*W[kz];

        int ibin;
        if(bLog){
            if(kk<=0.f) continue;
            ibin = int(std::log(kk)/std::log(iNyquist)*nBins);
        } else {
            ibin = int(kk/float(iNyquist)*nBins);
        }
        if(ibin>0 && ibin<nBins){
            ak(ibin)+= kk;
            pk(ibin)+= std::norm(*it);
            nk(ibin)++;
        }
    }

    Array<double,1> sum_ak(nBins), sum_pk(nBins);
    Array<long,1>   sum_nk(nBins);
    sum_ak=0.; sum_pk=0.; sum_nk=0;
    MPI_Reduce(ak.data(), sum_ak.data(), nBins, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(pk.data(), sum_pk.data(), nBins, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(nk.data(), sum_nk.data(), nBins, MPI_LONG,   MPI_SUM, 0,MPI_COMM_WORLD);

    if(irank==0){
        for(int i=1; i<nBins; i++){
            if(sum_nk(i)>0){
                double meanK = sum_ak(i)/sum_nk(i);
                double meanP = sum_pk(i)/sum_nk(i);
                std::cout<< meanK <<" "<< meanP <<" "<< sum_nk(i) <<"\n";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
