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
#include <algorithm>

#include "mpi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "blitz/array.h"
#include "fftw3-mpi.h"

// Your local headers for reading TIPSY and assignment weights:
#include "tipsy.h"
#include "aweights.h"

using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

struct FourFloats {
    float x, y, z, m;
};

// ------------------------------------------------------------------
// Utility: Assign masses onto a slab-with-ghosts
// ------------------------------------------------------------------
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
    int x_min = local_0_start - ghost_slabs;
    int x_max = local_0_start + local_n0 + ghost_slabs - 1; // inclusive

    auto local_xwrap = [&](int i) {
        if (i<0)       i += nGrid;
        else if (i>=nGrid) i -= nGrid;
        while (i < x_min) i += nGrid;
        while (i > x_max) i -= nGrid;
        return i;
    };

    auto yz_wrap = [&](int j) {
        if (j<0)       j += nGrid;
        else if (j>=nGrid) j -= nGrid;
        return j;
    };

#pragma omp parallel for
    for(int pn = R.lbound(0); pn<=R.ubound(0); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);

        // Use assignment weights from your aweights.h
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),
                                       Hy((y+0.5f)*nGrid),
                                       Hz((z+0.5f)*nGrid);

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

// ------------------------------------------------------------------
// Utility: Extract ghost data from the real-slab
// ------------------------------------------------------------------
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

// ------------------------------------------------------------------
// Utility: Merge ghost data
// ------------------------------------------------------------------
void merge_ghost_data(
    Array<float,3> &raw_slab,
    int local_0_start,
    int local_n0,
    int ghost_slabs,
    int nGrid,
    const std::vector<float> &ghostBuf
)
{
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
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif

    // -------------------------------------------------------
    // Parse command line:
    //   usage: assign [tipsyfilename] [nGrid[/[L]bin-count]] [order]
    // -------------------------------------------------------
    if(argc<=2) {
        if(irank==0) {
            std::cerr << "Usage: " << argv[0]
                      << " tipsyfile.std grid-size[/[L]bin-count] [order]\n";
        }
        MPI_Finalize();
        return 1;
    }

    int nGrid = std::atoi(argv[2]);
    int iNyquist = nGrid/2;
    int nBins = iNyquist;
    bool bLog=false;
    char* slash = std::strchr(argv[2],'/');
    if(slash) {
        *slash++ = '\0';
        if(*slash=='L') {
            bLog=true; 
            slash++;
        }
        nBins = std::atoi(slash);
    }

    int iOrder = 1;
    if(argc>3) {
        iOrder = std::atoi(argv[3]);
    }

    // -------------------------------------------------------
    // Step 1: Read Tipsy file
    // -------------------------------------------------------
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
        std::cerr<<"Reading "<< N <<" dark particles from "<< argv[1] <<"\n";
    }

    // distribute among ranks
    std::uint64_t nper = (N + nrank -1)/nrank;
    std::uint64_t beg  = nper*irank;
    std::uint64_t end  = nper*(irank+1);
    if(beg>N) beg=N;
    if(end>N) end=N;
    int localCount = (int)(end-beg);  

    // We'll store local data in an array of shape (localCount,4).
    blitz::Array<float,2> r_m(localCount,4);
    Array<float,2> r  = r_m(Range::all(), Range(0,2)); // positions x,y,z
    Array<float,1> mA = r_m(Range::all(), 3);          // masses

    // Seek to beg-th dark particle
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

    // -------------------------------------------------------
    // Step 2: MPI distribution for a 3D array (nGrid^3).
    // -------------------------------------------------------
    ptrdiff_t local_n0, local_0_start;
    ptrdiff_t local_n1, local_1_start; 
    ptrdiff_t complex_count;

    complex_count = fftwf_mpi_local_size_3d_transposed(
        nGrid, nGrid, nGrid,    
        MPI_COMM_WORLD,
        &local_n0, &local_0_start,     
        &local_n1, &local_1_start    
    );

    ptrdiff_t k_nz = nGrid/2 + 1; // real->complex dimension in z

    if(irank==0){
        std::cerr << "Rank=" << irank
                  << " local_n0=" << local_n0
                  << " local_0_start=" << local_0_start
                  << " local_n1=" << local_n1
                  << " local_1_start=" << local_1_start
                  << "\n";
    }

    // Create local slab array with ghost region on x dimension
    int ghost_slabs = (iOrder>1) ? (iOrder -1) : 0;
    int slab_count  = local_n0 + 2*ghost_slabs; 

    // We need space for either the "ghosted slab" or "2*complex_count floats"
    size_t needed_for_ghost = size_t(slab_count)*nGrid*(2*k_nz);
    size_t needed_for_fftw  = size_t(2)*complex_count;
    size_t final_needed     = std::max(needed_for_ghost, needed_for_fftw);

    float * slab_data = new (std::align_val_t(64)) float[final_needed];

    Array<float,3> raw_slab(
        slab_data,
        shape(slab_count, nGrid, 2*k_nz),
        deleteDataWhenDone
    );
    raw_slab.reindexSelf( TinyVector<int,3>( local_0_start-ghost_slabs, 0, 0) );
    raw_slab = 0.f;

    // subarray with no ghosts
    Array<float,3> slab = raw_slab(
        Range(local_0_start, local_0_start+local_n0-1),
        Range(0,nGrid-1),
        Range(0,2*k_nz-1)
    );

    // -------------------------------------------------------
    // Step 3: Assign mass (NGP/CIC/etc).
    // -------------------------------------------------------
    t0 = hrc::now();
    switch(iOrder){
        case 1: assign_mass_distributed<1>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 2: assign_mass_distributed<2>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 3: assign_mass_distributed<3>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        case 4: assign_mass_distributed<4>(raw_slab,nGrid,local_0_start,local_n0,ghost_slabs, r, mA); break;
        default:
            if(irank==0){
                std::cerr<<"Unsupported iOrder="<< iOrder <<"\n";
            }
            MPI_Abort(MPI_COMM_WORLD,1);
    }
    dt=hrc::now()-t0;
    float local_sum = sum(slab);
    float total_sum=0.f;
    MPI_Allreduce(&local_sum,&total_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    if(irank==0){
        std::cerr<<"Mass assignment took "<< dt.count()
                 <<" s. sum="<< total_sum <<"\n";
    }

    // -------------------------------------------------------
    // Step 4: Merge ghost zones
    // -------------------------------------------------------
    {
        std::vector<float> ghostBuf;
        extract_ghost_data(raw_slab, local_0_start, local_n0, ghost_slabs, nGrid, ghostBuf);

        int color1 = irank/2;
        int key1   = irank%2;
        MPI_Comm newcomm1;
        MPI_Comm_split(MPI_COMM_WORLD, color1, key1, &newcomm1);
        int newrank1; 
        MPI_Comm_rank(newcomm1, &newrank1);

        std::vector<float> recvBuf1(ghostBuf.size(), 0.f);
        MPI_Request req1;
        MPI_Ireduce(
            ghostBuf.data(), recvBuf1.data(),
            (int)ghostBuf.size(),
            MPI_FLOAT, MPI_SUM,
            1, newcomm1, &req1
        );

        int next = (irank+1)%nrank;
        int color2= (irank<next)? irank:next;
        int key2  = (irank==color2?0:1);
        MPI_Comm newcomm2;
        MPI_Comm_split(MPI_COMM_WORLD, color2, key2, &newcomm2);
        int newrank2=-1;
        if(newcomm2!=MPI_COMM_NULL) {
            MPI_Comm_rank(newcomm2,&newrank2);
        }
        std::vector<float> recvBuf2(ghostBuf.size(),0.f);
        MPI_Request req2 = MPI_REQUEST_NULL;
        if(newcomm2!=MPI_COMM_NULL){
            MPI_Ireduce(
                ghostBuf.data(), recvBuf2.data(),
                (int)ghostBuf.size(),
                MPI_FLOAT, MPI_SUM,
                1, newcomm2, &req2 
            );
        }
        MPI_Request reqs[2] = {req1, req2};
        MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

        if(newrank1==1){
            merge_ghost_data(raw_slab, local_0_start, local_n0, ghost_slabs, nGrid, recvBuf1);
        }
        if(newcomm2!=MPI_COMM_NULL && newrank2==1){
            merge_ghost_data(raw_slab, local_0_start, local_n0, ghost_slabs, nGrid, recvBuf2);
        }
        if(newcomm1!=MPI_COMM_NULL) MPI_Comm_free(&newcomm1);
        if(newcomm2!=MPI_COMM_NULL) MPI_Comm_free(&newcomm2);

        local_sum = sum(slab);
        MPI_Allreduce(&local_sum,&total_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        if(irank==0){
            std::cerr<<"After ring reduce, total mass="<< total_sum <<"\n";
        }
    }


    //Exercise4: 2D real->complex plan

    fftwf_plan plan_2d = nullptr;
    {
        int Nx2D = nGrid;
        int Ny2D = nGrid;
        // Temporary arrays just to create the plan
        std::vector<float> slab2D_in ( Nx2D * Ny2D, 0.f );
        std::vector<std::complex<float>> slab2D_out( Nx2D*(Ny2D/2 +1) );

        plan_2d = fftwf_plan_dft_r2c_2d(
            Nx2D, Ny2D,
            slab2D_in.data(),
            reinterpret_cast<fftwf_complex*>( slab2D_out.data() ),
            FFTW_ESTIMATE
        );
    }

    //Exercise5: Manual Transpose w/ MPI

    fftwf_plan plan_transpose = nullptr;
    ptrdiff_t local_n0_t, local_0_start_t;
    ptrdiff_t local_n1_t, local_1_start_t;
    fftwf_complex * transpose_in  = nullptr;
    fftwf_complex * transpose_out = nullptr;
    {
        ptrdiff_t howmany = k_nz;
        ptrdiff_t n2d[2] = { nGrid, nGrid };
        ptrdiff_t block0=1, block1=1;

        // Get local size in the transposed distribution
        ptrdiff_t local_size_transposed = fftwf_mpi_local_size_many_transposed(
            2,            // rank=2
            n2d,
            howmany,
            block0,
            block1,
            MPI_COMM_WORLD,
            &local_n0_t, &local_0_start_t,
            &local_n1_t, &local_1_start_t
        );

        // We allocate in/out once, for the final transpose
        // local_size_pre = local_n0 * nGrid * k_nz
        // local_size_post= local_n0_t * nGrid * k_nz
        ptrdiff_t local_size_pre  = local_n0 * nGrid * k_nz;
        ptrdiff_t local_size_post = local_n0_t * nGrid * k_nz;

        // We must allocate the max needed if plan is out-of-place
        ptrdiff_t max_needed = std::max(local_size_pre, local_size_post);

        transpose_in  = fftwf_alloc_complex( max_needed );
        transpose_out = fftwf_alloc_complex( max_needed );

        // Create an out-of-place plan referencing these pointers:
        plan_transpose = fftwf_mpi_plan_many_transpose(
            nGrid,       // n0
            nGrid,       // n1
            howmany,     // howmany
            block0,
            block1,
            reinterpret_cast<float*>(transpose_in),
            reinterpret_cast<float*>(transpose_out),
            MPI_COMM_WORLD,
            FFTW_ESTIMATE
        );
    }

    //Exercise6: Pencil 1D complex->complex plan
  
    fftwf_plan plan_1d = nullptr;
    {
        constexpr int rank_1d = 1;
        int n_1d[1] = { nGrid };

        int howmany_1d = int(local_n0_t * k_nz);

        int stride_1d = 1;
        int dist_1d   = nGrid; 

        // Just to create the plan (dummy arrays):
        fftwf_complex* tmp_in  = fftwf_alloc_complex(howmany_1d*nGrid);
        fftwf_complex* tmp_out = fftwf_alloc_complex(howmany_1d*nGrid);

        plan_1d = fftwf_plan_many_dft(
            rank_1d,
            n_1d,
            howmany_1d,
            tmp_in,
            n_1d,
            stride_1d,
            dist_1d,
            tmp_out,
            n_1d,
            stride_1d,
            dist_1d,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        fftwf_free(tmp_in);
        fftwf_free(tmp_out);
    }

    // Exercise7: Replace single 3D FFT call with loops + transpose + loops

    // (A) 2D transforms on real slabs (r2c)
    {
        int Nx2D = nGrid;
        int Ny2D = nGrid;
        std::vector<float> slab2D_in ( Nx2D * Ny2D );
        std::vector<std::complex<float>> slab2D_out( Nx2D*(Ny2D/2 +1) );

        auto t2d0 = hrc::now();

        // Loop over local x-slab
        for(int ix = local_0_start; ix < local_0_start + local_n0; ix++){
            // copy from raw_slab => slab2D_in
            for(int iy=0; iy<nGrid; iy++){
                for(int iz=0; iz<nGrid; iz++){
                    slab2D_in[iy*Ny2D + iz] = raw_slab(ix, iy, iz);
                }
            }
            // now do the transform: new-array execute
            fftwf_execute_dft_r2c(
                plan_2d,
                slab2D_in.data(),
                reinterpret_cast<fftwf_complex*>(slab2D_out.data())
            );

            // store back to raw_slab
            for(int iy=0; iy<nGrid; iy++){
                for(int iz=0; iz<(nGrid/2+1); iz++){
                    float re = slab2D_out[iy*(nGrid/2+1) + iz].real();
                    float im = slab2D_out[iy*(nGrid/2+1) + iz].imag();
                    raw_slab(ix, iy, 2*iz  ) = re;
                    raw_slab(ix, iy, 2*iz+1) = im;
                }
            }
        }

        auto t2d1 = hrc::now();
        if(irank==0){
            double dt2d = duration(t2d1 - t2d0).count();
            std::cerr<<"Exercise4/7: 2D slab transforms took "<< dt2d <<" s.\n";
        }
    }

    // (B) Transpose (x,y)->(y,x) with MPI - out-of-place
    {
        // Copy from raw_slab => transpose_in
        ptrdiff_t local_size_pre = local_n0 * nGrid * k_nz;
        {
            ptrdiff_t idx = 0;
            for(int ix=local_0_start; ix<local_0_start+local_n0; ix++){
                for(int iy=0; iy<nGrid; iy++){
                    for(int iz=0; iz<k_nz; iz++){
                        float re = raw_slab(ix, iy, 2*iz  );
                        float im = raw_slab(ix, iy, 2*iz+1);
                        transpose_in[idx][0] = re;
                        transpose_in[idx][1] = im;
                        idx++;
                    }
                }
            }
        }

        // do the transpose
        fftwf_execute(plan_transpose);

        // Now transpose_out has (y,x,z).
        // copy back to raw_slab for demonstration
        ptrdiff_t local_size_post = local_n0_t * nGrid * k_nz;
        {
            ptrdiff_t idx=0;
            for(int iy=local_0_start_t; iy<local_0_start_t+local_n0_t; iy++){
                for(int ix=0; ix<nGrid; ix++){
                    for(int iz=0; iz<k_nz; iz++){
                        float re = transpose_out[idx][0];
                        float im = transpose_out[idx][1];
                        raw_slab(iy, ix, 2*iz  ) = re;
                        raw_slab(iy, ix, 2*iz+1) = im;
                        idx++;
                    }
                }
            }
        }
    }

    // (C) 1D pencils: c2c along x
    {
        std::vector<std::complex<float>> pencil_in(nGrid);
        std::vector<std::complex<float>> pencil_out(nGrid);

        auto t1d0 = hrc::now();

        for(int iy = local_0_start_t; iy<local_0_start_t + local_n0_t; iy++){
            for(int iz=0; iz<k_nz; iz++){
                // copy in
                for(int ix=0; ix<nGrid; ix++){
                    float re = raw_slab(iy, ix, 2*iz);
                    float im = raw_slab(iy, ix, 2*iz+1);
                    pencil_in[ix] = std::complex<float>(re, im);
                }
                // new-array execute for plan_1d
                fftwf_execute_dft(
                    plan_1d,
                    reinterpret_cast<fftwf_complex*>(pencil_in.data()),
                    reinterpret_cast<fftwf_complex*>(pencil_out.data())
                );
                // copy out
                for(int ix=0; ix<nGrid; ix++){
                    float re = pencil_out[ix].real();
                    float im = pencil_out[ix].imag();
                    raw_slab(iy, ix, 2*iz  ) = re;
                    raw_slab(iy, ix, 2*iz+1) = im;
                }
            }
        }

        auto t1d1 = hrc::now();
        if(irank==0){
            double dt1d = duration(t1d1 - t1d0).count();
            std::cerr<<"Exercise6/7: 1D pencil transforms took "<< dt1d <<" s.\n";
        }
    }


    {
        // Use the existing AssignmentWindow from aweights.h
        AssignmentWindow wWin(nGrid, iOrder);  

        Array<double,1> ak(nBins), pk(nBins);
        Array<long,1>   nk(nBins);
        ak=0.; pk=0.; nk=0;

        for(int iy=local_0_start_t; iy<local_0_start_t+local_n0_t; iy++){
            for(int ix=0; ix<nGrid; ix++){
                for(int iz=0; iz<k_nz; iz++){
                    float re = raw_slab(iy, ix, 2*iz  );
                    float im = raw_slab(iy, ix, 2*iz+1);
                    std::complex<float> val(re, im);

                    auto fold = [&](int k){return (k<=nGrid/2) ? k : (k-nGrid);};

                    int ky = fold(iy);
                    int kx = fold(ix);
                    int kz = iz;
                    float kk = std::sqrt(kx*(float)kx + ky*(float)ky + kz*(float)kz);

                    // multiply by window deconvolution
                    val *= wWin[std::abs(kx)] * wWin[std::abs(ky)] * wWin[kz];

                    int ibin;
                    if(bLog){
                        if(kk<=0.f) continue;
                        ibin = int(std::log(kk)/std::log(float(iNyquist))*nBins);
                    } else {
                        ibin = int( kk/ float(iNyquist) * nBins );
                    }
                    if(ibin>=0 && ibin<nBins){
                        ak(ibin) += kk;
                        pk(ibin) += std::norm(val);
                        nk(ibin)++;
                    }
                }
            }
        }

        // reduce across ranks
        Array<double,1> sum_ak(nBins), sum_pk(nBins);
        Array<long,1>   sum_nk(nBins);
        sum_ak=0.; sum_pk=0.; sum_nk=0;

        MPI_Reduce(ak.data(), sum_ak.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(pk.data(), sum_pk.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(nk.data(), sum_nk.data(), nBins, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD);

        if(irank==0){
            for(int i=1; i<nBins; i++){
                if(sum_nk(i)>0){
                    double meanK = sum_ak(i)/sum_nk(i);
                    double meanP = sum_pk(i)/sum_nk(i);
                    std::cout << meanK <<"  "<< meanP <<"  "<< sum_nk(i) << "\n";
                }
            }
        }
    }

    // Cleanup
    fftwf_destroy_plan(plan_2d);
    fftwf_destroy_plan(plan_transpose);
    fftwf_destroy_plan(plan_1d);

    fftwf_free(transpose_in);
    fftwf_free(transpose_out);

    MPI_Finalize();
    return 0;
}
