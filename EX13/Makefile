# You might need to specify the include paths for blitz and/or boost
BLITZ_INCLUDE=-I /capstor/store/cscs/uzh/uzh8/blitz/include/
FFTW_INCLUDE=$(shell pkg-config --cflags fftw3f)
MPICH_INCLUDE=$(shell pkg-config --cflags mpich)
#BOOST_INCLUDE=-I $(BOOST_ROOT)/include
BOOST_INCLUDE = -I /capstor/scratch/cscs/yyuqi/spack/opt/spack/linux-sles15-neoverse_v2/gcc-13.3.0/boost-1.86.0-cjrbpcadr24szno4p6uqyssu4n3luu6j/include

CUDA_INCLUDE=-I $(CUDA_HOME)/include
CUDA_LIB = -L/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/cuda-12.6.0-s7tbmug4jjlhgnp42s5vimkqs2qywpvw/targets/sbsa-linux/lib
CUDA_RT_LIB = -L/user-environment/linux-sles15-neoverse_v2/gcc-13.3.0/cuda-12.6.0-s7tbmug4jjlhgnp42s5vimkqs2qywpvw/lib64

LDFFTW=$(shell pkg-config --libs fftw3f) -lfftw3f_threads -lfftw3f_mpi
LDMPICH=$(shell pkg-config --libs mpich)

CXX	= g++
OMP = -fopenmp
LDOMP = -fopenmp

INCLUDES = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(MPICH_INCLUDE) $(FFTW_INCLUDE) $(CUDA_INCLUDE)
CPPFLAGS = $(INCLUDES) $(OMP)  -O3 -std=c++17
LDLIBS = $(LDMPICH) $(LDFFTW)

all	: assign fft11 fft11-gpu

assign	: assign.o gpufft.o tipsy.h Makefile
	$(CXX) -o assign assign.o gpufft.o $(LDOMP) $(LDFFTW) $(LDFLAGS) $(LDLIBS) $(CUDA_LIB) $(CUDA_RT_LIB) -lcufft -lcudart

assign.o	: assign.cpp tipsy.h aweights.h Makefile

fft11	: fft11.o
	$(CXX) -o fft11 fft11.o $(LDFFTW) $(LDFLAGS) $(LDLIBS)

fft11.o	: fft11.cpp

fft11-gpu	: fft11-gpu.o tipsy.h gpufft.o Makefile
	$(CXX) -o fft11-gpu fft11-gpu.o gpufft.o $(LDFFTW) $(LDFLAGS) $(LDLIBS) $(CUDA_LIB) $(CUDA_RT_LIB) -lcufft -lcudart

fft11-gpu.o	: fft11-gpu.cu
	nvcc -Xcompiler -fopenmp -c -o fft11-gpu.o fft11-gpu.cu $(INCLUDES) -O3 -std=c++17


gpufft.o : gpufft.cu gpufft.h
	nvcc -c -o gpufft.o gpufft.cu $(INCLUDES) -O3 -std=c++17


clean:
	rm -f assign assign.o fft11 fft11.o fft11-gpu fft11-gpu.o gpufft.o
