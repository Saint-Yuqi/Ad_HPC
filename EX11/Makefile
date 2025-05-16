################################################################################
# Example Makefile: builds two executables:
#   1) assign_cpu  (CPU-only 2D FFT via FFTW)
#   2) assign_gpu  (GPU 2D FFT via cuFFT)
#
# Also includes your existing fft11-gpu target for reference.
################################################################################

# --- Load your environment (example) ---
# module load daint-gpu cudatoolkit cray-fftw Boost
#
# You might need to specify the include paths for blitz and/or boost

BLITZ_INCLUDE = -I /capstor/store/cscs/uzh/uzh8/blitz/include/
BOOST_INCLUDE = -I /capstor/scratch/cscs/yyuqi/spack/opt/spack/linux-sles15-neoverse_v2/gcc-13.3.0/boost-1.86.0-cjrbpcadr24szno4p6uqyssu4n3luu6j/include
FFTW_INCLUDE  = $(shell pkg-config --cflags fftw3f)
MPICH_INCLUDE = $(shell pkg-config --cflags mpich)
CUDA_INCLUDE  = -I $(CUDA_HOME)/include

LDFFTW  = $(shell pkg-config --libs fftw3f) -lfftw3f_threads -lfftw3f_mpi
LDMPICH = $(shell pkg-config --libs mpich)

CXX     = g++
OMP     = -fopenmp
LDOMP   = -fopenmp

# Collect all includes in one variable
INCLUDES = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) \
           $(MPICH_INCLUDE) $(FFTW_INCLUDE) $(CUDA_INCLUDE)

# Common compiler flags
CPPFLAGS = $(INCLUDES) $(OMP) -O3 -std=c++17
# Libraries to link with
LDLIBS   = $(LDMPICH) $(LDFFTW)

################################################################################
# Default target: build everything
################################################################################
all: assign_cpu assign_gpu fft11-gpu
# You can also add "fft11" or other CPU versions if needed.

################################################################################
# 1) assign_cpu: CPU version using FFTW for 2D FFT
################################################################################
assign_cpu : assign_cpu.o tipsy.h aweights.h
	$(CXX) -o assign_cpu assign_cpu.o $(LDOMP) $(LDFLAGS) $(LDLIBS)

assign_cpu.o : assign_cpu.cpp tipsy.h aweights.h Makefile
	$(CXX) -c assign_cpu.cpp $(CPPFLAGS) -o assign_cpu.o

################################################################################
# 2) assign_gpu: GPU version using cuFFT for 2D FFT
################################################################################
assign_gpu : assign_gpu.o tipsy.h aweights.h
	nvcc -Xcompiler -fopenmp -o assign_gpu assign_gpu.o \
	     $(INCLUDES) $(LDFFTW) $(LDFLAGS) $(LDLIBS) -lcufft -lcudart

assign_gpu.o : assign_gpu.cpp tipsy.h aweights.h Makefile
	nvcc -Xcompiler -fopenmp -c -o assign_gpu.o assign_gpu.cpp \
	     $(INCLUDES) -O3 -std=c++17

################################################################################
# 3) fft11-gpu: (from your reference) GPU example
################################################################################
fft11-gpu : fft11-gpu.o tipsy.h Makefile
	nvcc -Xcompiler -fopenmp -o fft11-gpu fft11-gpu.o \
	     $(INCLUDES) $(LDFFTW) $(LDFLAGS) $(LDLIBS) -lcufft

fft11-gpu.o : fft11-gpu.cu
	nvcc -Xcompiler -fopenmp -c -o fft11-gpu.o fft11-gpu.cu \
	     $(INCLUDES) -O3 -std=c++17

################################################################################
# Clean
################################################################################
clean:
	rm -f assign_cpu assign_cpu.o \
	      assign_gpu assign_gpu.o \
	      fft11-gpu fft11-gpu.o
