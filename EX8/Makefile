# Use system GNU compiler (via module cpeGNU)
CXX = mpic++

# Paths for user-installed Blitz++ and Boost
BLITZ_INCLUDE = -I$(HOME)/ESC412/tools/include
BOOST_INCLUDE = -I$(HOME)/ESC412/tools/include

# Paths from loaded cray-fftw module
FFTW_INCLUDE = -I$(FFTW_INC)
FFTW_LIB     = -L$(FFTW_DIR)

# OpenMP flags for GCC
OMP     = -fopenmp
LDOMP   = -fopenmp

# Link FFTW float (single-precision) with OpenMP
LDFFTW  = $(FFTW_LIB) -lfftw3f_threads -lfftw3f -lfftw3f_omp -lfftw3f_mpi

# Final compilation flags
CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(FFTW_INCLUDE) $(OMP) \
           -O3 -std=c++17 -Wno-deprecated-declarations
assign: assign.o tipsy.h Makefile
	$(CXX) -o assign assign.o $(LDOMP) $(LDFFTW)

assign.o: assign.cpp tipsy.h aweights.h Makefile
	$(CXX) -c $(CXXFLAGS) assign.cpp

clean:
	rm -f assign assign.o
