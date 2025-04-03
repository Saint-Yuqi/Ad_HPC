CXX = /opt/homebrew/Cellar/llvm/19.1.7_1/bin/clang++

BLITZ_INCLUDE = -I /opt/homebrew/Cellar/blitz/1.0.2/include
BOOST_INCLUDE = -I /opt/homebrew/Cellar/boost/1.87.0/include
FFTW_INCLUDE  = -I /opt/homebrew/Cellar/fftw/3.3.10_2/include
OMP_INCLUDE   = -I /opt/homebrew/opt/llvm/include

# OpenMP flags (Homebrew LLVM supports it natively)
OMP     = -fopenmp
LDOMP   = -L /opt/homebrew/opt/llvm/lib -lomp
LDFFTW  = -L /opt/homebrew/Cellar/fftw/3.3.10_2/lib -lfftw3f_omp -lfftw3f

CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(FFTW_INCLUDE) $(OMP_INCLUDE) \
           $(OMP) -O3 -std=c++17 -Wno-deprecated-declarations

assign: assign.o tipsy.h Makefile
	$(CXX) -o assign assign.o $(LDOMP) $(LDFFTW)

assign.o: assign.cpp tipsy.h aweights.h Makefile
	$(CXX) -c $(CXXFLAGS) assign.cpp

clean:
	rm -f assign assign.o
