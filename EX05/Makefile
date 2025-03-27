CXX = g++

# Blitz, Boost
BLITZ_INCLUDE = -I /opt/homebrew/Cellar/blitz/1.0.2/include
BOOST_INCLUDE = -I /opt/homebrew/Cellar/boost/1.87.0/include

# FFTW single-precision
FFTW_INCLUDE  = -I /opt/homebrew/Cellar/fftw/3.3.10_2/include
LDFFTW        = -L /opt/homebrew/Cellar/fftw/3.3.10_2/lib -lfftw3f -lm

# MacOS OpenMP
OMP   = -Xpreprocessor -fopenmp
LDOMP = -L$(shell brew --prefix libomp)/lib -lomp

# Compiler flags
CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(FFTW_INCLUDE) \
           $(OMP) -Wno-deprecated-declarations -O3 -std=c++17

# Libraries to link
LDLIBS = $(LDFFTW) $(LDOMP)

all: assign test_weights

assign: assign_4.o tipsy.h Makefile
	$(CXX) -o assign assign_4.o $(LDLIBS) $(LDFLAGS)

assign_4.o: assign_4.cpp tipsy.h aweights.h Makefile
	$(CXX) $(CXXFLAGS) -c assign_4.cpp

test_weights: test_weights.cpp weights.h
	$(CXX) -std=c++11 -O3 -o test_weights test_weights.cpp

clean:
	rm -f assign assign_4.o test_weights density.dat
