
CXX       = g++ -pipe
CXX_FLAGS = -mtune=native -march=native -m64 -O3 -fPIC -fopenmp
CXX_INC   = -I/usr/local/include -I.


all: convergence

convergence:
	$(CXX) -o well_convergence.exe lab1.cpp $(CXX_FLAGS) $(CXX_INC) -DWELL_POTENTIAL -DCOSSINBASE
	$(CXX) -o fermi_convergence.exe lab1.cpp $(CXX_FLAGS) $(CXX_INC) -DWOODS_SAXON -DCOSSINBASE

wavefunction:
	$(CXX) -o well_wavefunctions.exe wavefunctions.cpp $(CXX_FLAGS) $(CXX_INC) -DWELL_POTENTIAL -DCOSSINBASE
	$(CXX) -o fermi_wavefunctions.exe wavefunctions.cpp $(CXX_FLAGS) $(CXX_INC) -DWOODS_SAXON -DCOSSINBASE

test:
	$(CXX) -o well_test.exe lab1.cpp $(CXX_FLAGS) $(CXX_INC) -DWELL_POTENTIAL -DCOSSINBASE -DDEBUG
	$(CXX) -o fremi_test.exe lab1.cpp $(CXX_FLAGS) $(CXX_INC) -DWOODS_SAXON -DCOSSINBASE -DDEBUG


clean:
	@rm *.exe
