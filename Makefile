CXX = g++
CXXFLAGS = -static -flto -march=x86-64  -mtune=native -O2 -std=c++11

all: lattice_membrane

clean:
	rm -f *.o lattice_membrane

lattice_membrane: lattice_membrane_main.o lattice_membrane_core.o
	$(CXX) -o lattice_membrane lattice_membrane_main.o lattice_membrane_core.o $(CXXFLAGS)

lattice_membrane_main.o: lattice_membrane_main.cpp lattice_membrane_core.h timer.h
	$(CXX) -c lattice_membrane_main.cpp $(CXXFLAGS)

lattice_membrane_core.o: lattice_membrane_core.cpp lattice_membrane_core.h
	$(CXX) -c lattice_membrane_core.cpp $(CXXFLAGS)
