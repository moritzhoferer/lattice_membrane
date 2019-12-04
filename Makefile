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


# if [ ! -e lmu/ ]; then
# 	mkdir -p lmu/
# fi

# if [ ! -e unimi/ ]; then
# 	mkdir -p unimi/
# fi

# # compiles the program for the asc/lmu cluster
# src/basic_ising_main.cpp src/basic_ising_core.cpp -o lmu/basic_ising -O2 -std=c++11
# # old
# #g++-5 src/basic_ising_main.cpp src/basic_ising_core.cpp -o lmu/basic_ising -static -flto -march=native  -mtune=native -O2 -std=c++11

# # compiles the program for the unimi cluster
# g++-5 src/basic_ising_main.cpp src/basic_ising_core.cpp -o unimi/basic_ising 