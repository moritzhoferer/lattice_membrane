/*============================================================================
Name        : Basic_Ising_Kawasaki.cpp
Author      : Moritz Hoferer
Date		   : 20.09.2016
Version     : 2.0
Copyright   :
Description :Basic Ising model with conserved dynamics
			 local Kawasaki algorithm
			 during the burn - in phase (with only type A and B particles) of the simulation a nonlocal Kawasaki algorithm
			 is used to get faster into the equilibrium state. Type O particles (obstacles) are added later (2nd burn in)
			 NEW way of counting simulation steps, that sim. steps are comparable to a real time
			 (1 steps is defined as size_i*size_j sim. steps)
			 New in 2.0:
			 no more screen shots of the whole system
			 new recorded quantities: cluster sizes with and without obstacle, obstacles position in the system with periodic and with infinite boundary conditions

Note		:External field function is missing (until now not needed),
			 also open boundary conditions
			 recording density can be used to record snapshot with some exposure time
============================================================================*/

#include <iostream>
#include <string>
#include <iomanip>

#include "lattice_membrane_core.h"
#include "timer.h"

using namespace std;

int main(int argc, char *argv[]) {
	double Parameters[25];

	// initialization of the timer
	timer control;
	// start the timer
	control.start();

	// Initialization of the output files
	ofstream file0, file1, file2, file3, file4, file5, file6;
	string str, identifier, out0, out1, out2, out3, out4, out5, out6;

	// define name of input file
	str = "Parameters";
	str.append(argv[1]);
	str.append(".txt");

	// read in the control parameter for lattice configuration and particle properties
	read_in_parameters(str, Parameters);

	// construct an identifier
	identifier = "_T_";
	identifier.append( do_to_str(Parameters[5]));
	if(Parameters[6] != 0 && Parameters[7] != 0 && Parameters[8] != 0){
		identifier.append("_E_AA_");
		identifier.append( do_to_str(Parameters[9]));
		identifier.append("_E_BB_");
		identifier.append(do_to_str(Parameters[10]));
		identifier.append("_E_OO_");
		identifier.append(do_to_str(Parameters[11]));
	}

	if(Parameters[6] != 0 && Parameters[7] != 0){
		identifier.append("_E_AB_");
		identifier.append( do_to_str(Parameters[12]));
	}
	if(Parameters[6] != 0 && Parameters[8] != 0){
		identifier.append("_E_AO_");
		identifier.append(do_to_str(Parameters[13]));
	}
	if(Parameters[7] != 0 && Parameters[8] != 0){
		identifier.append("_E_BO_");
		identifier.append(do_to_str(Parameters[14]));
	}

	// initialization of the test system
	Lattice sys(Parameters);

	// print_parameters
	sys.print_parameters();

	// simulate the system w/o record || BRUN IN
	// non local Kawasaki algorithm with A and B
	sys.burn_in_sim();

	// local Kawaski algorithm with A, B and O
	sys.simulate();

	// simulate the system w/ record || SIMULATION
	sys.simulate_rec(identifier);


	// cluster sizes | obstacle position | total energy | bounds || output
	// construction of the names of the output files
	// open the output file
	// write the lattice parameters and particle properties into the output files
	// close the output files

	out0 = "cluster_size" + identifier + ".txt";
	file0.open(out0.c_str());
	sys.write_cluster_size(file0);
	file0.close();

	out3 = "density" + identifier + ".txt";
	file3.open(out3.c_str());
	sys.write_density(file3);
	file3.close();

	if( sys.recording_obs_pos()){
		out4 = "obs_pos" + identifier + ".txt";
		file4.open(out4.c_str());
		sys.write_obs_pos(file4);
		file4.close();

		out5 = "abs_obs_pos" + identifier + ".txt";
		file5.open(out5.c_str());
		sys.write_abs_obs_pos(file5);
		file5.close();
	}

	if( sys.recording_energy() ){
		out1 = "energy" + identifier + ".txt";
		file1.open(out1.c_str());
		sys.write_energy(file1);
		file1.close();

		out2 = "bounds" + identifier + ".txt";
		file2.open(out2.c_str());
		sys.write_bounds(file2);
		file2.close();
	}

	if( sys.recording_density() ){
		out6 = "nn" + identifier + ".txt";
		file6.open(out6.c_str());
		sys.write_nn(file6);
		file6.close();
	}

	// output in the console after finishing the simulation

	// stop the time measurement
	control.stop();
	sys.print_final_info();
	control.check();

	return 0;
}
