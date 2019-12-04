/*============================================================================
Name        : Basic_Ising_Kawasaki.cpp
Author      : Moritz Hoferer
Version     : 2.0
Copyright   :
Description :Basic Ising model with conserved dynamics
			 local Kawasaki algorithm
			 for more information see basic_ising_main.cpp
============================================================================*/

#include <cmath>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>    // std::sort
#include <functional>
#include <deque>


class Lattice{
private:

	/*
	boolean (True VS False)
	periodic_boundary_conditions_i: 	periodic boundary conditions in vertical direction	(1)
										open boundary conditions in vertical direction		(0)
	periodic_boundary_conditions_j: 	periodic boundary conditions in horizontal direction	(1)
										open boundary conditions in horizontal direction		(0)
	record_Density:						average occupation of every site is record			(1)
										average occupation of every site is NOT recorded	(0)
	record_total_Energy:				total Energy and the number of different bounds ( A-A , B-B , A-B ) is recorded		(1)
										total Energy and the number of different bounds ( A-A , B-B , A-B ) is NOT recorded	(0)
	external_field_flag:				simulation with an external field		(1)
										simulation without an external filed	(0)
	tester:								used for control if to neighboring particles are of the same type
	*/
	bool  periodic_bc_i , periodic_bc_j , GlauberProb, record_total_Energy, ext_field_flag ;

	/*
	integer
	size_i:			number of sites in the vertical direction
	size_j:			number of sites in the horizontal direction
	position:		to determine positions in random loops
	position2:
	seed:			to make simulations reproducible
	state1:			state of the first selected particle
	state2: 		state of the second selected particle
	select_i        ( 0 , size_i -1)
	select_j		( 0 , size_j -1)
	select2_i
	select2_j
	select_rel:		relative position to [select_i][select_j]
					1 = above		2 = on the right
					3 = below		4 = on the left
	number_A		number of particles of type A (lipid 1)
	number_B		number of particles of type B (lipid 2)
	number_O		number of particles of type O (obstacles)
	*/
	int size_i, size_j , position , position2 , seed , state1, state2 , select_i , select_j , select2_i , select2_j ,  select_rel , number_A , number_B, number_O, max_dist;

	/*
	unsigned long long integer
	step_count:					number of executed simulation steps
	max_step_count:				after that number of executed simulation the simulation is finished
	burn_in_count:				recording functions start if the burn in count is hit
	executed_movement_count:	count the number of simulations step in which a movement of particles is executed
	*/
	unsigned long long int step_count, max_step_count, burn_in_count , diff_count ,  executed_movement_count;

	/*
	double
	kB:						Boltzmann constant
	temperature:				temperature of the system
	beta:					1/(kB * temperature)
	part_A:					part of particles of type A (movable)
	part_B:					part of particles of type B (movable)
	part_O:					part of particles of type 0 (movable)
	E_AA:					energy of an A-A nearest neighbor bound
	E_BB:					energy of a  B-B nearest neighbor bound
	E_OO					energy of a  O-O nearest neighbor bound
	E_AB:					energy of an A-B nearest neighbor bound
	E_AO:					energy of a  A-O nearest neighbor bound
	E_BO					energy of a  B-O nearest neighbor bound
	bound_E					energy of one single bound
	delta_energy			energy difference caused by the switching of two neighboring particles
	ext_field_max_value:	peak value of the external field
	ext_field_decay			decay of the field from one the the next site (exp{x})
	*/
	double  kB , temperature , beta , part_A , part_B , part_O ,  E_AA , E_BB ,E_OO , E_AB , E_AO , E_BO , bound_E , delta_E ,  ext_field_max_value , ext_field_decay ;


	/*
	ARRAYS/MATRICES
	state:					2d array of size : size_i * size_j and every position represents a position in the system
							if state[i][j] 	= +1 site is occupied by a particle of type A,
											= -1 site is occupied by a particle of type B,
											=  O site is occupied by a particle of type O
	m						2d array: used for the Hoshen-Kopelman-algorithm
	density					2d array: probability density of cluster w/o O (first column) and w/ O (second)
	external_field:			2d array: contains the field strength that is appealing to each single site
	total_energy:			1d array: contains the total energy of the system after every diff_count simulation steps
	number_of_even:			1d array: contains the total number of even bounds after every diff_count simulations steps
	number_of_odd:			1d array: contains the total number of odd bounds after every diff_count simulations steps

	*/
	int state[512][512];
	int m[512][512];
	unsigned long long int density[262155][2];
	// double external_field[max_size][max_size];
	std::deque<int> urn;
	std::deque<int> typeA_i ;
	std::deque<int> typeA_j ;
	std::deque<int> typeB_i ;
	std::deque<int> typeB_j ;
	std::deque<int> obs_pos_i;
	std::deque<int> obs_pos_j;
	std::deque<int> abs_obs_pos_i;
	std::deque<int> abs_obs_pos_j;
	std::deque<double> clsize;
	std::deque<double> clsize_obs;
	std::deque<double> total_energy;
	std::deque<unsigned long int> N;
	std::deque<unsigned long int> number_of_even;
	std::deque<unsigned long int> number_of_odd;
	std::deque<int> nn_i;
	std::deque<int> nn_j;
	std::deque<unsigned long long int> nn_t;

	// set random number generator
	std::mt19937 engine;

	// select_4			: to select the neighboring particles
	// select_doubele	: to select everything else
	std::uniform_int_distribution< > select_4{1,4};
	std::uniform_real_distribution< > select_double{0.,1.};


	// private functions
	// fills the lattice with particles regarding part_A, part_B
	void fill_sys();

	// replace some particles of type A with particles of type_O, that the ratios are fulfilled
	void add_typeO();

	// determine the i relative to an selected i depending on rel
	int rel_i(int pos_i, int pos_rel);

	// determine the j relative to an selected j depending on rel
	int rel_j(int pos_j, int pos_rel);

	// select two particles of different type of different type
	void select_2part();

	// select a pair of two neighboring particles of different type
	void select_pair();

	// determine the energy of a two neighboring particles
	double bound_energy( int pos_i , int pos_j , int pos_rel) ;


	// determine the delta energy that is caused switching the positions of two neighboring particles
	void delta_energy();

	// depending on the determined delta energy switch the positions of the two particles // using Metropolis
	void switching();

	// depending on the determined delta energy switch the positions of the two neighboring particles
	void movement();

	// record the evolution of the energy of the system
	void record_energy();

	// record the evolution of the number of bound including to particles of the same type and including particles of different type
	void record_bounds();

	// record the cluster sizes and track the (absolute) obstacles position
	void record_cluster_sizes();

	// determine proper label
	int classify(int &m);

	// determines the cluster label
	int find_cluster(int &m ) const;

	bool element_of(int &k, std::deque<int> &mnew)const;

	void pbc_to_inf_sys();

	// record the minority lipids surrounding the obstacle
	void nn_scan();

	// writes the used parameters in a file
	void write_parameter(std::ofstream &write) const;

	// writes state of the lattice and time into a file
	void write_state(std::ofstream &write) const;


public:

	// initialize Lattice
	// insert number of the relevant parameters
	Lattice(const double Parameters[21]);

	//destructor
	~Lattice();

	// burn in simulation with non local Kawasaki dynamics
	void burn_in_sim();

	// burn in simulation with local Kawasaki dynamics
	void simulate();

	// simulate and record with Kawasaki dynamics
	void simulate_rec(std::string &identifier);

	// print all input parameters
	void print_parameters() const;

	// print the current state of the system with length and time
	void print_state() const;

	// print the number of particles of every type in the system
	void print_part_count() const;

	// print same information at the end of the simulation
	void print_final_info() const;

	// True if there is at least on obstacle in the system
	bool recording_obs_pos() const;

	// gives record_total_Energy
	bool recording_energy() const;

	// recording density function angle and radius dependent (nn_scan)
	bool recording_density() const;

	// writes the cluster sizes with and without obstacle(s) into a file
	void write_cluster_size(std::ofstream &write) const;

	// writes the probability density of the cluster sizes with and without obstacle(s) into a file
	void write_density(std::ofstream &write) const;

	// writes the obstacle positions into a file
	void write_obs_pos(std::ofstream &write) const;

	// writes the obstacle positions into a file
	void write_abs_obs_pos(std::ofstream &write) const;

	// writes the array with total energies of the system into a file
	void write_energy(std::ofstream &write) const;

	// write the array with the number of bound into a file
	void write_bounds(std::ofstream &write) const;

	//write nn_scan: time_count, pos_i, pos_j
	void write_nn(std::ofstream &write) const;

	};


// used to convert a number (int, double, etc.) into a string
std::string do_to_str(const double &d);

// used to read in the configuration parameters out of the parameters[i].txt-file
void read_in_parameters(const std::string &str, double Params[]);
