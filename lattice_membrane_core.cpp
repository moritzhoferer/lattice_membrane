/*============================================================================
Name        : Basic_Ising_Kawasaki.cpp
Author      : Moritz Hoferer
Version     : 2.0
Copyright   :
Description :Basic Ising model with conserved dynamics
			 local Kawasaki algorithm
			 for more information see basic_ising_main.cpp
============================================================================*/

#include "lattice_membrane_core.h"

using namespace std;


Lattice::Lattice(const double Parameters[25]){

	step_count = 0;
	executed_movement_count = 0;

	total_energy.clear();
	number_of_even.clear();
	number_of_odd.clear();


	// set all parameters
	// read in from parameter file
	size_i				= int(Parameters[0]);
	size_j				= int(Parameters[1]);
	periodic_bc_i		= bool(Parameters[2]);
	periodic_bc_j		= bool(Parameters[3]);
	kB					= Parameters[4];
	temperature			= Parameters[5];
	beta				= 1 / (kB * temperature);

	part_A				= Parameters[6];
	part_B				= Parameters[7];
	part_O				= Parameters[8];
	E_AA				= Parameters[9];
	E_BB				= Parameters[10];
	E_OO				= Parameters[11];
	E_AB				= Parameters[12];
	E_AO				= Parameters[13];
	E_BO				= Parameters[14];

	ext_field_flag		= bool(Parameters[15]);
	ext_field_max_value	= Parameters[16];
	ext_field_decay		= Parameters[17];

	seed 				= int(Parameters[18]);
	GlauberProb			= bool(Parameters[19]);
	max_step_count		= pow(10,int(Parameters[20]))* size_i * size_j;
	burn_in_count		= pow(10,int(Parameters[21]))* size_i * size_j;

	if(Parameters[20] - Parameters[22] <= 7){
		diff_count = pow(10,int(Parameters[22]))* size_i * size_j;
	}
	else{
		diff_count = pow(10,int(Parameters[20]-7)) * size_i * size_j;
	}

	max_dist		  =  Parameters[23];
	record_total_Energy	= bool(Parameters[24]);

	// USE RANDOM SEED FROM PARAMETER FILE
	engine.seed(seed);

	for(long int i=0; i<size_i*size_j+1; i++){
		for(int j=0; j<2; j++){
			density[i][j] = 0;
		}
	}

	// fill the lattice with "exact rounded" spatial initial density
	fill_sys();

}


Lattice::~Lattice(){

	// clear all arrays
	// NOT NESSASSARY

}


void Lattice::burn_in_sim(){

	cout << "Burn-In simulation with non local Kawasaki algorithm was started." << endl;

	// store the positions of all particles of type A and B in vectors
	for(int i = 0; i < size_i; i++){
		for(int j = 0; j < size_j; j++ ){
			if( state[i][j] == 1){
				typeA_i.push_back(i);
				typeA_j.push_back(j); }
			else if( state[i][j] == -1){
				typeB_i.push_back(i);
				typeB_j.push_back(j);	}
		}
	}

	step_count = 0;

	while(step_count < burn_in_count){

		select_2part();

		delta_E = 0;

		for( int i = 1 ; i <= 4 ; i++){
			delta_E -= bound_energy(select_i, select_j, i ) + bound_energy(select2_i , select2_j , i);
		}

		// switch the particle's positions
		state[select_i][select_j] = -1;
		state[select2_i][select2_j] = 1;

		for( int i = 1 ; i <= 4 ; i++){
			delta_E += bound_energy(select_i, select_j, i ) + bound_energy(select2_i , select2_j , i);
		}

		// undo the switch
		state[select_i][select_j] = 1;
		state[select2_i][select2_j] = -1;

		switching();

		step_count++;

	}

	if( number_O != 0 ){ add_typeO(); }

}


void Lattice::simulate(){

	step_count = 0;

	cout << "Simulation with the local Kawasaki algorithm without recording was started." << endl;

	while(step_count < 0.1 * burn_in_count){

		// select two neighboring particles of different type
		select_pair();


		if( state[select_i][select_j] != state[rel_i(select_i,select_rel)][rel_j(select_j,select_rel)]){

			state1 = state[select_i][select_j];
			state2 = state[rel_i(select_i, select_rel)][rel_j(select_j,select_rel)];

			// determine the delta energy that would result by switching the positions of the two selected particles
			delta_energy();

			// maybe the switching is executed
			movement();

		}

		step_count++;

	}
}


void Lattice::simulate_rec(std::string &identifier){

	step_count = 0;
	executed_movement_count = 0 ;

	cout << "Simulation with recording was started." << endl;

	while(step_count <= max_step_count){

		// select two neighboring particles of different type
		select_pair();

		if( state[select_i][select_j] != state[rel_i(select_i,select_rel)][rel_j(select_j,select_rel)]){

			state1 = state[select_i][select_j];
			state2 = state[rel_i(select_i, select_rel)][rel_j(select_j,select_rel)];

			// determine the delta energy that would result by switching the positions of the two selected particles
			delta_energy();

			// maybe the switching is executed
			movement();

		}


		// record the total energy and the number of even/odd bounds and the cluster sizes
		if( step_count % diff_count == 0 ){

			if(record_total_Energy ){
				record_energy();
				record_bounds();
			}

			record_cluster_sizes();

			if(max_dist>0 && part_O == 1 ){
				nn_scan();
			}
		}

		// recording 101 snap shots out of every simulation
		if( step_count % (max_step_count/100) == 0){

			// record state
			ofstream file6;
			string out6;

			out6 = "state_";
			out6.append( do_to_str( step_count / (max_step_count / 100) ) );
			out6.append(identifier);
			out6.append(".txt");

			file6.open(out6.c_str());

			write_state(file6);

			file6.close();

			cout << double(step_count)/double(max_step_count)*100 << "% of the simulation are executed." << endl;

		}

		step_count++;

	}

	if(number_O == 1){
		pbc_to_inf_sys();
	}

}


// fill the system with the "exact" initial spatial density rho_0
void Lattice::fill_sys(){

	// determine the number of positions that are occupied by A, B or O, respectively
	// if part_X of any type is 0, this is a hard condition, that no particle of this type is in the system (prevent rounding errors)
	if( part_O != 0){
		number_A = trunc(size_i * size_j * ( (part_A ) / ( part_A + part_B + part_O ) ));
		number_B = trunc(size_i * size_j *( part_B / ( part_A + part_B + part_O ) ));
		number_O = size_i * size_j - number_A - number_B;
	}
	else{
		number_A = trunc(size_i * size_j * ( (part_A + part_O) / ( part_A + part_B + part_O ) ));
		number_B = size_i * size_j - number_A;
		number_O = 0;
	}

	// fill the whole lattice with 0
	for( int ij = 0 ; ij < size_i * size_j ; ij++ ){ state[ ij%size_i ][ ij/size_i ] = 0; }

	// determine the initial size of the urn
	urn.resize( size_i * size_j , 0 );

	// fill urn
	for( int ij = 0 ; ij < size_i * size_j ; ij++){urn[ij] = ij ; }


	// select positions that are occupied by a particle of type A
	for( int ij = 0 ; ij < number_A + number_O ; ij++){

		position = int(urn.size() * select_double(engine) );

		state[ urn[position] % size_i][ urn[position] / size_i ] = 1;

		urn.erase(urn.begin() + position);
	}

	// select positions that are occupied by a particle of type B
	for( int ij = 0 ; ij < number_B ; ij++){

		position = int(urn.size() * select_double(engine));

		state[ urn[position] % size_i][ urn[position] / size_i ] = -1;

		urn.erase(urn.begin() + position);
	}
}


void Lattice::add_typeO(){

	for(int i = 0 ; i < number_O ; i++){

		position = int( typeA_i.size() * select_double(engine) );

		state[typeA_i[position]][typeA_j[position]] = 0;

		typeA_i.erase(typeA_i.begin() + position);
		typeA_j.erase(typeA_j.begin() + position);

	}
}


// pos_rel: 1 (above), 2 (on the right), 3(below), 4 (on the left)
int Lattice::rel_i(int pos_i, int pos_rel){
	return (pos_i + size_i - 1 * (pos_rel == 1) + 1 * (pos_rel == 3))%size_i ;
}

int Lattice::rel_j(int pos_j, int pos_rel){
	return (pos_j + size_j - 1 * (pos_rel == 4) + 1 * (pos_rel == 2))%size_j ;
}


void Lattice::select_2part(){

	// select a particle of type A
	position = int(select_double(engine) * typeA_i.size());

	select_i = typeA_i[position];
	select_j = typeA_j[position];


	// select a particle of type B
	position2 = int(select_double(engine) * typeB_i.size());

	select2_i = typeB_i[position2];
	select2_j = typeB_j[position2];

}


void Lattice::select_pair(){

	// select one random particle in the lattice
	position =  int(size_i * size_j * select_double(engine));
	select_i = position / size_j;
	select_j = position % size_j;

	// select a neighboring particle
	// (1) above		(2) on the right
	// (3) below		(3) on the left
	select_rel = select_4(engine);

	// Think about some array that contains all possible connections and pick one random element from this array

}


double Lattice::bound_energy( int pos_i ,  int pos_j ,  int pos_rel) {


	// check if both particles have the same type
	if(state[pos_i][pos_j] == state[rel_i(pos_i , pos_rel)][rel_j(pos_j , pos_rel)]){
		 // A-A bound
		 if( state[pos_i][pos_j] ==  1 ){ bound_E = E_AA; }
		 // B-B bound
		 else if( state[pos_i][pos_j] == -1 ){ bound_E = E_BB; }
		 // O-O bound
		 else if( state[pos_i][pos_j] ==  0 ){ bound_E = E_OO; }
	}
	// both particles have not the same type
	else{
		if( state[pos_i][pos_j] == 1){
			// A-B bound
			if( state[rel_i(pos_i , pos_rel)][rel_j(pos_j , pos_rel)] == -1 ){ bound_E = E_AB; }
			// A-O bound
			else{ bound_E = E_AO; }
		}
		else if( state[pos_i][pos_j] == -1){
			// B-A bound
			if( state[rel_i(pos_i , pos_rel)][rel_j(pos_j , pos_rel)] == 1 ){ bound_E = E_AB; }
			// B-O bound
			else{ bound_E = E_BO; }
		}
		else if( state[pos_i][pos_j] == 0 ){
			// O-A bound
			if (state[rel_i(pos_i , pos_rel)][rel_j(pos_j , pos_rel)] == 1){ bound_E = E_AO; }
			// O-B bound
			else{ bound_E = E_BO; }
		}
	}

	return bound_E;

}


void Lattice::delta_energy(){

	// set delta energy to 0
	delta_E = 0;

	for(int rel = 1 ; rel <= 4 ; rel++){
		// delta_E = E_after - E_before
		// before switching
		delta_E -=   bound_energy( select_i , select_j , rel)
					+bound_energy( rel_i( select_i , select_rel ) , rel_j( select_j , select_rel ) , rel );

	}

	// switch particle positions;
	state[select_i][select_j] = state2;
	state[rel_i( select_i , select_rel )][rel_j( select_j , select_rel )] = state1;

	for(int rel = 1 ; rel <= 4 ; rel++){
		// after switching
		delta_E +=   bound_energy( select_i , select_j , rel)
					+bound_energy( rel_i ( select_i , select_rel ),  rel_j( select_j , select_rel ) ,  rel );
	}

	// switch back particle positions;
	state[select_i][select_j] = state1;
	state[rel_i( select_i , select_rel) ][rel_j( select_j , select_rel)] = state2;

}


void Lattice::switching(){

	if ( select_double(engine) < fmin( 1., exp(- beta * delta_E) ) ){

		// switch particle positions;
		state[select_i][select_j] = -1;
		state[select2_i][select2_j] = 1;

		// add the new positions to the arrays
		typeA_i.push_back(typeB_i[position2]) ; typeA_j.push_back(typeB_j[position2]);
		typeB_i.push_back(typeA_i[position]) ; typeB_j.push_back(typeA_j[position]);

		// remove the old position from the arrays
		typeA_i.erase(typeA_i.begin()+ position) ;typeA_j.erase(typeA_j.begin() + position);
		typeB_i.erase(typeB_i.begin()+ position2) ; typeB_j.erase(typeB_j.begin() + position2);

	}

}


void Lattice::movement(){

	if(GlauberProb){
		if ( select_double(engine) <= exp( - beta * delta_E) / ( 1 + exp( - beta * delta_E) ) ){
			// switch particle positions;
			state[select_i][select_j] = state2;
			state[rel_i(select_i , select_rel )][rel_j( select_j , select_rel )] = state1;
			executed_movement_count++;
		}
	}
	else{
		if ( select_double(engine) < fmin( 1., exp(- beta * delta_E) ) ){
			// switch particle positions;
			state[select_i][select_j] = state2;
			state[rel_i( select_i , select_rel )][rel_j( select_j , select_rel )] = state1;
			executed_movement_count++;
		}
	}


}


void Lattice::record_energy(){

	double tot_E = 0;

	for(int i = 0 ; i < size_i ; i++){

		for(int j = 0 ; j < size_j ; j++){

			// add up the bound energies of each particles to its right and lower neighbor
			tot_E += bound_energy( i , j , 2 )  + bound_energy( i , j , 3 );

		}
	}

	// append one extra entry at the end of the array
	total_energy.push_back(tot_E);



}


void Lattice::record_bounds(){
	int even = 0;
	int odd = 0;

	for(int i = 0 ; i < size_i ; i++){

		for(int j = 0 ; j < size_j ; j++){
			// check if there are two particles of the same type or different
			if( state[i][j] == state[rel_i(i,2)][rel_j(j,2)]){ even++; }
			else{ odd++;}
			// check if there are two particles of the same type or different
			if(state[i][j] == state[rel_i(i,3)][rel_j(j,3)]){ even++; }
			else{ odd++;}
		}

	}

	// append the results to the vector
	number_of_even.push_back(even);
	number_of_odd.push_back(odd);


}


void Lattice::record_cluster_sizes(){

	// set all values in m to 0
	for(int i = 0; i < size_i; i++){
		for(int j = 0; j < size_j; j++){
			m[i][j] = 0;
		}
	}

	N.clear();
	N.push_back(0);
	int k = 0;


	deque<int> mnew;
	int tmp;

	// execute the Hoshen-Kopelmann-Algorithm
	for(int i = 0; i < size_i; i++){
		for(int j = 0; j < size_j; j++){

			// if particles is NOT of type B
			if(state[i][j] != -1){

				// 0bstacle tracker
				if( state[i][j] == 0 ){
					obs_pos_i.push_back(i);
					obs_pos_j.push_back(j);
				}

				// if any of the four neighboring particles is already classify as type A or O
				if( m[rel_i(i,1)][rel_j(j,1)]!= 0 || m[rel_i(i,2)][rel_j(j,2)]!= 0 || m[rel_i(i,3)][rel_j(j,3)]!= 0 || m[rel_i(i,4)][rel_j(j,4)]!= 0){

					mnew.clear();

					for(int rel = 1; rel <= 4; rel++){

						if(m[rel_i(i,rel)][rel_j(j,rel)] != 0 && !(element_of(m[rel_i(i,rel)][rel_j(j,rel)],mnew)) ){
							tmp = classify( m[rel_i(i,rel)][rel_j(j,rel)] );
							if(!(element_of(tmp,mnew))){
								mnew.push_back(tmp);
							}
						}

					}

					// sort the element in ascending order
					sort(mnew.begin(),mnew.end());

					// m[i][j] is the smallest element out of mnew
					m[i][j] = mnew[0];

					// rename the other clusters and sum all new connected cluster up
					for(unsigned int u = 1 ; u < mnew.size(); u++){
						N[mnew[0]] += N[mnew[u]];
						N[mnew[u]] = -mnew[0];
					}

					N[mnew[0]] += 1;
				}
				// none of the surrounding particles is of type A or O and already classified
				else{
					k++;
					m[i][j] = k;
					N.push_back(1);
				}
			}
		}
	}

	// determine if the clusters contain an obstacle or not and average the size
	double clsize_tmp = 0, clsize_obs_tmp = 0;
	int n_cl = 0, n_cl_obs = 0;

	// search the number of clusters containing at least one obstacle and the sizes
	for(unsigned int i = obs_pos_i.size() - number_O ; i < obs_pos_i.size(); i++){

		int cl_no = find_cluster(m[obs_pos_i[i]][obs_pos_j[i]]);

		// used to prohibit over counting. every cluster containing more than one obstacle is just counted once.
		if(N[cl_no] != 0){
			clsize_obs_tmp += N[cl_no];
			density[N[cl_no]][1]++;
			density[0][1]++;
			N[cl_no] = 0;
			n_cl_obs++;
		}
	}

	// if there is at least one cluster with an obstacle divided the sum of the sizes by the number of cluster containing at least one obstacle
	if(n_cl_obs > 0){
		clsize_obs_tmp /= double(n_cl_obs);
	}

	// sum up the clusters without any obstacle
	for(unsigned int i = 0; i < N.size(); i++){
		if(N[i] > 0){
			clsize_tmp += N[i];
			density[N[i]][0]++;
			density[0][0]++;
			n_cl++;
		}
	}

	// if there is at least one cluster without any obstacle divided the sum of the sizes by the number of clusters containing no obstacle
	if(n_cl > 0){
		clsize_tmp /= double(n_cl);
	}

	// add the sizes to arrays, with the time series
	clsize.push_back(clsize_tmp);
	clsize_obs.push_back(clsize_obs_tmp);

	// delete m, N , k;
}


int Lattice::classify(int &m){
	int x =  m;
	int y = -N[x];


	if(y > 0){
		x =  y;
		y = -N[y];

		if(y > 0){
			while(y > 0){
				x = y;
				y = -N[y];
			}

			N[m] = -x;
		}
	}

	return x;
}


int Lattice::find_cluster(int &m)const{

	int x = m;
	int y = x;

	if(y > 0){
		while(y > 0){
			x = y;
			y = -N[y];
		}
	}

	return x;

}


bool Lattice::element_of(int &k, std::deque<int> &array)const{

	bool tester = false;
	unsigned int i = 0;

	while( i < array.size() && !(tester) ){
		if(k == array[i]){
			tester = true;
		}
		i++;
	}

	return tester;
}

void Lattice::pbc_to_inf_sys(){
	int sys_i = 0, sys_j = 0, tester[9],min_d ;
	int delta[3] = { -1 , 0 , 1 };


	abs_obs_pos_i.push_back(obs_pos_i[0]);
	abs_obs_pos_j.push_back(obs_pos_j[0]);

	for(unsigned int s = 1; s < obs_pos_i.size(); s++){

		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				tester[i * 3 + j] = pow(obs_pos_i[s] + delta[i] *size_i- obs_pos_i[s-1],2) + pow(obs_pos_j[s] + delta[j] *size_j- obs_pos_j[s-1],2) ;
			}
		}

		// MIN FUNC
		min_d = tester[0];
		for(int i = 1; i < 9; i++){
			if(tester[i]< min_d){min_d = tester[i];}
		}

		if(tester[1 * 3 + 1] != min_d){
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
					if (tester[i*3 + j] == min_d){
						sys_i += delta[i];
						sys_j += delta[j];
					}
				}
			}
		}

		abs_obs_pos_i.push_back(obs_pos_i[s] + sys_i * size_i);
		abs_obs_pos_j.push_back(obs_pos_j[s] + sys_j * size_j);

	}

}

void Lattice::nn_scan(){

	int dist_i, dist_j;

	for(int i=-max_dist; i < max_dist + 1; i++){
		for(int j=-max_dist; j < max_dist + 1; j++){

			if( pow(i,2) + pow(j,2) <= pow(max_dist,2) ){

				if( i + obs_pos_i.back() > -1 && i + obs_pos_i.back() < size_i){
					dist_i = i + obs_pos_i.back();
				}
				else if( i + obs_pos_i.back() < 0){
					dist_i = i + obs_pos_i.back() + size_i;
				}
				else{
					dist_i = i + obs_pos_i.back() - size_i;
				}

				if( j + obs_pos_j.back() > -1 && j + obs_pos_j.back() < size_j){
					dist_j = j + obs_pos_j.back();
				}
				else if( j + obs_pos_j.back() < 0){
					dist_j = j + obs_pos_j.back() + size_j;
				}
				else{
					dist_j = j + obs_pos_j.back() - size_j;
				}

				if(state[dist_i][dist_j]==1){
					nn_t.push_back(step_count / (size_i * size_j));
					nn_i.push_back( i );
					nn_j.push_back( j );
				}
			}
		}
	}
}


void Lattice::print_parameters() const{
	cout << "Simulation started with the following properties and parameters" << endl;
	cout << "Lattice properties" << endl;
	cout << "Horizontal system size:              " << size_i << endl;
	cout << "Vertical system size:                " << size_j << endl;
	cout << "Periodic BC in horizontal direction: " << periodic_bc_i << endl;
	cout << "Periodic BC in vertical direction:   " << periodic_bc_j << endl;
	cout << "Boltzmann constant:                  " << kB << endl;
	cout << "Temperature:                         " << temperature << endl;
	cout <<  endl;
	cout << "Particles properties" << endl;
	cout << "Part of particles of type A:         " << part_A << endl;
	cout << "Part of particles of type B:         " << part_B << endl;
	cout << "Part of particles of type O:         " << part_O << endl;
	cout << "Energy A-A bound:                    " << E_AA << endl;
	cout << "Energy B-B bound:                    " << E_BB << endl;
	cout << "Energy O-O bound:                    " << E_OO << endl;
	cout << "Energy A-B or B-A bound:             " << E_AB << endl;
	cout << "Energy A-O or O-A bound:             " << E_AO << endl;
	cout << "Energy B-O or O-B bound:             " << E_BO << endl;
	cout << endl;
	cout << "External field" << endl;
	cout << "External field On:                   " << ext_field_flag << endl;
	cout << "Maximum value of the ext field:      " << ext_field_max_value << endl;
	cout << "Decay of the external filed:         " << ext_field_decay << endl;
	cout <<  endl;
	cout << "Simulations properties" << endl;
	cout << "Seed:                                " << seed << endl;
	cout << "Glauber probability used:            " << GlauberProb << endl;
	cout << "Burn-In simulation steps:            " << burn_in_count / (size_i * size_j) << endl;
	cout << "Maximum number of simulation steps:  " << max_step_count / (size_i * size_j) << endl;
	cout << "Steps btw. recordings:               " << diff_count / (size_i * size_j) << endl;
	cout << "# Max. radius (nn_scan):             " << max_dist << endl;
	cout << "Energy/bounds recorded:              " << record_total_Energy << endl;

}


void Lattice::print_state() const{
	cout << "State of the system after " << step_count << ":" << endl;
	for(int i = 0 ; i < size_i ; i++){
		for(int j = 0 ; j < size_j ; j++){
			if(state[i][j] == 1){ cout << " 1\t";}
			else if (state[i][j] == -1 ){cout << "-1\t";}
			else{cout << " 0\t";}
			//cout << state[i][j] << "\t";
		}
		cout << endl;
	}
}


void Lattice::print_part_count() const {
	int a = 0;
	int b = 0;
	int c = 0;

	for(int i = 0; i< size_i; i++){
		for(int j =0 ; j < size_j ; j++){
			if(state[i][j] ==  1){ a++; }
			else if(state[i][j] == -1){ b++; }
			else if(state[i][j] ==  0){ c++; }
		}
	}

	cout << "# of particles A: " << a << endl;
	cout << "# of particles B: " << b << endl;
	cout << "# of particles O: " << c << endl;
}


void Lattice::print_final_info() const{
	cout << "During the simulation in " << setprecision(2) << double(executed_movement_count) / double(max_step_count) * 100. << "% of the steps a movement was executed." << endl;
}


void Lattice::write_parameter(std::ofstream &write) const{
	write << "# Simulation of a system with the following properties and parameters" << endl;
	write << "# Lattice properties" << endl;
	write << "# Horizontal system size:                 " << size_i << endl;
	write << "# Vertical system size:                   " << size_j << endl;
	write << "# Periodic BC in horizontal direction:    " << periodic_bc_i << endl;
	write << "# Periodic BC in vertical direction:      " << periodic_bc_j << endl;
	write << "# Boltzmann constant:                     " << kB << endl;
	write << "# Temperature:                            " << temperature << endl;
	write << "#" << endl;
	write << "# Particles properties" << endl;
	write << "# Part of particles of type A:            " << part_A << endl;
	write << "# Part of particles of type B:            " << part_B << endl;
	write << "# Part of particles of type O:            " << part_O << endl;
	write << "# Energy A-A bound:                       " << E_AA << endl;
	write << "# Energy B-B bound:                       " << E_BB << endl;
	write << "# Energy O-O bound:                       " << E_OO << endl;
	write << "# Energy A-B or B-A bound:                " << E_AB << endl;
	write << "# Energy A-O or O-A bound:                " << E_AO << endl;
	write << "# Energy B-O or O-B bound:                " << E_BO << endl;
	write << "#" << endl;
	write << "# External field" << endl;
	write << "# External field On:                      " << ext_field_flag << endl;
	write << "# Maximum value of the ext field:         " << ext_field_max_value << endl;
	write << "# Decay of the external filed:            " << ext_field_decay << endl;
	write << "#" << endl;
	write << "# Simulations properties" << endl;
	write << "# Seed:                                   " << seed << endl;
	write << "# Glauber probability used:               " << GlauberProb << endl;
	write << "# Burn-In simulation steps:               " << burn_in_count / (size_i * size_j) << endl;
	write << "# Maximum number of simulation steps:     " << max_step_count / (size_i * size_j) << endl;
	write << "# Steps btw. recordings:                  " << diff_count / (size_i * size_j) << endl;
	write << "# Max. radius (nn_scan):                  " << max_dist << endl;
	write << "# Energy/bounds recorded:                 " << record_total_Energy << endl;
	write << "#" << endl;

}


void Lattice::write_state(std::ofstream &write) const{
	write_parameter(write);
	write << "# State of the lattice after " << step_count / (size_i * size_j) << " time steps:" << endl;
	for(int i = 0 ; i < size_i ; i++){
		for(int j = 0 ; j < size_j ; j++){
			write << state[i][j] << "\t";
		}
		write << endl;
	}
}


bool Lattice::recording_obs_pos()const{
	bool tester;
	if(number_O !=0){tester = true;}
	else{tester = false;}
	return tester;
}


void Lattice::write_obs_pos(std::ofstream &write) const{
	write_parameter(write);
	write << "# Position of the obstacle the system recorded every " << diff_count/(size_i*size_j) << " time steps:" << endl;
	for (unsigned long int i = 0; i < obs_pos_i.size(); i++){
		write << obs_pos_i[i] << "\t" << obs_pos_j[i] << endl ;
	}
}


void Lattice::write_abs_obs_pos(std::ofstream &write) const{
	write_parameter(write);
	write << "# Absolute position of the obstacle the system recorded every " << diff_count/(size_i*size_j) << " time steps:" << endl;
	for (unsigned int i = 0; i < abs_obs_pos_i.size() ; i++){
		write << abs_obs_pos_i[i] << "\t" << abs_obs_pos_j[i] << endl ;
	}
}


void Lattice::write_cluster_size(std::ofstream &write) const{
	write_parameter(write);
	write << "# Average cluster sizes without and with obstacle the system recorded every " << diff_count/(size_i*size_j) << " time steps:" << endl;
	if( number_O == 0 ){
		for(unsigned int i = 0; i < clsize.size(); i++){
			write << setprecision(18) << scientific << clsize[i] << "\t" << 0 << endl;
		}
	}
	else{
		for(unsigned int i = 0; i < clsize.size(); i++){
			write << setprecision(18) << scientific <<clsize[i] << "\t" << clsize_obs[i] << endl ;
		}
	}

}

void Lattice::write_density(std::ofstream &write) const{
	write_parameter(write);
	for(long int i=0; i<size_i*size_j+1; i++){
		write << density[i][0] << "\t" << density[i][1] << endl;
	}
}


bool Lattice::recording_energy()const{
	return record_total_Energy;
}

bool Lattice::recording_density()const{
	bool out;
	if(max_dist>0 && part_O ==1){out = true; }
	else{ out = false; }

	return out;
}

void Lattice::write_energy(std::ofstream &write) const{
	if( record_total_Energy == true ){
		write_parameter(write);
		write << "# Total energy of the system recorded every " << diff_count/(size_i*size_j) << " time steps:" << endl;
		write << setprecision(18) << scientific;
		for(unsigned int i = 0; i < max_step_count/diff_count ; i++){
			write << total_energy[i] << endl;
		}
	}
}


void Lattice::write_bounds(std::ofstream &write) const{
	if( record_total_Energy == true ){
		write_parameter(write);
		write << "# Total number of even and odd bounds in the system recorded every " << diff_count/(size_i*size_j) << " time steps:" << endl;
		for(unsigned int i = 0; i < max_step_count/diff_count ; i++){
			write <<  number_of_even[i] << "\t" << number_of_odd[i] << endl;
		}

		write << endl;
	}
}

void Lattice::write_nn(std::ofstream &write) const{
	write_parameter(write);
	write << "# Relative position of lipid to the obstacle: time count, pos_i, pos_j" << endl;
	for(unsigned long long int i = 0; i < nn_t.size(); i++){
		write << nn_t[i] << "\t" << nn_i[i] << "\t" << nn_j[i] << endl;
	}
}


// used to convert a number (int, double, etc.) into a string
string do_to_str(const double &d){
	stringstream ss;
	ss << d;
	return(ss.str());
}


// used to read in the configuration parameters out of the parameters[i].txt-file
void read_in_parameters(const string &str, double Params[]){
	ifstream in(str.c_str());
	double data;
	int i = 0;


	while(in >> data){
		Params[i]=data;
		i++;
	}
	in.close();
}
