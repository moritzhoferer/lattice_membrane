#!/usr/bin/ipython

# use this program to construct the inputfiles for basic_ising (v1.0)
# creates also a file "option.txt". this is neccessary for running jobs
# with several task on the small cluster @ UniMI

import os
import numpy as np
import sys

os.makedirs('parameter_files/', exist_ok=True)

# system configuartion
# verticle size of the system (# of rows)
size_i = int(sys.argv[1])

# horizontal size of the system (# of columns)
size_j = size_i

# periodic boundary condtions in vertical direction
periodic_bc_i = 1

# periodic boundary condtions in horizontal direction
periodic_bc_j = 1

# Boltzmann constan
kB = 1

# temperatur
T 						=  1


# system properties

# part of particles of type A / +1
part_A = .125 * (size_i*size_j) - int(sys.argv[2])

# part of particles of type B / -1
part_B = .875 * (size_i*size_j)

# part of particles of type O / 0
part_O = int(sys.argv[2])

# energy of to neighbouring particles of type A
E_AA = -1

# energy of to neighbouring particles of type B
E_BB = -1 # E_AA

# energy of to neighbouring particles of type O
E_OO = 0 

# energy of to neighbour particles one with type A and one with type B
E_AB = 1 # -E_AA

# energy of to neighbour particles one with type A and one with type O
E_AO = 0 # -E_OO

# energy of to neighbour particles one with type B and one with type O
E_BO = 0 #  E_OO



# external may automatically positioned in the middle of the lattice
# external field on or off
ext_field_flag = 0

# max value
ext_field_max_value	= 0 

# homogenous ( = 0 ) or with decay ( > 0 )
ext_field_decay	= 0 


# simulation properties
# seed used to make simulations reproducible
seed 					= 13
np.random.seed(seed)

# selection methode: if =1 Glauber and if = 0 Metropolis is used
GlauberProb = 1

# maximum number of executed simulations steps
# pow(10, max_step_count) * size_i * size_j || JUST ENTER THE EXPONENT
max_step_count = int(4)

# number of simulation steps to executed before the average occupation recording is activated
# pow(10, burn_in_count) * size_i * size_j || JUST ENTER THE EXPONENT
burn_in_count = int(6)

# number of simulation steps to pass by untill the next snap shot (used for total energy and the number of different bounds)
# pow(10, diff_count) * size_i * size_j || JUST ENTER THE EXPONENT
diff_count				= int(0)

# used for angle and distance dependent distribution function ( 0 off, if > 0 cut off radius )
record_Density = 10

# the total energy and the number of different bound is recorded
record_Energy = 0


# choosen : "E_BO"
if(part_O == 0):
	var = np.array([0])
else:
	var_min = .4 
	var_max = 1.4
	var_delta = .1
	var = np.arange( var_min , var_max + var_delta/2 , var_delta )
	#var = 2.**np.arange(-1,9)
	#var = np.array([ .5, 1.])
	#var = np.linspace( .05,.4 , 8)

# choosen : "E_AB"
if(part_O == 0):
	var2_min = .35
	var2_max = .55
	var2_delta = 1/200.
	var2 = np.arange( var2_min , var2_max + var2_delta/2 , var2_delta )
	#var2 = np.arange( 100, 151, 1)/50.
else:
	var2_min = .05
	var2_max = .4
	var2_delta = .05
	var2 = np.arange( var2_min , var2_max + var2_delta/2 , var2_delta )
	#var2 = np.array([ .05, .2, .35])
	#var2 = np.linspace( .05,.4 , 8)
	

file= open('./parameter_files/option.txt','w+')
for i in range(len(var)*len(var2)):
	file.write('{} \n'.format(i+1))

file.close()

c = 1
for j in range(len(var2)):
	# !!! var2 to be changed !!!
	#T = var2[j]
	E_AA = -var2[j]
	E_BB = -var2[j]
	E_AB =  var2[j]
	#E_BO =  var2[j]

	for i in range(len(var)):
		# !!! var to be changed !!!
		#E_AA = -var2[i]
		#E_BB = -var2[i]
		#E_AB =  var2[i]
		E_OO = -var[i]
		E_AO = -var[i]
		E_BO =  var[i]
		
		sim_seed = np.random.randint(1,1000)

		# open new file
		file = open('./parameter_files/Parameters{}.txt'.format(c),'w+')

		# lattice properties
		file.write('{} \n'.format(size_i))
		file.write('{} \n'.format(size_j))
		file.write('{} \n'.format(periodic_bc_i))
		file.write('{} \n'.format(periodic_bc_j))
		file.write('{} \n'.format(kB))
		file.write('{} \n'.format(T))

		# particles properties
		file.write('{} \n'.format(part_A))
		file.write('{} \n'.format(part_B))
		file.write('{} \n'.format(part_O))
		file.write('{} \n'.format(E_AA))
		file.write('{} \n'.format(E_BB))
		file.write('{} \n'.format(E_OO))
		file.write('{} \n'.format(E_AB))
		file.write('{} \n'.format(E_AO))
		file.write('{} \n'.format(E_BO))

		# external field
		file.write('{} \n'.format(ext_field_flag))
		file.write('{} \n'.format(ext_field_max_value))
		file.write('{} \n'.format(ext_field_decay))

		# simulation properties
		file.write('{} \n'.format(sim_seed))
		file.write('{} \n'.format(GlauberProb))
		file.write('{} \n'.format(max_step_count))
		file.write('{} \n'.format(burn_in_count))
		file.write('{} \n'.format(diff_count))
		file.write('{} \n'.format(record_Density))
		file.write('{}'.format(record_Energy))

		file.close()
		c += 1
