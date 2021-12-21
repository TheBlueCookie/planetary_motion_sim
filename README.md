# planetary_motion_sim
Fortran 90 based simulation of planetary motion using Newtonian mechanics and Verlet algorithm. Visualization implemented in Python 3.9

By Elias Ankerhold

COMPILATION:
Compile the main program using one of three options:

(1) run "make" in the main directory where Makefile is located -->  /src/fortran_sim/sim.out will be created

(2) manually compile to modules, then combine in final compilation step. run the following commands in /src/fortran_sim:

gfortran -o vector_framework.o -c vector_framework.f90
gfortran -o i_o_framework.o -c i_o_framework.f90
gfortran -o plan_sim_framework.o -c plan_sim_framework.f90
gfortran -o simulation_starting.o -c simulation_starting.f90
gfortran -o start_sim.out vector_framework.o i_o_framework.o plan_sim_framework.o simulation_starting.o

(3) compile all at once, run the following command in /src/fortran_sim:

gfortran -o sim.out vector_framework.f90 i_o_framework.f90 plan_sim_framework.f90 simulation_starting.f90

To compile the initial value helper coordinate_conversion.f90, run the following command in /src/fortran_sim:

gfortran -o init_params.out coordinate_conversion.f90

EXECUTION:

Make sure the input file [my_file.txt] complies with the format specifications described in report.pdf and is located in /run/ 
Navigate to /src/fortran_sim and xecute the program using:

./start_sim.out my_file.txt
