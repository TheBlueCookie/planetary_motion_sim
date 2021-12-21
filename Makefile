# compiler
FC=gfortran
# compiler flags
FFLAGS=-fcheck=all
PNAME1=simulation
OUTNAME1=start_sim.out

all: $(PNAME1) clean

SRC1=src/fortran_sim/vector_framework.f90 src/fortran_sim/i_o_framework.f90 src/fortran_sim/plan_sim_framework.f90 src/fortran_sim/simulation_starting.f90

OBJ1=${SRC1:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(PNAME1): $(OBJ1)
	$(FC) $(FFLAGS) -o src/fortran_sim/$(OUTNAME1) $(OBJ1)

clean:
	@rm -f *.mod src/fortran_sim/*.o $(PNAME1)
