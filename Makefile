FC=ifx #fortran compiler
CC=icpx #C++/SYCL compiler

COREFLAGS=-qopenmp -fopenmp-targets=spir64 -fsycl
CFFLAGS=-i8 -DMKL_ILP64 -i_tapi -free -fpp $(COREFLAGS) 
LFFLAGS=-i_tapi -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $(CFFLAGS)
#LFFLAGS=-i_tapi -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $(COREFLAGS)

CCFLAGS=-g -O2 -DSYCL_DEVICES_gpu -fsycl -DMKL_ILP64  -qmkl=sequential

SRC_F=solve_complex_matrix.f90 matrix_inverse_omp_onemkl.f90 matrix_inverse_omp_onemkl_batch.f90 matrix_inverse_wrapper_c.f90 init_from_file.f90 init_from_file_batch.f90 init_random.f90 check_inverse.f90 check_inverse_batch.f90 

OBJ_C=getri.o
OBJ_F=${SRC_F:.f90=.o} #substitute .f90 with .o

OBJ_solver=solve_complex_matrix.o init_from_file.o matrix_inverse_omp_onemkl.o check_inverse.o
#OBJ_solver=solve_complex_matrix.o init_from_file.o matrix_inverse_wrapper_c.o check_inverse.o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(CFFLAGS) -o $@ -c $<
%.o: %.cpp #wildcard rule, creation of *.o depends on *.f90
	$(CC) $(CCFLAGS) -o $@ -c $<

all_objs: $(OBJ_F) $(OBJ_C)	
	echo "make all objs: F obj="$(OBJ_F)"C obj="$(OBJ_C) 

solve_complex_matrix: $(OBJ_solver) $(OBJ_C)
	@echo "F obj="$(OBJ_solver)", C obj="$(OBJ_C) 
	$(FC) $(LFFLAGS) -o $@ $(OBJ_solver) $(OBJ_C)

## example of use explicit rule for individual source file, by default, use the same wildcard rule for all source files
#matrix_inverse_omp_onemkl_batch.o: matrix_inverse_omp_onemkl_batch.f90 
#	$(FC) $(CFFLAGS) -o $@ -c $<
#getri.o: getri.cpp
#	$(CC) $(CCFLAGS) -o $@ -c $< 

clean: #cleans all the old compilation files
	rm -f *.modmic *.mod *.o solve_complex_matrix 
