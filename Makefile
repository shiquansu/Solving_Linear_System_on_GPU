FC=ifx #fortran compiler
CC=icpx #C++/SYCL compiler

COREFLAGS=-qopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device pvc" -fsycl #AOT
#COREFLAGS=-qopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend "-device 0x0bd6" -fsycl #AOT
#COREFLAGS=-qopenmp -fopenmp-targets=spir64 -fsycl #default, JIT
CFFLAGS= -O3 -i8 -DMKL_ILP64 -free -fpp $(COREFLAGS) # OpenMP must have -i8 to allocate host memory correctly
LFFLAGS= -O3 -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $(COREFLAGS)
#LFFLAGS=     -L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl $(CFFLAGS)
CCFLAGS= -O3 -DSYCL_DEVICES_gpu -fsycl -DMKL_ILP64  -qmkl=sequential

SRC_F=solve_complex_matrix.f90 matrix_inverse_omp_onemkl.f90 matrix_inverse_omp_onemkl_batch.f90 matrix_inverse_wrapper_c.f90 init_from_file.f90 init_random.f90 init_trivial.f90 make_batch_identical.f90 check_inverse.f90 check_inverse_batch.f90 

OBJ_C=getri.o
OBJ_F=${SRC_F:.f90=.o} #substitute .f90 with .o


OBJ_solver_c_wrapper=solve_complex_matrix.o init_trivial.o matrix_inverse_wrapper_c.o check_inverse.o
OBJ_solver_omp_dispatch=solve_complex_matrix.o init_trivial.o matrix_inverse_omp_onemkl.o check_inverse.o
OBJ_solver_omp_dispatch_batch=solve_complex_matrix.o init_trivial.o make_batch_identical.o matrix_inverse_omp_onemkl_batch.o check_inverse_batch.o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(CFFLAGS) -o $@ -c $<
%.o: %.cpp #wildcard rule, creation of *.o depends on *.f90
	$(CC) $(CCFLAGS) -o $@ -c $<
#getri.o: getri.cpp #wildcard rule, creation of *.o depends on *.f90
#	$(CC) $(CCFLAGS) -o getri.o -c getri.cpp

all_objs: $(OBJ_F) $(OBJ_C)	
	echo "compile all objs: F obj="$(OBJ_F)"C obj="$(OBJ_C) 

#scm: $(OBJ_solver_c_wrapper) $(OBJ_C)
#scm: $(OBJ_solver_omp_dispatch) $(OBJ_C)
scm: $(OBJ_solver_omp_dispatch_batch) $(OBJ_C)
	@echo "linking: F and C obj="$^ 
	$(FC) $(LFFLAGS) -o $@ $^

## example of use explicit rule for individual source file, by default, use the same wildcard rule for all source files
#matrix_inverse_omp_onemkl_batch.o: matrix_inverse_omp_onemkl_batch.f90 
#	$(FC) $(CFFLAGS) -o $@ -c $<
#getri.o: getri.cpp
#	$(CC) $(CCFLAGS) -o $@ -c $< 

clean: #cleans all the old compilation files
	rm -f *.modmic *.mod *.o core.* scm 
