MKLROOT=/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl

all: a.out


a.out: mod_model2.o mod_dynamics_sqc.o sqc.o 
	ifort -o a.out mod_model2.o mod_dynamics_sqc.o sqc.o -qopt-matmul -ipo -O3  -no-prec-div -static-intel -fp-model fast=2 -xHost -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include   ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

%.o: %.f90
	ifort  -c $< -qopt-matmul -ipo -O3  -no-prec-div -static-intel -fp-model fast=2 -xHost

quick:
	gfortran -o a.out mod_model2.f90 mod_dynamics_sqc.f90 sqc.f90

quicker: 
	gfortran -o a.out mod_model2.f90 mod_dynamics_sqc.f90 sqc.f90 ~/lapack-3.8.0/liblapack.a ~/lapack-3.8.0/librefblas.a

qlocal:
	gfortran -o a.out mod_model2.f90 mod_dynamics_sqc.f90 sqc.f90 -L/usr/local/lib -llapack 

clean:
	rm *.o *.mod a.out

