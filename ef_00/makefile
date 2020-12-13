all: a.out

a.out: mod_model1.o mod_dynamics_ef.o ehrenfest.o 
	gfortran -o a.out mod_model1.o mod_dynamics_ef.o ehrenfest.o 

%.o: %.f90
	gfortran -c $<

clean:
	rm *.o *.mod a.out
	rm *.100 *.102 *.103 *.200