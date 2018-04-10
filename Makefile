#Objet
OBJ= fonction.o Sequentiel.o
OBJ_parra= parallele.f90
optn= -g -fbounds-check -o0 -Wall

#Compile et crée l'executable
sequentiel: $(OBJ)
	gfortran $(optn) -o exe $(OBJ)

parrallele:
	mpif90 -o para $(OBJ_parra)

%.o: %.f90
	gfortran -c $<

#Cleaner
clean:
	rm -r *~ *.o *.mod Data/
