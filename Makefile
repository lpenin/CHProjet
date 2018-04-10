#Objet
OBJ= fonction.o Sequentiel.o
OBJ_para= fonction_parallele.o parallele.o
optn= -g -fbounds-check -o0 -Wall

#Compile et crée l'executable
sequentiel: $(OBJ)
	gfortran $(optn) -o exe $(OBJ)

parallele:
	mpif90 -c fonction_parallele.f90 parallele.f90
	mpif90 -o para $(OBJ_para)

%.o: %.f90
	gfortran -c $<

#Cleaner
clean:
	rm -r *~ *.o *.mod Data/
