# Yu Li
# Lake Problem 
# Oct. 21, 2014 at Cornell University


CC2 = mpicc  

main: LakeProblem.exe

LakeProblem.exe: ../../../mt19937ar.c ../../../borg-ms-code/borgms.c LakeModel.cpp LakeProblem_4obj_1_const_control.cpp
	$(CC2) -Wall -o LakeProblem.exe ../../../mt19937ar.c ../../../borg-ms-code/borgms.c \
	LakeModel.cpp LakeProblem_4obj_1_const_control.cpp

clean : 
	rm LakeProblem.exe 

	
	