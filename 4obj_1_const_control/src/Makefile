# Yu Li
# Lake Problem 
# Oct. 21, 2014 at Cornell University


CC  = gcc  
CC2 = g++  

main: LakeProblem.exe

LakeProblem.exe: mt19937ar.o borg.o LakeModel.o LakeProblem_4obj_1_const_control.cpp
	$(CC2) -Wall -o LakeProblem.exe mt19937ar.o borg.o LakeModel.o LakeProblem_4obj_1_const_control.cpp
	rm *.o
	
mt19937ar.o: ../../../mt19937ar.c ../../../mt19937ar.h
	$(CC) -Wall -c ../../../mt19937ar.c -o mt19937ar.o
	
borg.o : ../../../borg.c  ../../../borg.h
	$(CC) -Wall -c ../../../borg.c -o borg.o 
	
LakeModel.o: LakeModel.cpp LakeModel.h
	$(CC2) -Wall -o LakeModel.o -c LakeModel.cpp	

clean :
	rm LakeProblem.exe
	
	