# Yu Li
# Lake Problem 
# Oct. 21, 2014 at Cornell University


CC 	 = gcc
CC2  = g++

main: LakeProblem4obj_control BorgExec

LakeProblem4obj_control : LakeModel.o moeaframework.o LakeProblem_4obj_1_const_control.cpp
	$(CC2) -Wall -o LakeProblem4obj_control.exe LakeModel.o moeaframework.o LakeProblem_4obj_1_const_control.cpp -lm
	rm *.o 

BorgExec : ../../../frontend.c	mt19937ar.o borg.o
	$(CC) -Wall -o BorgExec.exe ../../../frontend.c borg.o mt19937ar.o -lm
	rm borg.o mt19937ar.o

LakeModel.o : LakeModel.cpp LakeModel.h
	$(CC2) -Wall -o LakeModel.o -c LakeModel.cpp -lm

moeaframework.o: ../../../moeaframework.c ../../../moeaframework.h
	$(CC) -Wall -o moeaframework.o -c ../../../moeaframework.c -lm
	
borg.o : ../../../borg.c ../../../borg.h
	$(CC) -Wall -o borg.o -c ../../../borg.c -lm

mt19937ar.o : ../../../mt19937ar.c ../../../mt19937ar.h 
	$(CC) -Wall -o mt19937ar.o -c ../../../mt19937ar.c 	-lm
	
clean:
	rm BorgExec LakeProblem4obj_control 	