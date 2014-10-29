Lake Problem Training（C++)
===========
Lake Problem based on Carpenter et al 1999.

Intended for use with MOEAFramework and Borg MOEA. For downloading those libraries please visit: `http://borgmoea.org`;

##Branches: 
* **master**: the main branch of the lake problem model, using the embedded coding approach combining the declaration of the 4-objective stochstic lake problems and the Borg algorithm within one file. 
* **diagnostic**: this version is intended for comparative study and hence coded independently from the Borg algorithm, which means one will get a standalone executable that runs given specfied input variables. To run the optimization requires the interconnection with series Borg in the command line. See 'To compile and run' session for more details.
* **lakeproblem_par**: the parallel version using embedded form of coding. To compile this version requires the Borg Master-Slave version. 

**Note**: 
Three branches take the same form of folder structure (see 'Contents'), and the differences between three branches are only within the 4-objective control problem. The single objective problem is same across all three branches. 

##Contents:
###Single objective lake problem
* `LakeProblem_1obj_control.cpp`: C++ source code for the 1 objective formulation.
* `LakeProblem.h`: header file for declaration of `Lake` class.
* `LakeProblem.cpp`: source code for `Lake` class.

###Four objective one control problem formulation
* `LakeProblem_4obj_1_const_control.cpp`: C++ source code for the 4 objectives formulation.
* `LakeProblem.h`: header file for declaration of `Lake` class.
* `LakeProblem.cpp`: source code for `Lake` class. 
* `Makefile`: makefile for compilation. 

##To compile and run:
To compile and run the program requires the Borg source codes, which can be obtained from here: `http://borgmoea.org`;

* Type the following command `make -f Makefile` to compile the model; for single objective formulation, type `g++ -o lake_executable LakeProblem_1obj_control.cpp LakeProblem.cpp` ;
* Before you can use the lake model, you will need SOW (state of world) file named `SOWs_Type6.txt`.  This file is essentially a matrix containing [10000 x 100] samples generated from a specified lognormal distribution function. One can use any method or software, e.g. Matlb, to generate this state of world and simply save it as 'SOWs_Type6.txt';
* There are two way to run the program to optimize the lake problem with Borg. The first version refers to the separate coding style, where one have to compile the Borg execetuable and the problem seperately and interconnect them via command line (see below). The last two arguements (.42 and 0.0) specify the decay rate and initial concentration of the lake.

`./Borg_executable.exe -v 100 -o 4 -c 1 -l 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 -u 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 -e 0.01,0.01,0.0001,0.0001 -n 10000 -f outputfile.set -- ./lake_executable.exe .42 0.0`

* for embedded form, simply use `./name_of_executable.exe .42 0.0`. You can append `> output.txt` command to pipe out your output into a file named `output.txt` ;
* 
**Note**: please check carefully the path when including the `borg.h` or `moeaframework.h` library; For any question, please send me the email: yli@elet.polimi.it ;
