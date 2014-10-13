/*
* Lake Problem <another version>
* Yu Li, October 2014
* Cornell University
* <yl2537@cornell.edu>
*/

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "./LakeModel.h"
#include "../borg.h"

using namespace std;


int nvars = 100 ;
int nobjs = 1   ;

const int no_years = 100 ;

double param_b, param_Xo ;

// Problem definition
void Lake_Problem(double * vars, double * objs, double * consts) 
{	
	vector<double> pol_flow  ;

	for(int i=0; i<nvars; i++) pol_flow.push_back(vars[i]) ;

	// define the type of the lake
	Lake Yu_Lake(param_b, param_Xo) ; 

	Yu_Lake.Lake_Setup(pol_flow) ;
	Yu_Lake.Lake_Sim(no_years) ;
	Yu_Lake.Util_Cal() ;

	objs[0] = Yu_Lake.npv_util ;
}

// Connect to Borg with defined parameter values 
int main(int argc, char* argv[])
{
	// add error check message in future version
	
	// char* pEnd ;
	param_b  = strtod(argv[1], NULL) ;
	param_Xo = strtod(argv[2], NULL) ;	
			
	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, 0, Lake_Problem);

	for (int i=0; i<nvars; i++) {
		BORG_Problem_set_bounds(problem, i, 0.0, 0.1);
	}

	for (int i=0; i<nobjs; i++) {
		BORG_Problem_set_epsilon(problem, i, 0.01);
	}

	BORG_Archive result = BORG_Algorithm_run(problem, 10000);

	BORG_Archive_print(result, stdout);

	BORG_Archive_destroy(result) ;
	BORG_Problem_destroy(problem);

	return EXIT_SUCCESS;	
}
