/*
* Lake Problem <Parallel version>
* Yu Li, October 2014
* Cornell University
* <yl2537@cornell.edu>
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <mpi.h>

// notice here I include the header of master-slave borg
// instead of series borg version (borg.h)
#include "./LakeModel.h"
#include "../../../borg-ms-code/borgms.h" 

using namespace std ;

#define alpha 0.4
#define delta 0.98

#define precis  3
#define samples 100   	// no. of samples for calculating objectives 
#define inertia_thres  (-0.02)  // decision inertia threshold 
#define reliab_thres 	 0.85   // reliability threshold 

int nvars   ; // no. of decision variables 
int nobjs   ; // no. of objectives 
int nconsts ; // no. of constraints 

const int no_years = 100 ; // simulation time horizon
double param_b, param_Xo ; // parameters for defining lake type and initial condition
vector<vector<double> > nat_pol_flow(10000, vector<double>(no_years)) ; // stochastic natural polluted flow 

// Problem definition
void Stoch_Lake_Problem(double * vars, double * objs, double * consts) 
{
	// opt 1 - get randomized index of sample-to-use 
    vector<int> linetouse(samples);
	srand (time(NULL)); //gives a random seed based on run-time, using determined value for random seed analysis
	for (int i=0; i<samples; i++) {
	    //pick a random number based on time
	    //choose 100 of 10,000 available inflow value lines
	    linetouse.at(i) = rand() % 10000;
  	}
	
	// opt 2 - initialize the indicators for calculating final objectives  	
  	double acc_benefit      = 0.0 ; // accumulated benefit over all samples 
  	double acc_reliability  = 0.0 ; // accumulated reliability over all samples 
  	double acc_prob_inertia = 0.0 ; // accumulated probability of maintaining inertia over all samples 
  	vector<double> acc_lake_state(no_years, 0.0) ; // time-series of lake state accumulated over all samples 
	
	// opt 3 - simulation over all samples 
	for (int sample = 0; sample < samples; sample++){
	  int inertia_ctr         = 0   ; // counter for decision inertia ; accumulating over all samples       
      vector<double> nat_flow(no_years, 0.0) ; // initialize natural flow 
      
	  int index = linetouse.at(sample);	  
	  double benefit = 0; // initialize step indicator for benefit 

      vector<double> pol_flow(no_years, 0.0)     ;  // pol_flow = authr_pol_flow + nat_pol_flow ; 	
      vector<double> change_dec(no_years-1, 0.0) ;
      	  
	  for (int i=0; i<no_years; i++){
      	nat_flow.at(i) = nat_pol_flow[index][i]	  ; // draw one sample from natural polluted flow samples 
      	pol_flow.at(i) = vars[i] + nat_flow.at(i) ;
		pol_flow.at(i) = round( pol_flow.at(i)*pow(10,(double)precis) )
						/(pow(10,(double)precis)); //round the value to defined precision

      	if( i>0 ) { // calculate decision change 
      		change_dec.at(i-1) = vars[i] - vars[i-1] ;      
      		change_dec.at(i-1) = round( change_dec.at(i-1)*pow(10,(double)precis) ) 
								/(pow(10,(double)precis)) ;	//round the value to defined precision
      	}
		
		// benefit is independent on the lake state, and thus can be calculated 
		// a priori to the lake simulation 
		double dum_var = vars[i] ;
      	dum_var 	= round( dum_var*pow(10,(double)precis) )/pow(10,(double)precis); 
      	benefit 	= alpha*dum_var ;
      	acc_benefit = acc_benefit + benefit*pow(delta,(i)) ;
      } 
      
	  Lake Yu_Lake(param_b, param_Xo) ;  // create "lake" instance with defined parameters 

	  Yu_Lake.Lake_Setup(pol_flow) ; 
	  Yu_Lake.Lake_Sim(no_years)   ;	  
	  Yu_Lake.Util_Cal() ;

      for (int i=0; i<no_years-1; i++){
  		if (change_dec.at(i) > inertia_thres) inertia_ctr++ ;  		
	  } 
	  acc_prob_inertia = acc_prob_inertia + double(inertia_ctr) / double(no_years-1) ; // cumulate inertia
	  acc_reliability  = acc_reliability + Yu_Lake.reliability ; // cumulate the reliability  	  	  
	  Yu_Lake.Lake_Prober(acc_lake_state) ; // cumulate current lake state vector
	}

    if ( (acc_reliability/samples) > reliab_thres ){
     consts[0]= 0.0;
    }else{
     consts[0] = reliab_thres - (acc_reliability/samples);
    }
  
	double dum_max_ele = *max_element(acc_lake_state.begin(), acc_lake_state.end()) ; // get maximum P value

	// compute the final objectives 
	objs[0] =  dum_max_ele/samples   ;    // minimize the maximum Phosphorous concentration 
	objs[1] = -acc_benefit/samples   ;    // maximize the gain from pollution loads
	objs[2] = -acc_prob_inertia/samples ; // maximize the probability of maintaining inertia
	objs[3] = -acc_reliability/samples  ; // maximize the reliability
}

// Function for reading input data 
void Read_Nat_Flow(string filename, vector<vector<double> > & output )
{
  FILE * myfile ;
  myfile = fopen(filename.c_str(), "r") ;
  
  int linenum = 0 	 ;
  int maxSize = 5000 ; 
  
  if (myfile==NULL){
    perror("Error opening file");
  }else{
      char buffer [maxSize];
      while ( fgets(buffer, maxSize, myfile)!=NULL) 
	     { linenum++;	     	
	         if (buffer[0]!='#')
	       {
	           char *pEnd;
	           char *testbuffer = new char [maxSize];
	             for (int i=0; i <maxSize; i++)
		              testbuffer[i] = buffer[i];
  
	           for (int cols =0;cols<no_years;cols++) // use nDays not nvars, since now they are different
		         {
		            output[linenum-1][cols] = strtod(testbuffer, &pEnd);
		            testbuffer  = pEnd;	
		          }				
	        }
	      }
    }

  fclose(myfile);
}

// Connect to Borg with defined parameter values 
int main(int argc, char* argv[])
{
	param_b  = strtod(argv[1], NULL) ;
	param_Xo = strtod(argv[2], NULL) ;	
	
	// Declare the global variables 
	nobjs 	= 4 ;
	nvars 	= no_years ;
	nconsts = 1 ;

	// Load input data of natual polluted flow
	Read_Nat_Flow("SOWs_Type6.txt", nat_pol_flow);
	
	int rank;
	char runtime[256], outputFilename[256];

	// All master-slave runs need to call startup and set the runtime
	// limits (in hour).
	BORG_Algorithm_ms_startup(&argc, &argv);
	BORG_Algorithm_ms_max_time(0.5);
	BORG_Algorithm_ms_max_evaluations(100000);

	// Create Borg problem
	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconsts, Stoch_Lake_Problem);

	// Set upper and lower bounds 
	for (int i=0; i<nvars; i++) {
		BORG_Problem_set_bounds(problem, i, 0.0, 0.1);
	}

	// Set epsilon value 	
	BORG_Problem_set_epsilon(problem, 0, 0.01);
	BORG_Problem_set_epsilon(problem, 1, 0.01);
	BORG_Problem_set_epsilon(problem, 2, 0.0001);
	BORG_Problem_set_epsilon(problem, 3, 0.0001);	

	// Get the rank of this process.  The rank is used to ensure each
	// parallel process uses a different random seed.
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// When running experiments, we want to run the algorithm multiple
	// times over different seeds and average the results.
	for (int i=0; i<100; i++) { // 100 seeds to be used 
		
		// Save runtime dynamics to a file for diagnostic purpose (optional)
		// Only the master node will write to this file.  
		// Note how we create separate files for each run.
		sprintf(runtime, "runtime_%d.txt", i); 			 // assign the name for runtime output file
		sprintf(outputFilename, "end_of_run_%d.txt", i); // assign the name for end-run output file

		BORG_Algorithm_output_runtime(runtime);
		BORG_Algorithm_output_frequency(500);

		// Seed the random number generator.
		BORG_Random_seed(37*i*(rank+1));

		// Run the master-slave Borg MOEA on the problem.
		BORG_Archive result = BORG_Algorithm_ms_run(problem);

		// Only the master process will return a non-NULL result.
		// Print the Pareto optimal solutions to the screen.
		if (result != NULL) {
			FILE * outputFile = fopen(outputFilename, "w") ;

			if (!outputFile) {
				BORG_Debug("Unable to open final output file\n");
			}

			BORG_Archive_print(result, outputFile);
			BORG_Archive_destroy(result);
			fclose(outputFile) ;
		}
	}

	// Shutdown the parallel processes and exit.
	BORG_Algorithm_ms_shutdown();
	BORG_Problem_destroy(problem);

	return EXIT_SUCCESS;	
}
