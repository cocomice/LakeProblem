/*
* Lake Problem <another version>
* Yu Li, October 2014
* Cornell University 
* <yl2537@cornell.edu>
*/


#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#define pcrit 		 	0.5
#define reliab_thres 	0.85   // reliability threshold 

using namespace std; 

class Lake {

  private: 
    double b, q, beta ;    //alpha,delta
    double Xo     ;

	vector<double> lake_state ;
	vector<double> pol_flow ;

  public:
	Lake()  ;
	Lake(const double param_b, const double param_Xo) ;
	virtual ~Lake() ;


	//vector<double> utility_series  ;
	double  reliability; //benefit  npv_util ; 

	void Lake_Setup(const vector<double> & input_var) ;
	void Lake_Sim(const int & nYrs) ;
    void Util_Cal() ;
    void Lake_Prober(vector<double> & output) ;
} ;
