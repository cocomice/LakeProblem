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

using namespace std; 

class Lake {

  private: 
    double b      ;
    double q      ;
    double alpha  ;
    double beta   ;
    double delta  ;
    double Xo     ;

	vector<double> lake_state ;
	vector<double> authro_pol_flow ;

  public:
	Lake()  ;
	Lake(const double param_b, const double param_Xo) ;
	virtual ~Lake() ;


	vector<double> utility_series  ;
	double npv_util ;

	void Lake_Setup(const vector<double> & input_var) ;

	void Lake_Sim(const int & nYrs) ;

    void Util_Cal() ;
} ;
