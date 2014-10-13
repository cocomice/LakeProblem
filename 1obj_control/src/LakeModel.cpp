/*
* Lake Problem <another version>
* Yu Li, October 2014
* Cornell University 
* <yl2537@cornell.edu>
*/

#include "LakeModel.h"
using namespace std ;

Lake:: Lake()
{
    b      = 0.62 ;
    q      = 2.0  ;
    alpha  = 0.4  ;
    beta   = 0.08 ;
    delta  = 0.98 ;
    Xo     = 0.0  ;

    npv_util = 0  ;
}

Lake:: Lake(const double param_b, const double param_Xo)
{
	b      = param_b  ;
	Xo     = param_Xo ;
    	
    q      = 2.0  ;
    alpha  = 0.4  ;
    beta   = 0.08 ;
    delta  = 0.98 ;	   

    npv_util = 0  ;
}

Lake:: ~Lake()
{

}

void Lake:: Lake_Setup(const vector<double> & input_var)
{
  int no_ele ;
  no_ele = input_var.size();

  for (int i = 0;  i<no_ele; i++) {
	authro_pol_flow.push_back( input_var.at(i) ) ;
  }
}

void Lake:: Lake_Sim(const int & nYrs)
{

  for (int i=0; i<nYrs; i++){
   // update the state variables 
   if (i==0){
      lake_state.push_back(
    		  Xo*(1-b) + pow(Xo,q)/(1+pow(Xo,q))+ authro_pol_flow.at(i)
    		  );
    }else{
      lake_state.push_back(
    		  lake_state.at(i-1)*(1-b) + ( pow(lake_state.at(i-1), q))/(1 + pow( lake_state.at(i-1), q) )
    		   + authro_pol_flow.at(i)
    		   );
    }	    
  }
}

void Lake:: Util_Cal()
{
	int no_ele ;
	no_ele = lake_state.size();

	for (int i = 0; i < no_ele; ++i){
	    utility_series.push_back(
	    		alpha * authro_pol_flow.at(i) - beta * pow(lake_state.at(i), 2.0 )
	    						);
	    npv_util = npv_util + pow( delta, double(i) ) * utility_series.at(i);
	}
   npv_util = - npv_util ;
}
