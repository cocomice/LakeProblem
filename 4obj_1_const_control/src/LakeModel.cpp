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
    beta   = 0.08 ;
    Xo     = 0.0  ;
	//delta  = 0.98 ;
	//alpha  = 0.4  ;

	reliability = 0 ;
    //benefit     = 0 ;    
    //npv_util    = 0 ;
}

Lake:: Lake(const double param_b, const double param_Xo)
// can be improved to modify all parameters 
{
    b      = param_b  ;
	Xo     = param_Xo ;
	
    q      = 2.0  ;
	beta   = 0.08 ;
    //alpha  = 0.4  ;
    //delta  = 0.98 ;	  

	reliability = 0 ;
    //benefit     = 0 ;
    //npv_util    = 0 ;
}

Lake:: ~Lake()
{

}

void Lake:: Lake_Setup(const vector<double> & input_var)
{
  int no_ele ;
  no_ele = input_var.size();

  for (int i = 0;  i<no_ele; i++) {
	pol_flow.push_back( input_var.at(i) ) ;
  }
}

void Lake:: Lake_Sim(const int & nYrs)
{

  for (int i=0; i<nYrs; i++){
   // update the state variables 
   if (i==0){
      lake_state.push_back(
    		  Xo*(1-b) + pow(Xo,q)/(1+pow(Xo,q))+ pol_flow.at(i)
    		  );
    }else{
      lake_state.push_back(
    		  lake_state.at(i-1)*(1-b) + ( pow(lake_state.at(i-1), q))/(1 + pow(lake_state.at(i-1),q) )
    		   + pol_flow.at(i)
    		   );
    }	
  }
}

void Lake:: Util_Cal()
{
  int no_ele = lake_state.size();
  int ctr    = 0  ; 	

	for (int i = 0; i < no_ele; ++i){
	    //utility_series.push_back( alpha * pol_flow.at(i) ) ; 	    						
	    //benefit = benefit + utility_series.at(i) * pow(delta, (double)i);

      if( lake_state.at(i) < pcrit ) ctr++ ;      
	}   
  reliability = double(ctr) / double(no_ele) ;
  // update constraints 

}

void Lake:: Lake_Prober(vector<double> & output)
{
  const int no_ele = lake_state.size() ;
  
  if (output.size()!=(unsigned)no_ele){
    cout << "argument of 'Lake_Prober' must be " << no_ele << endl ;
  }

  for (int i = 0; i < no_ele; ++i)
  {
    output.at(i) = output.at(i) + lake_state.at(i) ;
  }
}
