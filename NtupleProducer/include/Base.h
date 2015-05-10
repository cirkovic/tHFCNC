#ifndef BASE_H
#define BASE_H

#include "TObject.h"
#include "TLorentzVector.h"

extern unsigned int idx;

class Base : public TLorentzVector
{
 public:
   
   Base();
   virtual ~Base();
   
   template< class BR >
     bool CHECK(BR br)
       {
	  bool res = false;   
	  if( br )
	    {		    
	       if( idx < br->size() )
		 res = true;
	    }	       
	  return res;
       }
   
//   std::pair<int,bool> truthGetParent(int p);
   float GetDPhi(float phi1,float phi2);
   float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
   
 private:

   ClassDef(Base,1)
};   

#endif
