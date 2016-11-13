#include "b0j3.h"
#include "b0j4.h"
#include "b1j3.h"
#include "b1j4.h"
#include "b2j3.h"
#include "b2j4.h"
#include "b3j3.h"
#include "b3j4.h"
#include "b4j4.h"

int main(int argc, char *argv[])
{
   std::string imva = std::string(argv[1]);

   if( imva == "HutSTallb0j3" || imva == "HutSTeb0j3" || imva == "HutSTmb0j3" ||
       imva == "HctSTallb0j3" || imva == "HctSTeb0j3" || imva == "HctSTmb0j3" ||
       imva == "HutTTallb0j3" || imva == "HutTTeb0j3" || imva == "HutTTmb0j3" ||
       imva == "HctTTallb0j3" || imva == "HctTTeb0j3" || imva == "HctTTmb0j3")
     b0j3(imva);

   if( imva == "HutSTallb0j4" || imva == "HutSTeb0j4" || imva == "HutSTmb0j4" ||
       imva == "HctSTallb0j4" || imva == "HctSTeb0j4" || imva == "HctSTmb0j4" ||
       imva == "HutTTallb0j4" || imva == "HutTTeb0j4" || imva == "HutTTmb0j4" ||
       imva == "HctTTallb0j4" || imva == "HctTTeb0j4" || imva == "HctTTmb0j4")
     b0j4(imva);
   
   if( imva == "HutSTallb1j3" || imva == "HutSTeb1j3" || imva == "HutSTmb1j3" ||
       imva == "HctSTallb1j3" || imva == "HctSTeb1j3" || imva == "HctSTmb1j3" ||
       imva == "HutTTallb1j3" || imva == "HutTTeb1j3" || imva == "HutTTmb1j3" ||
       imva == "HctTTallb1j3" || imva == "HctTTeb1j3" || imva == "HctTTmb1j3")
     b1j3(imva);

   if( imva == "HutSTallb1j4" || imva == "HutSTeb1j4" || imva == "HutSTmb1j4" ||
       imva == "HctSTallb1j4" || imva == "HctSTeb1j4" || imva == "HctSTmb1j4" ||
       imva == "HutTTallb1j4" || imva == "HutTTeb1j4" || imva == "HutTTmb1j4" ||
       imva == "HctTTallb1j4" || imva == "HctTTeb1j4" || imva == "HctTTmb1j4")
     b1j4(imva);

   if( imva == "HutSTallb2j3" || imva == "HutSTeb2j3" || imva == "HutSTmb2j3" ||
       imva == "HctSTallb2j3" || imva == "HctSTeb2j3" || imva == "HctSTmb2j3" ||
       imva == "HutTTallb2j3" || imva == "HutTTeb2j3" || imva == "HutTTmb2j3" ||
       imva == "HctTTallb2j3" || imva == "HctTTeb2j3" || imva == "HctTTmb2j3")
     b2j3(imva);

   if( imva == "HutSTallb2j4" || imva == "HutSTeb2j4" || imva == "HutSTmb2j4" ||
       imva == "HctSTallb2j4" || imva == "HctSTeb2j4" || imva == "HctSTmb2j4" ||
       imva == "HutTTallb2j4" || imva == "HutTTeb2j4" || imva == "HutTTmb2j4" ||
       imva == "HctTTallb2j4" || imva == "HctTTeb2j4" || imva == "HctTTmb2j4")
     b2j4(imva);

   if( imva == "HutSTallb3j3" || imva == "HutSTeb3j3" || imva == "HutSTmb3j3" ||
       imva == "HctSTallb3j3" || imva == "HctSTeb3j3" || imva == "HctSTmb3j3" ||
       imva == "HutTTallb3j3" || imva == "HutTTeb3j3" || imva == "HutTTmb3j3" ||
       imva == "HctTTallb3j3" || imva == "HctTTeb3j3" || imva == "HctTTmb3j3")
     b3j3(imva);

   if( imva == "HutSTallb3j4" || imva == "HutSTeb3j4" || imva == "HutSTmb3j4" ||
       imva == "HctSTallb3j4" || imva == "HctSTeb3j4" || imva == "HctSTmb3j4" ||
       imva == "HutTTallb3j4" || imva == "HutTTeb3j4" || imva == "HutTTmb3j4" ||
       imva == "HctTTallb3j4" || imva == "HctTTeb3j4" || imva == "HctTTmb3j4")
     b3j4(imva);

   if( imva == "HutSTallb4j4" || imva == "HutSTeb4j4" || imva == "HutSTmb4j4" ||
       imva == "HctSTallb4j4" || imva == "HctSTeb4j4" || imva == "HctSTmb4j4" ||
       imva == "HutTTallb4j4" || imva == "HutTTeb4j4" || imva == "HutTTmb4j4" ||
       imva == "HctTTallb4j4" || imva == "HctTTeb4j4" || imva == "HctTTmb4j4")
     b4j4(imva);
}
