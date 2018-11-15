/* $Id: phg_rungekutta3.C,v 1.1.1.1 2015/12/24 07:00:56 zjcao Exp $ */
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <map>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <map.h>
#endif


#include "phg.h"

void phg_rungekutta3(const FLOAT dt,
		     DOF **dofs0,DOF **dofs1,DOF **dofsrhs,
		     const int rk3)
{
//   cout<<dofs0[0]->name<<", "<<dofs1[0]->name<<", "<<dofsrhs[0]->name<<endl;
   const FLOAT F1o6=1.0/6,HLF=0.5,TWO=2.0;
   if(rk3==0)
   {
      for(int i=0;i<VN;i++)
      {
        phgDofAXPBY(1.0,dofsrhs[i],0     ,&dofs1[i]);
        phgDofAXPBY(1.0,dofs0[i]  ,HLF*dt,&dofs1[i]);
      }
   }
   else if(rk3==1)
   {
      for(int i=0;i<VN;i++)
      {
        phgDofAXPBY( 4.0,  dofs1[i],1.0   ,&dofsrhs[i]);
        phgDofAXPBY(-1.0,dofsrhs[i],6.0   ,&dofs1[i]);
        phgDofAXPBY( 1.0,  dofs0[i],dt    ,&dofs1[i]);
      }
   }
   else if(rk3==2)
   {
      for(int i=0;i<VN;i++)
      {
        phgDofAXPBY(1.0,dofsrhs[i],1.0    ,&dofs1[i]);
        phgDofAXPBY(1.0,dofs0[i]  ,F1o6*dt,&dofs1[i]);
      }
   }
   else
   {
      cout<<__FILE__<<"("<<__LINE__<<"), "<<__func__
	  <<": something is wrong in RK3 counting!!"<<endl;
      assert(false);
   }
//   cout<<dofs0[0]->name<<", "<<dofs1[0]->name<<", "<<dofsrhs[0]->name<<endl;
}
