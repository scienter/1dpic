#include <stdio.h>
#include <stdlib.h>
#include "laser.h"
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void boostLoadLaser(Domain *D,Laser *L)
{
   float x1,x0,rU,rD,longitudinal;
   int i;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*L->divisionLambda*D->dx;
   rD=L->rD*L->divisionLambda*D->dx;
   x0=-2*rU;

   for(i=1; i<=D->nxSub; i++)
   {
      x1=(i+D->minXSub-1)*D->dx;
      if(x1>=2*x0)
      {
         longitudinal=L->amplitude*sin(2*pi*(x1-x0))*exp(-(x1-x0)*(x1-x0)/rU/rU);
         if(L->polarity==2)
         {
            D->field[i].Pr=longitudinal;            
            D->field[i].Pl=longitudinal;
         }           
      }
   }
}

void loadLaser(Domain *D,Laser *L,External *Ext,double t)
{
   float rU,rD,longitudinal,t0;
   int position,rank;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*L->divisionLambda*D->dt;
   rD=L->rD*L->divisionLambda*D->dt;
   t0=2*rU;
   if(t<=2*rU) {
      t0=2*rU;
      longitudinal=L->amplitude*sin(2.0*pi*(t-t0))*exp(-(t-t0)*(t-t0)/rU/rU);
   }
   else if(t>=2*rU && t<2*rU+2*rD) {
      t0=2*rU;
      longitudinal=L->amplitude*sin(2.0*pi*(t-t0))*exp(-(t-t0)*(t-t0)/rD/rD);
   }
   else if(t>=2*rU+2*rD) {
      longitudinal=0.0;
   }

   if(L->loadPoint>=D->minXSub && L->loadPoint<D->maxXSub)
      rank=myrank;
   else
      rank=nTasks+1;

   position=L->loadPoint+1-D->minXSub;
         

   if(myrank==rank)
   {
      if(L->polarity==1) { 
         D->field[position].E1=longitudinal+Ext->E1;
      }
      else if(L->polarity==2)
      {
         D->field[position].Pr=longitudinal+Ext->Pr;            
         D->field[position].Pl=longitudinal+Ext->Pl;
      }           
      else if(L->polarity==3)
      {
         D->field[position].Sl=longitudinal+Ext->Sl;            
         D->field[position].Sr=longitudinal+Ext->Sr; 
      }
   }
}
