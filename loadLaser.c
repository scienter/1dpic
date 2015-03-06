#include <stdio.h>
#include <stdlib.h>
#include "laser.h"
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void boostLoadLaser(Domain *D,LaserList *L)
{
   float x1,x0,rU,rD,longitudinal,kfactor;
   int i;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*D->divisionLambda*D->dx;
   rD=L->rD*D->divisionLambda*D->dx;
   x0=-2*rU;
   kfactor=L->lambda/D->lambda;

   for(i=1; i<=D->nxSub; i++)
   {
      x1=(i+D->minXSub-1)*D->dx;
      if(x1>=2*x0)
      {
         longitudinal=L->amplitude*sin(2*pi*kfactor*(x1-x0))*exp(-(x1-x0)*(x1-x0)/rU/rU);
         if(L->polarity==2)
         {
            D->field[i].Pr=longitudinal;            
            D->field[i].Pl=longitudinal;
         }           
      }
   }
}

void loadLaser(Domain *D,LaserList *L,External *Ext,double t)
{
   float rU,rD,longitudinal,t0,fratio,flat;
   int position,rank,i;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

/*
   rU=L->rU*L->lambda/D->lambda;
   rD=L->rD*L->lambda/D->lambda;
   flat=L->flat*L->lambda/D->lambda;
*/
  // 20150220 mshur
   rU=L->rU;
   rD=L->rD;
   flat=L->flat;
  //--------------
  //
   fratio=L->omega/D->omega;

   if(t<=3*rU) {
      t0=3*rU; 
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0))*exp(-(t-t0)*(t-t0)/rU/rU);
   }
   else if(t>=3*rU && t<3*rU+flat)  {
      t0=3*rU;
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0));
   }
   else if(t>=3*rU+flat && t<3*rU+3*rD+flat) {
      t0=3*rU+flat;
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0))*exp(-(t-t0)*(t-t0)/rD/rD);
   }
   else if(t>=3*rU+3*rD+flat) 
      longitudinal=0.0;


   if(L->loadPoint>=D->minXSub && L->loadPoint<=D->maxXSub)
      rank=myrank;
   else
      rank=nTasks+1;
   position=L->loadPoint-D->minXSub+1;

   if(myrank==rank)
   {
      if(L->polarity==1) { 
         D->field[position].E1=longitudinal;
      }
      else if(L->polarity==2)
      {
         D->field[position].Pr=longitudinal;            
         D->field[position].Pl=longitudinal;
      }           
      else if(L->polarity==3)
      {
         D->field[position].Sl=longitudinal;            
         D->field[position].Sr=longitudinal; 
      }
      else if(L->polarity==4)	//right circular
      {
         D->field[position].Pr=longitudinal;            
         D->field[position].Pl=longitudinal; 
         D->field[position].Sl=longitudinal/tan(2.0*pi*fratio*(t-t0));            
         D->field[position].Sr=longitudinal/tan(2.0*pi*fratio*(t-t0)); 
      }
      else if(L->polarity==5)	//left circular
      {
         D->field[position].Pr=longitudinal;            
         D->field[position].Pl=longitudinal; 
         D->field[position].Sl=-longitudinal/tan(2.0*pi*fratio*(t-t0));            
         D->field[position].Sr=-longitudinal/tan(2.0*pi*fratio*(t-t0)); 
      }
   }
}

void addLoadLaser(Domain *D,LaserList *L,External *Ext,double t)
{
   float rU,rD,longitudinal,t0,fratio,flat,phs;
   int position,rank,i;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

//   rU=L->rU*L->lambda/D->lambda;
//   rD=L->rD*L->lambda/D->lambda;
//   flat=L->flat*L->lambda/D->lambda;
  // 20150220 mshur
   rU=L->rU;
   rD=L->rD;
   flat=L->flat;
  //----------------
   fratio=L->omega/D->omega;
   phs=L->phs;  // 20150219 mshur

   if(t<=3*rU) {
      t0=3*rU; 
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0)+phs)*exp(-(t-t0)*(t-t0)/rU/rU);
   }
   else if(t>=3*rU && t<3*rU+flat)  {
      t0=3*rU;
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0)+phs);
   }
   else if(t>=3*rU+flat && t<3*rU+3*rD+flat) {
      t0=3*rU+flat;
      longitudinal=L->amplitude*sin(2.0*pi*fratio*(t-t0)+phs)*exp(-(t-t0)*(t-t0)/rD/rD);
   }
   else if(t>=3*rU+3*rD+flat) 
      longitudinal=0.0;


   if(L->loadPoint>=D->minXSub && L->loadPoint<=D->maxXSub)
      rank=myrank;
   else
      rank=nTasks+1;
   position=L->loadPoint-D->minXSub+1;

   if(myrank==rank)
   {
      if(L->polarity==1) { 
         D->field[position].E1+=longitudinal;
      }
      else if(L->polarity==2)
      {
         D->field[position].Pr+=longitudinal;            
         D->field[position].Pl+=longitudinal;
      }           
      else if(L->polarity==3)
      {
         D->field[position].Sl+=longitudinal;            
         D->field[position].Sr+=longitudinal; 
      }
      else if(L->polarity==4)	//right circular
      {
         D->field[position].Pr+=longitudinal;            
         D->field[position].Pl+=longitudinal; 
         D->field[position].Sl+=longitudinal/tan(2.0*pi*fratio*(t-t0));         
         D->field[position].Sr+=longitudinal/tan(2.0*pi*fratio*(t-t0)); 
      }
      else if(L->polarity==5)	//left circular
      {
         D->field[position].Pr+=longitudinal;            
         D->field[position].Pl+=longitudinal; 
         D->field[position].Sl+=-longitudinal/tan(2.0*pi*fratio*(t-t0));  
         D->field[position].Sr+=-longitudinal/tan(2.0*pi*fratio*(t-t0)); 
      }
   }
}
