#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"

void boundary(Domain *D,External *Ext)
{
     int i,j,s,rank,remain,tmp,sub;
     float min,max;
     int myrank, nTasks;
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     D->nxSub=D->nx/nTasks;
     sub=D->nxSub;
     remain=D->nx%nTasks;
     min=max=0;
     rank=0;
     while(rank<nTasks)
     {
        if(rank<remain)   tmp=sub+1;
        else              tmp=sub;
        min=max;
        max=min+tmp;
        if(myrank==rank)
        {
           D->minXSub=min+D->minXSub;
           D->maxXSub=max+D->minXSub;
           D->nxSub=tmp;
        }
        rank++;
     }

     D->field = (FieldElement *)malloc((D->nxSub+5)*sizeof(FieldElement ));     
     for(i=0; i<D->nxSub+5; i++)
     {
        D->field[i].E1=0.0;
        D->field[i].Pr=0.0;
        D->field[i].Pl=0.0;
        D->field[i].Sr=0.0;
        D->field[i].Sl=0.0;
        D->field[i].J1=0.0;        
        D->field[i].J2=0.0; 
        D->field[i].J3=0.0; 
     }

     D->particle = (Particle *)malloc((D->nxSub+3)*sizeof(Particle ));     
     for(i=0; i<D->nxSub+3; i++)
     {
        D->particle[i].rho=0.0;  
        D->particle[i].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
        for(s=0; s<D->nSpecies; s++)
        {
          D->particle[i].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
          D->particle[i].head[s]->pt = NULL;
        }
     }     

/*
     D->boost = (Boost *)malloc((D->nxSub+2)*sizeof(Boost ));
     for(i=0; i<D->nxSub+2; i++)
     {
        D->boost[i].E1=0.0;
        D->boost[i].Pr=0.0;
        D->boost[i].Pl=0.0;
        D->boost[i].Sr=0.0;
        D->boost[i].Sl=0.0;
     }
*/

     D->probe = (Probe **)malloc(D->probeNum*sizeof(Probe *));     
     for(i=0; i<D->probeNum; i++)
       D->probe[i] = (Probe *)malloc(D->maxStep*sizeof(Probe ));     
     for(i=0; i<D->probeNum; i++)
       for(j=0; j<D->maxStep; j++)
       {
         D->probe[i][j].E1=0;
         D->probe[i][j].Pr=0;
         D->probe[i][j].Pl=0;
         D->probe[i][j].B1=0;
         D->probe[i][j].Sr=0;
         D->probe[i][j].Sl=0;

         D->probe[i][j].x=0;
         D->probe[i][j].p1=0;
         D->probe[i][j].p2=0;
         D->probe[i][j].p3=0;

         D->probe[i][j].J1=0;
         D->probe[i][j].J2=0;
         D->probe[i][j].J3=0;
       }
      
}
                       
