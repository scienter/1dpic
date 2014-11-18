#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "laser.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void saveCurrent(Domain *D,int iteration)
{
    int i;
    char name[100];
    float x,J1,J2,J3,factor;
    FILE *out;
    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"current%d_%d",iteration,myrank);
    out = fopen(name,"w");

    factor=D->gamma*(1+D->beta);
    for(i=2; i<D->nxSub+2; i++)
    {
       x=(i-2+D->minXSub)*D->dx*D->lambda;
       J1=D->field[i].J1;    
       J2=D->field[i].J2;
       J3=D->field[i].J3;
       fprintf(out,"%g %g %g %g\n",x,J1,J2,J3);
    }             
    fclose(out);
}

void saveProbe(Domain *D,int iteration)
{
    int i,n;
    char name[100];
    float t,Ex,Ey,Ez,Bx,By,Bz,Pr,Pl,Sr,Sl,x;
    float px,p1,p2,p3,J1,J2,J3;
    float omega,frequency,dt;
    FILE *out;
    int myrank, nprocs; 
    Probe **probe;   
    probe=D->probe;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    omega=2*pi*velocityC/D->lambda;
    frequency=omega/2.0/pi;   
    dt=1.0/frequency*D->dt;

    for(n=0; n<D->probeNum; n++)
    { 
//      if(D->probeX[n]>=D->minXSub && D->probeX[n]<D->maxXSub)
//      {
        x=D->probeX[n]*D->dx*D->lambda;
        sprintf(name,"probeRaman%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Pr=probe[n][i].Pr;    
          Pl=probe[n][i].Pl;
          Sr=probe[n][i].Sr;
          Sl=probe[n][i].Sl;    
          fprintf(out,"%g %g %g %g %g %g\n",t,Pr,Pl,Sr,Sl,x);
        }             
        fclose(out);

        sprintf(name,"probe%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Ex=probe[n][i].E1;
          Bx=probe[n][i].B1;
          Ey=probe[n][i].Pr+probe[n][i].Pl;    
          Ez=probe[n][i].Sr+probe[n][i].Sl;
          By=probe[n][i].Sl-probe[n][i].Sr;
          Bz=probe[n][i].Pr-probe[n][i].Pl;    
          fprintf(out,"%g %g %g %g %g %g %g %g\n",t,Ex,Ey,Ez,Bx,By,Bz,x);
        }
        fclose(out);
             
        sprintf(name,"probeParticle%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          px=probe[n][i].x*D->dx*D->lambda;
          p1=probe[n][i].p1;
          p2=probe[n][i].p2;    
          p3=probe[n][i].p3;
          fprintf(out,"%g %g %g %g %g\n",t,px,p1,p2,p3);
        }             
        fclose(out);

        sprintf(name,"probeJ%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          J1=probe[n][i].J1;
          J2=probe[n][i].J2;    
          J3=probe[n][i].J3;
          fprintf(out,"%g %g %g %g\n",t,J1,J2,J3);
        }             
        fclose(out);
//      }
    }
}

void saveRho(Domain *D,int iteration)
{
    int i,s;
    char name[100];
    double x1,tmp;
    Particle *particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->charge*LL->density/((float)LL->numberInCell);
       LL=LL->next;
       s++;
    }

    //initializing density
    for(i=0; i<D->nxSub+3; i++)
       particle[i].rho=0;

    for(s=0; s<D->nSpecies; s++)
    {
       for(i=2; i<D->nxSub+2; i++)
       {
          p=particle[i].head[s]->pt;
          while(p)
          {
             x1=p->x1;
             particle[i].rho+=(1-x1)*rho0[s];
             particle[i+1].rho+=x1*rho0[s];
             p=p->next;
          }
       }

       sprintf(name,"rho%d_%d",iteration,myrank);
       out = fopen(name,"w");    
       for(i=2; i<D->nxSub+2; i++)
       {
          x1=(i-2+D->minXSub)*D->dx*D->lambda;
          fprintf(out,"%g %g\n",x1,particle[i].rho);    
       }           
       fclose(out);
    }	//End of for(s)
}

void boostSaveField(Domain *D,int labSaveStep)
{
    int i,j;
    char name[100];
    float x1,e1,pr,pl,sr,sl;
    float factor,dx;
    FILE *out;
    int myrank, nprocs;    
    FieldElement *field;
    Boost *boost;
    field=D->field;
    boost=D->boost;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"bField%d_%d",labSaveStep,myrank);
    out = fopen(name,"w");

    for(i=1; i<=D->nxSub; i++)
    {
       x1=boost[i].X1;
       e1=boost[i].E1;
       pr=boost[i].Pr;
       pl=boost[i].Pl;
       sr=boost[i].Sr;
       sl=boost[i].Sl;
       fprintf(out,"%g %g %g %g %g %g\n",x1,e1,pr,pl,sr,sl);
    }             
    fclose(out);
    
    if(myrank==0)
       printf("bField%d is saved.\n",labSaveStep);
}

void saveField(Domain *D,int iteration)
{
    int i;
    char name[100];
    float x,Ex,Ey,Ez,Bx,By,Bz,factor;
    FILE *out;
    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"field%d_%d",iteration,myrank);
    out = fopen(name,"w");

    factor=D->gamma*(1+D->beta);
    for(i=2; i<D->nxSub+2; i++)
    {
       x=(i-2+D->minXSub+0.5)*D->dx*D->lambda;
       Ex=D->field[i].E1/factor;    
       Ey=D->field[i].Pr+D->field[i].Pl;
       Ez=D->field[i].Sr+D->field[i].Sl;
       Bx=0;
       By=D->field[i].Sl-D->field[i].Sr;
       Bz=D->field[i].Pr-D->field[i].Pl;       
       fprintf(out,"%g %g %g %g %g %g %g\n",x,Ex,Ey,Ez,Bx,By,Bz);
    }             
    fclose(out);
}

void saveRaman(Domain *D,int iteration)
{
    int i;
    char name[100];
    float x,Pr,Pl,Sr,Sl;
    FILE *out;
    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"raman%d_%d",iteration,myrank);
    out = fopen(name,"w");
    for(i=2; i<D->nxSub+2; i++)
    {
       x=(i-2+D->minXSub+0.5)*D->dx*D->lambda;
       Pr=D->field[i].Pr;
       Pl=D->field[i].Pl;
       Sr=D->field[i].Sr;
       Sl=D->field[i].Sl;
       fprintf(out,"%g %g %g %g %g\n",x,Pr,Pl,Sr,Sl);
    }             
    fclose(out);
}

void saveParticle(Domain *D,int iteration)
{
    Particle *particle;
    particle=D->particle;
    int i,s;
    char name[100];
    float x1,p1,p2,p3,index,gamma,mc;
    ptclList *p;
    FILE *out;
    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for(s=0; s<D->nSpecies; s++)
    {
       sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
       out = fopen(name,"w");    
       for(i=2; i<D->nxSub+2; i++)
       {
          p=particle[i].head[s]->pt;
          while(p)
          {
             x1=((i-2+D->minXSub)+p->x1)*D->dx*D->lambda; 
             p1=p->p1;
             p2=p->p2;    
             p3=p->p3;
             gamma=sqrt(1+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
             index=p->index;
             fprintf(out,"%g %g %g %g %g %g\n",x1,p1,p2,p3,gamma,index);               
//             fprintf(out,"%g %g %g %g %g %g\n",x1,p->E1,p->Pr,p->Pl,p->Sr,p->Sl);               
             p=p->next;
          }
       }
       fclose(out);
    }	//End of for(s)
}
