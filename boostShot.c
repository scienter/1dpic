#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "laser.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void boostShot(Domain *D,int iteration,int labSaveStep)
{
    int i,j,labStep;
    char name[100];
    float x1,e1,e2,e3,b1,b2,b3;
    float X1,E1,E2,E3,B1,B2,B3;
    float tmp,factor,factor1;
    int myrank, nprocs;    
    FieldElement *field;
    Boost *boost;
    field=D->field;
    boost=D->boost;

    factor=D->gamma*(1+D->beta);
    factor1=1+D->beta;
    i=(int)(labSaveStep/D->gamma/D->beta/factor-iteration/D->beta-D->minXSub+1);
    if(i>=1 && i<=D->nxSub)
    {       
       //boost frame data
       x1=(i-1+D->minXSub)*D->dx*D->lambda;
       e1=field[i].E1;
       e2=field[i].Pr+field[i].Pl;
       e3=field[i].Sr+field[i].Sl;
       b1=0;
       b2=field[i].Sl-field[i].Sr;
       b3=field[i].Pr-field[i].Pl;

       //lab frame data
       X1=D->gamma*(x1+D->beta*iteration*D->dt*D->lambda);
       E1=e1/factor;
       E2=(e2+D->beta*b3)/factor1;
       E3=(e3-D->beta*b2)/factor1;
       B2=(b2-D->beta*b3)/factor1;
       B3=(b3+D->beta*e2)/factor1;

       boost[i].X1=X1;   
       boost[i].E1=E1;
       boost[i].Pr=0.5*(E2+B3);
       boost[i].Pl=0.5*(E2-B3);
       boost[i].Sr=0.5*(E3-B2);
       boost[i].Sl=0.5*(E3+B2);
    }
}

