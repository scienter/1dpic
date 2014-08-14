#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"

void solveField(Domain *D)
{
    int i;  
    float dx,dt,nowPr,nowSr,prevPr,prevSr;
    FieldElement *field;
    field=D->field;

    dx=D->dx;
    dt=D->dt;
    nowPr=field[1].Pr;
    nowSr=field[1].Sr;
    for(i=2; i<D->nxSub+2; i++)
    {             

       field[i].E1=field[i].E1-2*pi*dt*field[i].J1;       
       prevPr=nowPr;
       nowPr=field[i].Pr;
       field[i].Pr=prevPr-pi*dt*field[i].J2;
       field[i-1].Pl=field[i].Pl-pi*dt*field[i].J2;
       prevSr=nowSr;
       nowSr=field[i].Sr;
       field[i].Sr=prevSr-pi*dt*field[i].J3;
       field[i-1].Sl=field[i].Sl-pi*dt*field[i].J3;
    }
        
}


