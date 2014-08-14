#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void movingDomain(Domain *D)
{
    int i,s;
    FieldElement *field;
    field=D->field;
    ptclList *p,*tmp,*New;      
   
    for(i=1; i<D->nxSub+2; i++)
    {
       field[i].E1=field[i+1].E1;
       field[i].Pr=field[i+1].Pr;
       field[i].Pl=field[i+1].Pl;
       field[i].Sr=field[i+1].Sr;
       field[i].Sl=field[i+1].Sl;
       field[i].J1=field[i+1].J1;
       field[i].J2=field[i+1].J2;
       field[i].J3=field[i+1].J3;
    }

    for(i=2; i<D->nxSub+2; i++)
    {
       for(s=0; s<D->nSpecies; s++)
       {
          p=D->particle[i].head[s]->pt;
          while(p)
          {
             p->x1-=1.0;
             p->oldX1-=1.0;
             p=p->next;
          }
       } 	//End of for(s)          
    }		//End of for(i)

    D->minXSub+=1;
    D->maxXSub+=1;
}  
     
