#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void probe(Domain *D,int iteration)
{
    int i,n,iter,s,ii,index;
    ptclList *p;

    iter=iteration;
    for(n=0; n<D->probeNum; n++)
    {
      i=D->probeX[n]-D->minXSub+2;

      D->probe[n][iter].E1=D->field[i].E1;    
      D->probe[n][iter].Pr=D->field[i].Pr;
      D->probe[n][iter].Pl=D->field[i].Pl;
      D->probe[n][iter].B1=0;
      D->probe[n][iter].Sr=D->field[i].Sr;
      D->probe[n][iter].Sl=D->field[i].Sl;      

      D->probe[n][iter].J1=D->field[i].J1;
      D->probe[n][iter].J2=D->field[i].J2;
      D->probe[n][iter].J3=D->field[i].J3;

      for(s=0; s<D->nSpecies; s++)
        for(ii=2; ii<D->nxSub+2; ii++)
        {
          p=D->particle[ii].head[s]->pt;
          while(p)
          {
            index=p->index;
            if(index==D->probeIndex[n])
            {
              D->probe[n][iter].x=p->x1+ii+D->minXSub;
              D->probe[n][iter].p1=p->p1;
              D->probe[n][iter].p2=p->p2;
              D->probe[n][iter].p3=p->p3;
            }
            p=p->next;
          }
        }     //End of ii
    }

}

void findProbeParticle(Domain *D)
{
    int i,n,s;
    ptclList *p;

    for(n=0; n<D->probeNum; n++)
    {
      i=D->probeX[n]-D->minXSub+2;
      for(s=0; s<D->nSpecies; s++)
      {
        p=D->particle[i].head[s]->pt;
        if(p)
          D->probeIndex[n]=p->index;
      }
    }
}

