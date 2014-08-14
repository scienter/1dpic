#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include "laser.h"

void clean1D(Domain *D)
{
//    FieldElement *field;
//    field=D->field;
    int i,s;
    ptclList *p,*tmp;
    LoadList *LL,*tmpLL;
    LaserList *L, *tmpL;

    //remove particles
    for(i=0; i<D->nxSub+3; i++)
    {
      for(s=0; s<D->nSpecies; s++)
      {
        p=D->particle[i].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          D->particle[i].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=D->particle[i].head[s]->pt;
        }
        free(D->particle[i].head[s]);
      }
      free(D->particle[i].head);
    }		
    //remove field
    free(D->particle);
    free(D->field);
      
    LL=D->loadList;
    while(LL)
    {	
      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;
      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);

    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);

    //remove trans
//    free(D->btJ);
//    free(D->upJ);


/*
    //remove probe
    if(D->probeNum>0)
    {
      for(i=0; i<D->probeNum; i++)
       free(D->probe[i]);
      free(D->probe);
      free(D->probeX);
    }
*/
}

