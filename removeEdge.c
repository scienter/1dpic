#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"


void removeEdge(Domain *D)
{
    int i,s;
    Particle *particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    for(i=0; i<2; i++)
    {
      for(s=0; s<D->nSpecies; s++)
      {   
         p=particle[i].head[s]->pt;
         while(p)
         { 	
           tmp=p->next;
           particle[i].head[s]->pt=tmp; 
           p->next=NULL;
           free(p);
           p=particle[i].head[s]->pt;
         }
      }
    }

    i=D->nxSub+2;
    for(s=0; s<D->nSpecies; s++)
    {   
       p=particle[i].head[s]->pt;
       while(p)
       {	
          tmp=p->next;
          particle[i].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i].head[s]->pt;
       }
    }

}
