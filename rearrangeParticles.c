#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void rearrangeParticles(Domain *D)
{
    Particle *particle;
    particle=D->particle;

    int i,s,intX,cnt,deleteFlag=0;
    float x;
    ptclList *p,*New,*prev,*tmp;

    for(i=2; i<D->nxSub+2; i++)
      for(s=0; s<D->nSpecies; s++)
      {
        cnt=1;
        p=particle[i].head[s]->pt;
        while(p)
        {
          x=p->x1;
          if(x>=1.0)  {
            intX=(int)x;
            x-=intX;
            deleteFlag=1;
          }
          else if(x<0) {              
            intX=(int)(x-1);
            x-=intX;
            deleteFlag=1;
          } 
          else   intX=0;

          if(deleteFlag==1)
          {
            New = (ptclList *)malloc(sizeof(ptclList)); 
            New->next = particle[i+intX].head[s]->pt;
            particle[i+intX].head[s]->pt = New;
            New->x1=x;    New->oldX1=p->oldX1;
            New->p1=p->p1;  New->p2=p->p2;  New->p3=p->p3;
            New->index=p->index;

            if(cnt==1)
            {
              tmp=p->next;
              particle[i].head[s]->pt=tmp; 
              p->next=NULL;
              free(p);
              p=particle[i].head[s]->pt; 
              cnt=1;
            }
            else
            {
              prev->next=p->next;
              p->next=NULL;
              free(p);
              p=prev->next; 
              cnt++;
            }
          }		//End of if(deleteFlag==1)
          else
          {
            prev=p;
            p=p->next;
            cnt++;
          }

          deleteFlag=0;
        }
      }		//End of for(s)
}

