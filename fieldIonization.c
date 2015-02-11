#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void fieldIonization(Domain *D,Ionization *I)
{
    int i,m,s,startS,endS,numMat;
    float tmp,Efield;
    LoadList *LL;   
    Particle *particle;
    particle=D->particle;

    ptclList *p,*p0,*p1;
    int Z[D->nSpecies];
 printf("totalS=%d\n",D->nSpecies);   

    LL=D->loadList;
    s=0;
    while(LL->next)
    {
      Z[s]=LL->charge+1;
      LL=LL->next;
      s++;
    }
    for(s=0; s<D->nSpecies; s++)
      printf("Z[%d]=%d\n",s,Z[s]);

    numMat=I->numberOfMaterials;
    for(i=2; i<D->nxSub+2; i++)
    {
      for(m=0; m<numMat; m++)
      {
        startS=I->material[m][0];
        endS=I->material[m][1];
        for(s=startS; s<=endS; s++)
        {
          p=particle[i].head[s]->pt;
          while(p)
          {
            tmp=p->E1*p->E1;
            tmp+=(p->Pr+p->Pl)*(p->Pr+p->Pl);
            tmp+=(p->Sr+p->Sl)*(p->Sr+p->Sl);
            Efield=sqrt(tmp);

            p=p->next; 
          }      
        }
      }		//End of for(s)
    }
}

