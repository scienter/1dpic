#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void fieldIonization(Domain *D,Ionization *I)
{
    int i,m,s,startS,endS,numMat,intRand,randRange=1e5,cnt;
    float tmp,Efield,constRate,constE,possibility,rate,random;
    LoadList *LL;   
    float invE; // 20150219 mshur
    Particle *particle;
    particle=D->particle;

    ptclList *p,*p0,*p1,*tmpP,*prev;
    int index[D->nSpecies];
    float Z[D->nSpecies], effN[D->nSpecies], ioniEnergy[D->nSpecies];
//    constRate=4.134e16/(D->omega*0.5/pi);
    constRate=1.52e15/(D->omega*0.5/pi);  // 20150219 mshur
//    constE=ioniEfield*eCharge/eMass/D->omega/velocityC;
    constE=eMass*D->omega*velocityC/eCharge;  // 20150219 mshur

    LL=D->loadList;
    s=0;
    while(LL->next)
    {
      Z[s]=LL->charge+1;
      effN[s]=Z[s]/sqrt(LL->ioniEnergy);
      //ioniEnergy[s]=LL->ioniEnergy;
      ioniEnergy[s]=LL->ioniEnergy*13.6;  // 20150219,mshur, [eV] unit
      index[s]=LL->index;
      LL=LL->next;
      s++;
    }

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
          cnt=1;
          while(p)
          {
            intRand=rand()%randRange;
            random=((float)intRand)/((float)randRange);

            tmp=0; //p->E1*p->E1;
            tmp+=(p->Pr+p->Pl)*(p->Pr+p->Pl);
            tmp+=(p->Sr+p->Sl)*(p->Sr+p->Sl);
/*
            Efield=sqrt(tmp)/constE;
            if(Efield>0)
              rate=constRate/8.0/pi/Z[s]*Efield*pow(4*eulerNum*pow(Z[s],3.0)/pow(effN[s],4.0)/Efield,2*effN[s])*exp(-2.0/3.0/Efield*pow(ioniEnergy[s],1.5));
            else rate=0.0;
*/
          // Modified by mshur 20150219, formula from Bruhwiler's paper
            Efield=constE*sqrt(tmp)*1e-9; //make it [GV/m]
            if(Efield>1e-3) { // 1e-3 is just a arbitrarily small number
               invE=pow(ioniEnergy[s],1.5)/Efield;
               rate = constRate*pow(4.0,effN[s])*ioniEnergy[s]
                     /(effN[s]*tgamma(2.0*effN[s]))
                     *pow(20.5*invE,2.0*effN[s]-1.0)*exp(-6.83*invE);
            }
            else rate=0.0;
           //-------------------------------

            possibility=rate*D->dt;

            if(cnt==1)
            {
              if(possibility>random)
              {
                //generating electrons
                p0=(ptclList *)malloc(sizeof(ptclList)); 
                p0->next = particle[i].head[startS-1]->pt;
                particle[i].head[startS-1]->pt = p0;
                p0->x1=p->x1;
                p0->p1=p0->p2=p0->p3=0.0;
                p0->E1=p0->Pr=p0->Pl=p0->Sr=p0->Sl=0.0;
                index[startS-1]+=1;
                p0->index=index[startS-1];
  
                //generating ions
                p1=(ptclList *)malloc(sizeof(ptclList)); 
                p1->next = particle[i].head[s+1]->pt;
                particle[i].head[s+1]->pt = p1;
                p1->x1=p->x1;
                p1->p1=p1->p2=p1->p3=0.0;
                p1->E1=p1->Pr=p1->Pl=p1->Sr=p1->Sl=0.0;
                index[s+1]+=1;
                p1->index=index[s+1];
  
                //deleting present species
                tmpP=p->next;
                particle[i].head[s]->pt=tmpP;
                p->next=NULL;
                free(p);
                p=particle[i].head[s]->pt;
                cnt=1;
              }
              else
              {
                prev=p;
                p=p->next;
                cnt++;
              }
            }		//End of if(cnt==1)

            else
            {
              if(possibility>random)
              {
                //generating electrons
                p0=(ptclList *)malloc(sizeof(ptclList)); 
                p0->next = particle[i].head[startS-1]->pt;
                particle[i].head[startS-1]->pt = p0;
                p0->x1=p->x1;
                p0->p1=p0->p2=p0->p3=0.0;
                p0->E1=p0->Pr=p0->Pl=p0->Sr=p0->Sl=0.0;
                index[startS-1]+=1;
                p0->index=index[startS-1];
  
                //generating ions
                p1=(ptclList *)malloc(sizeof(ptclList)); 
                p1->next = particle[i].head[s+1]->pt;
                particle[i].head[s+1]->pt = p1;
                p1->x1=p->x1;
                p1->p1=p1->p2=p1->p3=0.0;
                p1->E1=p1->Pr=p1->Pl=p1->Sr=p1->Sl=0.0;
                index[s+1]+=1;
                p1->index=index[s+1];
  
                //deleting present species
                prev->next=p->next;
                p->next=NULL;
                free(p);
                p=prev->next;
                cnt++;
              }
              else
              {
                prev=p;
                p=p->next;
                cnt++;
              }
            }		//End of else

          }     	//End of while(p) 
        }		//End of for(s)
      }			//End of for(m)
    }

    //saving index information
    LL=D->loadList;
    s=0;
    while(LL->next)
    {
      LL->index=index[s];
      LL=LL->next;
      s++;
    }
}

