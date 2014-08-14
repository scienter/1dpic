#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include <math.h>

void particlePush(Domain *D)
{
    int i,l,m,s,shift,cnt;
    float x,shiftX=0,dt,dx,gamma,sqrT,coef;
    float pMinus[3],T[3],S[3],operate[3][3],pPlus[3];
    float E1,B1,Pr,Pl,Sr,Sl;
    Particle *particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;
    void deleteNode();
 
    for(i=0; i<3; i++)   {
       pMinus[i]=0.0;
       T[i]=0.0;
       S[i]=0.0;
       pPlus[i]=0.0;
    }

    dt=D->dt;
    dx=D->dx;
    B1=0.0;

    float mass[D->nSpecies];
    int charge[D->nSpecies];
    LL=D->loadList;
    s=0;
    while(LL->next)
    {
       mass[s]=LL->mass;
       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    for(i=2; i<D->nxSub+2; i++)
    {
       for(s=0; s<D->nSpecies; s++)
       {
          p=particle[i].head[s]->pt;     
          while(p)
          {    
             coef=pi*charge[s]/mass[s]*dt;

             //Calculate vector P- 
             pMinus[0]=p->p1+coef*(p->E1);
             pMinus[1]=p->p2+coef*(p->Pr+p->Pl);    
             pMinus[2]=p->p3+coef*(p->Sr+p->Sl);
 
             //Calculate vector T 
             gamma=sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
             T[0]=coef/gamma*B1;   //Bx=0;
             T[1]=coef/gamma*(p->Sl-p->Sr);
             T[2]=coef/gamma*(p->Pr-p->Pl);

             //Calculate vector S
             sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
             for(l=0; l<3; l++)  
                S[l]=2.0*T[l]/sqrT;
  
             //Calculate operator A from P+=A.P-
             operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
             operate[0][1]=S[1]*T[0]+S[2];    
             operate[0][2]=S[2]*T[0]-S[1];     
             operate[1][0]=S[0]*T[1]-S[2];        
             operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
             operate[1][2]=S[2]*T[1]+S[0];         
             operate[2][0]=S[0]*T[2]+S[1];    
             operate[2][1]=S[1]*T[2]-S[0];    
             operate[2][2]=1-S[0]*T[0]-S[1]*T[1]; 
             //Calculate vector P+
             for(l=0; l<3; l++)  {
                pPlus[l]=0.0;
                for(m=0; m<3; m++)   
                   pPlus[l]+=operate[l][m]*pMinus[m];
                }
             //Updated momentum              
             p->p1=pPlus[0]+coef*(p->E1); 
             p->p2=pPlus[1]+coef*(p->Pr+p->Pl);    
             p->p3=pPlus[2]+coef*(p->Sr+p->Sl); 
    
             //Translation
             gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
             shiftX=p->p1/gamma;    //*dt is ignored because of dx=dt=1 in cell.
             if(shiftX>=1)  {
                printf("particle's movement exceeds C velocity");
                exit(0);
             } 
             p->oldX1=i+p->x1;
             p->x1+=shiftX;
           
             p=p->next;
          }	//End of while(p)
       }	//End of for(s)
    }      	//End of for(i)

}

             
