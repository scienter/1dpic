#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "laser.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>

void loadMovingPlasma_crystal(Domain *D)
{
   int i,s,l,intNum,cnt,np;
   float space,position,positionX,x,ne,minX;
   float leftIndex,rightIndex,nc;
   float v1,v2,v3,gamma,mass;
   double maxwellianVelocity();
   Particle *particle;
   LoadList *LL;
   particle=D->particle;
   ptclList *New,*p;         

   i=D->nxSub+1;
   LL=D->loadList;
   s=0;
   while(LL->next)
   {  
     mass=LL->mass*eMass;
     for(l=0; l<LL->lnodes-1; l++)
     {
       position=i+D->minXSub;
 
       if(position>=LL->lpoint[l] && position<LL->lpoint[l+1])
       {
         nc=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
            *(position-LL->lpoint[l])+LL->ln[l]);
         nc*=LL->numberInCell;	//it is the float number of superparticles.
         space=1.0/nc;
         p=particle[D->nxSub+1].head[s]->pt;
         x=-2;
         while(p)
         {
           if(p->x1>x)  x=p->x1;
           p=p->next;
         }
   
         while(x<1)
         {
           x+=space;
           positionX=x;

           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i].head[s]->pt;
           particle[i].head[s]->pt = New;

           New->x1 = positionX;   New->oldX1=i+positionX;
           New->E1=New->Pr=New->Pl=New->Sr=New->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v2=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v3=maxwellianVelocity(LL->temperature,mass)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           New->p1=-D->gamma*D->beta+gamma*v1;
           New->p2=gamma*v2;
           New->p3=gamma*v3;
           LL->index++;
           New->index=LL->index;            
         }
       }
     }
     LL=LL->next;
     s++;
   }		//End of while(LL)  
}

void loadPlasma_crystal(Domain *D)
{
   int i,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   float space,positionX,x,n0,n1,nc1,xL;
   float wp,pDt,v1,v2,v3,gamma,mass,refX;
   double maxwellianVelocity();
   LoadList *LL;
   Particle *particle;
   particle=D->particle;
   ptclList *New,*p;   

   LL=D->loadList;
   s=0;
   while(LL->next)
   {
      mass=LL->mass*eMass;
      for(l=0; l<LL->lnodes-1; l++)
      {

         if(LL->ln[l+1]-LL->ln[l]>0 || LL->ln[l+1]-LL->ln[l]<0)
         {
            n0=LL->ln[l];
            n1=LL->ln[l+1];
            nc=LL->numberInCell;
            xL=LL->lpoint[l+1]-LL->lpoint[l];

            np=((int)((n1+n0)*xL*0.5*nc));
            cnt=0;
            x=(sqrt(n0*n0*(np-1)*(np-1)+(n1+n0)*(n1-n0)*cnt*(np-1))-n0*(np-1))/(n1-n0)/(np-1)*xL+LL->lpoint[l]-D->minXSub;
            i=(int)(x);
            refX=0;
               p = particle[i].head[s]->pt;
               while(p)  {
                  if(refX<=p->x1)
                     refX=p->x1;
                  p=p->next;
               }
            space=1.0/LL->numberInCell*LL->ln[l];
            refX+=space;
            refX-=(int)refX;              

            while(cnt<np)
            {
               x=(sqrt(n0*n0*(np-1)*(np-1)+(n1+n0)*(n1-n0)*cnt*(np-1))-n0*(np-1))/(n1-n0)/(np-1)*xL+LL->lpoint[l]-D->minXSub+refX;
               i=(int)(x);
               if(i>=2 && i<=D->nxSub+1)
               {
                  positionX=x-i;
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;

                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  New->E1=New->Pr=New->Pl=New->Sr=New->Sl=0.0;
                  v1=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  v2=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  v3=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
                  New->p1=-D->gamma*D->beta+gamma*v1;
                  New->p2=gamma*v2;
                  New->p3=gamma*v3;
                  LL->index++;
                  New->index=LL->index;            
               }
               cnt++;
            }
         }

         else
         {
            nc=LL->numberInCell*LL->ln[l];
            space=1.0/((float)nc);
            leftIndex=(int)(LL->lpoint[l]-D->minXSub);
            rightIndex=(int)(LL->lpoint[l+1]-D->minXSub);
            if(rightIndex>D->nxSub+1)
               rightIndex=D->nxSub+2; 
            np=(int)((rightIndex-leftIndex)*nc);
            cnt=0;
            refX=0;
            x=leftIndex+space*(cnt+0.5);
            i=(int)x-1;
              p = particle[i].head[s]->pt;
              while(p)  {
                if(refX<=p->x1)
                  refX=p->x1;
                 p=p->next;
              }
            refX+=space*1.5;
            refX=fabs(1-refX);
            if(refX==space) { 
               cnt=1;
            }

            while(cnt<np)
            {
               x=leftIndex+space*(cnt+0.5)+refX;
               i=((int)(x));
               if(i>=2 && i<=D->nxSub+1)
               {
                  positionX=x-i;
                  
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;

                  New->x1 = positionX;   New->oldX1=i+positionX;
                  New->E1=New->Pr=New->Pl=New->Sr=New->Sl=0.0;
                  v1=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  v2=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  v3=maxwellianVelocity(LL->temperature,mass)/velocityC;
                  gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
                  New->p1=-D->gamma*D->beta+gamma*v1;
                  New->p2=gamma*v2;
                  New->p3=gamma*v3;
                  LL->index++;
                  New->index=LL->index;            
               }
               cnt++;
            }
         }

      }

      LL=LL->next;
      s++;
   }		//End of while(LL)   

}

