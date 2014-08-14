#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>

void loadMovingPlasma_random(Domain *D)
{
   int i,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int position;
   float tmp,dx,cx,mass;
   float wp,pDt,v1,v2,v3,gamma;

   float ne,randTest=0,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();

   LoadList *LL;
   ptclList *New,*p;   
   Particle *particle;
   particle=D->particle;

   dx=D->dx;

   i=D->nxSub+1;
   //position define      
     LL=D->loadList;
     s=0;
     while(LL->next)
     {
       for(l=0; l<LL->lnodes-1; l++)
       {
         position=i+D->minXSub;
 
         if(position>=LL->lpoint[l] && position<LL->lpoint[l+1])
         {
           ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
              *(position-LL->lpoint[l])+LL->ln[l]);
           ne*=LL->numberInCell;	//it is the float number of superparticles.
           intNum=(int)ne;
           randTest=ne-intNum;
           cnt=0;
           while(cnt<intNum)
           {               
             positionX=randomValue((int)LL->numberInCell);

             if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
             {
                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s]->pt;
                particle[i].head[s]->pt = New;
                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            
             } 
             else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
             {
                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s]->pt;
                particle[i].head[s]->pt = New;

                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            

                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s+1]->pt;
                particle[i].head[s+1]->pt = New;

                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            
             } 
             
             cnt++;
           }		//end of while(cnt)

           if(randTest>randomValue((int)LL->numberInCell))
           {
             positionX=randomValue((int)LL->numberInCell);

             if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
             {
                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s]->pt;
                particle[i].head[s]->pt = New;

                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            
             } 
             else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
             {
                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s]->pt;
                particle[i].head[s]->pt = New;

                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            

                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i].head[s+1]->pt;
                particle[i].head[s+1]->pt = New;

                New->x1 = positionX;
                New->oldX1=i+positionX;
                LL->index++;
                New->index=LL->index;            
             } 
           }		//end of if(randTest)

         }
       } 			//end of for(lnodes)  

       LL=LL->next;
       s++;
     }				//End of while(LL)   


   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
      mass=LL->mass*eMass;
         p=particle[i].head[s]->pt;
         while(p)
         {
           p->E1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v2=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v3=maxwellianVelocity(LL->temperature,mass)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+gamma*v1;
           p->p2=gamma*v2;
           p->p3=gamma*v3;

           p=p->next;
         } 
       LL=LL->next;
       s++;
   }					//End of while(LL)

}


void loadPlasma_random(Domain *D)
{
   int i;
   int s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int position;
   float cx,tmp,dx,mass;
   float wp,pDt,v1,v2,v3,gamma;
   float ne,randTest,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();
   Particle *particle;
   particle=D->particle;

   LoadList *LL;
   ptclList *New,*p;   

   dx=D->dx;


   //position define      
   for(i=2; i<D->nxSub+2; i++)
   {
       LL=D->loadList;
       s=0;
       while(LL->next)
       {
         for(l=0; l<LL->lnodes-1; l++)
         {
           position=i+D->minXSub;
 
           if(position>=LL->lpoint[l] && position<LL->lpoint[l+1])
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])*(position-LL->lpoint[l])+LL->ln[l]);
             ne*=LL->numberInCell;	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue((int)LL->numberInCell);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;

                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;
 
                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s+1]->pt;
                  particle[i].head[s+1]->pt = New;

                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            
               }              
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue((int)LL->numberInCell))
             {
               positionX=randomValue((int)LL->numberInCell);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;

                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s]->pt;
                  particle[i].head[s]->pt = New;
 
                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = particle[i].head[s+1]->pt;
                  particle[i].head[s+1]->pt = New;

                  New->x1 = positionX;
                  New->oldX1=i+positionX;
                  LL->index++;
                  New->index=LL->index;            
               } 
             }		//end of if(randTest)
           }
         } 			//end of for(lnodes)  

         LL=LL->next;
         s++;
       }				//End of while(LL)   
     }				//End of for(i,j)
    
   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
     mass=LL->mass*eMass;

     for(i=2; i<D->nxSub+2; i++)
     {
         p=particle[i].head[s]->pt;
         while(p)
         {
           p->E1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v2=maxwellianVelocity(LL->temperature,mass)/velocityC;
           v3=maxwellianVelocity(LL->temperature,mass)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+gamma*v1;
           p->p2=gamma*v2;
           p->p3=gamma*v3;

           p=p->next;
         } 
       }				//end of for(i,j)
       LL=LL->next;
       s++;
   }					//End of while(LL)
}

double maxwellianVelocity(double temperature,double mass)
{
   float vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/mass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((float)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((float)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

double randomValue(int numberInCell)
{
   double r;
   int intRand, randRange=numberInCell*100;

   intRand = rand() % randRange;
   r = ((float)intRand)/randRange;

   return r;
}
