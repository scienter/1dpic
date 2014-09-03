#include <stdio.h>
#include "mesh.h"


void interpolation_2nd(Domain *D,External *Ext)
{
   int i,ii,s;
   float E1,Pr,Pl,Sr,Sl,alpha,Wx1,Wx2,Wx3;
   float E1_i,E1_ip,Pr_i,Pr_ip,Pl_i,Pl_ip,Sr_i,Sr_ip,Sl_i,Sl_ip;
   float extE1,extPr,extPl,extSr,extSl,extB1;
   ptclList *p;
   FieldElement *field;
   field=D->field;

   extE1=Ext->E1;   
   extB1=Ext->B1;   
   extPr=Ext->Pr;   
   extPl=Ext->Pl;   
   extSr=Ext->Sr;   
   extSl=Ext->Sl;
   
   for(i=2; i<D->nxSub+2; i++)
   {
      for(s=0; s<D->nSpecies; s++)
      {
         p=D->particle[i].head[s]->pt;
         while(p)
         {
            alpha=p->x1;
            ii=i;
            Wx1=0.5*(1-alpha)*(1-alpha);
            Wx2=0.75-(alpha-0.5)*(alpha-0.5);          
            Wx3=0.5*(alpha)*(alpha);          
            E1=Wx1*field[ii-1].E1+Wx2*field[ii].E1+Wx3*field[ii+1].E1;
            Pr=Wx1*field[ii-1].Pr+Wx2*field[ii].Pr+Wx3*field[ii+1].Pr;
            Pl=Wx1*field[ii-1].Pl+Wx2*field[ii].Pl+Wx3*field[ii+1].Pl;
            Sr=Wx1*field[ii-1].Sr+Wx2*field[ii].Sr+Wx3*field[ii+1].Sr;
            Sl=Wx1*field[ii-1].Sl+Wx2*field[ii].Sl+Wx3*field[ii+1].Sl;
            p->E1=E1+extE1;
            p->Pr=Pr+extPr;  p->Pl=Pl+extPl;
            p->Sr=Sr+extSr;  p->Sl=Sl+extSl;

            p=p->next;
         }
     }        
   }
}

void interpolation_1st(Domain *D,External *Ext)
{
   int i,ii,s;
   float E1,Pr,Pl,Sr,Sl,alpha;
   float E1_i,E1_ip,Pr_i,Pr_ip,Pl_i,Pl_ip,Sr_i,Sr_ip,Sl_i,Sl_ip;
   float extE1,extPr,extPl,extSr,extSl,extB1;
   ptclList *p;
   FieldElement *field;
   field=D->field;

   extE1=Ext->E1;   
   extB1=Ext->B1;   
   extPr=Ext->Pr;   
   extPl=Ext->Pl;   
   extSr=Ext->Sr;   
   extSl=Ext->Sl;
   
   for(i=2; i<D->nxSub+2; i++)
   {
      for(s=0; s<D->nSpecies; s++)
      {
         p=D->particle[i].head[s]->pt;
         while(p)
         {
            alpha=p->x1+0.5-((int)(p->x1+0.5));
            ii=(int)(p->x1+i+0.5)-1;
            E1=(1-alpha)*field[ii].E1+alpha*field[ii+1].E1;
            Pr=(1-alpha)*field[ii].Pr+alpha*field[ii+1].Pr;
            Pl=(1-alpha)*field[ii].Pl+alpha*field[ii+1].Pl;
            Sr=(1-alpha)*field[ii].Sr+alpha*field[ii+1].Sr;
            Sl=(1-alpha)*field[ii].Sl+alpha*field[ii+1].Sl;
            p->E1=E1+extE1;
            p->Pr=Pr+extPr;  p->Pl=Pl+extPl;
            p->Sr=Sr+extSr;  p->Sl=Sl+extSl;

            p=p->next;
         }
      }        
   }
}
