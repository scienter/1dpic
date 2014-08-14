#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include <math.h>
#include "mpi.h"


void updateCurrent_3rd(Domain *D)
{
    int i,ii,s,n,i1,i2,intXc;
    int nxSub;
    double x,inverDt,xold,xnew,xr,xxc;
    double xx1,xx2,xx3,tmp;
    double Fx1,Fx2,Fx3,Wx1,Wx2,Wx3,dx,dy,dt;
    double vz,vy,gamma,xc,wx,xcc;
    ptclList *p;
    FieldElement *field; 
    LoadList *LL;   
    field=D->field;

    double maximum();
    double minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    nxSub=D->nxSub;
   
    dt=D->dt;
    dx=D->dx; 
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
//       coeff[s]=LL->charge*LL->superP/LL->criticalDensity/D->lambda/D->lambda/dx/dy;
//       charge[s]=LL->charge;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
    {
      field[i].J1=0.0;
      field[i].J2=0.0;
      field[i].J3=0.0;
    }

    for(i=2; i<nxSub+2; i++)
    {
      for(s=0; s<D->nSpecies; s++)
      {
        p=D->particle[i].head[s]->pt;     
         
        while(p) 
        {
           gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
           vy=p->p2/gamma;
           vz=p->p3/gamma;
           xnew=p->x1+i;
           xold=p->oldX1; 

           xc=0.5*(xold+xnew); 
           intXc=(int)xc; 
           xcc=xc-intXc; 

           i1=(int)(xold);
           i2=(int)(xnew);
           if(i1==i2) 
             xr=0.5*(xold+xnew);
           else 
//             xr=(i1+i2)*0.5;
             xr=maximum(i1*1.0,i2*1.0);

           tmp=0.5*(xold+xr);                       
           x=tmp-(int)tmp;
           xx1=0.5+x;
           xx2=0.5-x;
           xx3=1.5-x;
           Fx1=(xr-xold)*0.5*(xx1-1.5)*(xx1-1.5);
           Fx2=(xr-xold)*(0.75-xx2*xx2);
           Fx3=(xr-xold)*0.5*(1.5-xx3)*(1.5-xx3);
           field[i1-1].J1+=Fx1*coeff[s];
           field[i1+0].J1+=Fx2*coeff[s];
           field[i1+1].J1+=Fx3*coeff[s];

           tmp=0.5*(xnew+xr);                       
           x=tmp-(int)tmp;
           xx1=0.5+x;
           xx2=0.5-x;
           xx3=1.5-x;
           Fx1=(xnew-xr)*0.5*(xx1-1.5)*(xx1-1.5);
           Fx2=(xnew-xr)*(0.75-xx2*xx2);
           Fx3=(xnew-xr)*0.5*(1.5-xx3)*(1.5-xx3);
           field[i2-1].J1+=Fx1*coeff[s];
           field[i2+0].J1+=Fx2*coeff[s];
           field[i2+1].J1+=Fx3*coeff[s];

           intXc=(int)(xc+0.5);
           xxc=xc-intXc;
           Wx1=0.5*(xxc-0.5)*(xxc-0.5);      
           Wx2=0.75-xxc*xxc;      
           Wx3=0.5*(0.5+xxc)*(0.5+xxc);      
           field[intXc-1].J2+=Wx1*vy*coeff[s];
           field[intXc+0].J2+=Wx2*vy*coeff[s];
           field[intXc+1].J2+=Wx3*vy*coeff[s];
           field[intXc-1].J3+=Wx1*vz*coeff[s];
           field[intXc+0].J3+=Wx2*vz*coeff[s];
           field[intXc+1].J3+=Wx3*vz*coeff[s];

           p=p->next;
        }	//End of while(p)
      }	   //End of for(s)     
    }		//End of for(i,j)


}

void updateCurrent_2nd(Domain *D)
{
    int i,s,n,i1,i2,intXc;
    int nxSub;
    double inverDt,xold,xnew,xr;
    double Fx1,Fx2,Wx1,Wx2,Wx3,dx,dy,dt;
    double vz,vy,gamma,xc,wx,xcc,xxc;
    ptclList *p;
    FieldElement *field; 
    LoadList *LL;   
    field=D->field;

    double maximum();
    double minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    nxSub=D->nxSub;
   
    dt=D->dt;
    dx=D->dx; 
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
//       coeff[s]=LL->charge*LL->superP/LL->criticalDensity/D->lambda/D->lambda/dx/dy;
//       charge[s]=LL->charge;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
    {
      field[i].J1=0.0;
      field[i].J2=0.0;
      field[i].J3=0.0;
    }

    for(i=2; i<nxSub+2; i++)
    {
      for(s=0; s<D->nSpecies; s++)
      {
        p=D->particle[i].head[s]->pt;     
         
        while(p) 
        {
           gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
           vy=p->p2/gamma;
           vz=p->p3/gamma;
           xnew=p->x1+i;
           xold=p->oldX1; 

           xc=0.5*(xold+xnew); 
           intXc=(int)xc; 
           xcc=xc-intXc; 

           i1=(int)(xold+0.5);
           i2=(int)(xnew+0.5);
           if(i1==i2) 
             xr=0.5*(xold+xnew);
           else 
             xr=(i1+i2)*0.5;
//             xr=maximum(i1*1.0,i2*1.0);

           Wx1=0.5*(xold+xr)-i1;
           Fx1=(xr-xold)*(0.5-Wx1);
           Fx2=(xr-xold)*(0.5+Wx1);
           field[i1-1].J1+=Fx1*coeff[s];
           field[i1+0].J1+=Fx2*coeff[s];

           i1=(int)(xr+0.5);
           Wx1=0.5*(xnew+xr)-i1;
           Fx1=(xnew-xr)*(0.5-Wx1);
           Fx2=(xnew-xr)*(0.5+Wx1);
           field[i1-1].J1+=Fx1*coeff[s];
           field[i1+0].J1+=Fx2*coeff[s];

           intXc=(int)(xc+0.5);
           xxc=xc-intXc;
           Wx1=0.5*(0.5-xxc)*(0.5-xxc);      
           Wx2=0.75-xxc*xxc;      
           Wx3=0.5*(0.5+xxc)*(0.5+xxc);      
           field[intXc-1].J2+=Wx1*vy*coeff[s];
           field[intXc+0].J2+=Wx2*vy*coeff[s];
           field[intXc+1].J2+=Wx3*vy*coeff[s];
           field[intXc-1].J3+=Wx1*vz*coeff[s];
           field[intXc+0].J3+=Wx2*vz*coeff[s];
           field[intXc+1].J3+=Wx3*vz*coeff[s];

           p=p->next;
        }	//End of while(p)
      }	   //End of for(s)     
    }		//End of for(i,j)


}

void updateCurrent_1st(Domain *D)
{
    int i,s,n,i1,i2,intXc;
    int nxSub;
    double inverDt,x1,x2,xr;
    double Fx1,Fx2,Wx1,Wx2,dx,dy,dt;
    double vz,vy,gamma,xc,wx,xcc;
    ptclList *p;
    FieldElement *field; 
    LoadList *LL;   
    field=D->field;

    double maximum();
    double minimum();

    double coeff[D->nSpecies];
    int charge[D->nSpecies];

    nxSub=D->nxSub;
   
    dt=D->dt;
    dx=D->dx; 
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       coeff[s]=LL->charge*LL->density/LL->criticalDensity/LL->numberInCell;
//       coeff[s]=LL->charge*LL->superP/LL->criticalDensity/D->lambda/D->lambda/dx/dy;
//       charge[s]=LL->charge;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //initialize J
    for(i=0; i<nxSub+5; i++)
    {
      field[i].J1=0.0;
      field[i].J2=0.0;
      field[i].J3=0.0;
    }

    for(i=2; i<nxSub+2; i++)
    {
      for(s=0; s<D->nSpecies; s++)
      {
        p=D->particle[i].head[s]->pt;     
         
        while(p) 
        {
           gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
           vy=p->p2/gamma;
           vz=p->p3/gamma;
           x2=p->x1+i;
           x1=p->oldX1; 

           xc=0.5*(x1+x2); 
           intXc=(int)xc; 
           xcc=xc-intXc; 
           i1=(int)x1;
           i2=(int)x2;
           if(i1==i2) 
             xr=0.5*(x1+x2);
           else 
             xr=maximum(i1*1.0,i2*1.0);

           Fx1=(xr-x1);
           Fx2=(x2-xr);

           field[i1].J1+=Fx1*coeff[s];
           field[i2].J1+=Fx2*coeff[s];

           wx=1-xcc;            
           field[intXc].J2+=wx*vy*coeff[s];
           field[intXc].J3+=wx*vz*coeff[s];
           wx=xcc;            
           field[intXc+1].J2+=wx*vy*coeff[s];
           field[intXc+1].J3+=wx*vz*coeff[s];

           p=p->next;
        }	//End of while(p)
      }	   //End of for(s)     
    }		//End of for(i,j)


}


double findR(double x1, double x2,double x3, double x4)
{
  double minimum();
  double maximum();
  double result,result1,result2,result3;

  result1=minimum(x1-0.5,x2-0.5);
  result2=maximum(x1-1.5,x2-1.5);
  result3=maximum(result2,(x3+x4)*0.5);
  result=minimum(result1,result3);

  return result;
}

double maximum(double x1,double x2)
{
   double result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

double minimum(double x1,double x2)
{
   double result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}

int intmaximum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

int intminimum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}
