#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "laser.h"
#include "plasma.h"
#include "mesh.h"
#include "constants.h"
#include <math.h>


void parameterSetting(Domain *D,External *Ext, char *input)
{
   int FindParameters();
   int findLoadParameters();
   int findLaserParameters();
   float minX, maxX, position,factor,pMinX,pMaxX,pPosition;
   float normalB,normalE,E1,E2,E3,B1,B2,B3;
   char str[100],name[100];
   int i,rank,minT,maxT,boostSaveStep,boostSaveStart,boostSaveEnd,tmp,probeNum,fail=0;
   float lambda;
   LoadList *LL, *New;
   LaserList *L, *LNew;

   if(FindParameters("Domain",1,"CurrentType",input,str)) D->currentType=atoi(str);
   else D->currentType=1;
   if(FindParameters("Domain",1,"InterpolationType",input,str)) D->interpolationType=atoi(str);
   else D->interpolationType=1;

   //Boost frame
   if(FindParameters("Domain",1,"boostGamma",input,str)) D->gamma=atof(str);
   else D->gamma=1;
   if(D->gamma>1)   D->boostOn=1;
   else             D->boostOn=0;
   if(FindParameters("Domain",1,"boostSaveStart",input,str)) boostSaveStart=atoi(str);
   else  {
      printf("Boost saveStart must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"boostSaveEnd",input,str)) boostSaveEnd=atoi(str);
   else  {
      printf("Boost saveEnd must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"boostSaveStep",input,str)) boostSaveStep=atoi(str);
   else  {
      printf("Boost saveStep must be defined.\n");
      fail=1;
   }
   D->beta=sqrt(1-1.0/D->gamma/D->gamma);

   //Domain parameter setting
   if(FindParameters("Domain",1,"maxStep",input,str)) D->maxStep=atoi(str);
   else  {
      printf("in [Domain], maxStep=?  [number of iteration]\n");
      fail=1;
   }

   if(FindParameters("Domain",1,"saveStep",input,str)) D->saveStep=atoi(str);
   else  {
      printf("in [Domain], saveStep=? [number of iteration]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"saveStart",input,str)) D->saveStart=atoi(str);
   else  {
      printf("in [Domain], saveStart=? [number of iteration]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"fieldSave",input,str)) D->fieldSave=atoi(str);
   else  D->fieldSave=1;
   if(FindParameters("Domain",1,"particleSave",input,str)) D->particleSave=atoi(str);
   else  D->particleSave=1;
   if(FindParameters("Domain",1,"rhoSave",input,str)) D->rhoSave=atoi(str);
   else  D->rhoSave=1;
   if(FindParameters("Domain",1,"currentSave",input,str)) D->currentSave=atoi(str);
   else  D->currentSave=1;
   if(FindParameters("Domain",1,"dumpSave",input,str)) D->dumpSave=atoi(str);
   else  {
      printf("in [Domain], dumpSave=?  [0:off, 1:on]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dumpStart",input,str)) D->dumpStart=atoi(str);
   else  D->dumpStart=D->saveStart;

   if(FindParameters("Domain",1,"minX",input,str)) 
   {
      minX=atof(str);
      minX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("in [Domain], minX=? [m] \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str)) 
   {
      maxX=atof(str);
      maxX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("in [Domain], maxX=? [m]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"moving",input,str)) D->moving=atoi(str);
   else  {
      printf("in [Domain], moving=? \n");
      printf("1:moving domain on,  0:moving domain off\n");
      fail=1;
   }   
   if(FindParameters("Domain",1,"lambda",input,str)) 
   {
      D->lambda=atof(str);
      D->lambda*=D->gamma*(1+D->beta);
   }
   else  {
      printf("in [Domain], lambda=? [m]\n");
      printf("basic parameter for dx, dt.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"divisionLambda",input,str)) 
      D->divisionLambda=atof(str);
   else  {
      printf("In [Domain], divisionLambda=? [number of devided wavelength]\n");
      fail=1;
   }

   //External field parameter setting
   if(FindParameters("External",1,"E1",input,str)) Ext->E1=atof(str);
   else  {
      printf("in [External], E1=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"E2",input,str)) E2=atof(str);
   else  {
      printf("in [External], E2=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"E3",input,str)) E3=atof(str);
   else  {
      printf("in [External], E3=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B1",input,str)) Ext->B1=atof(str);
   else  {
      printf("in [External], B1=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B2",input,str)) B2=atof(str);
   else  {
      printf("in [External], B2=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B3",input,str)) B3=atof(str);
   else  {
      printf("in [External], B3=? [Tesla]\n");
      fail=1;
   }


   //additional Domain parameters  
   D->nx=((int)((maxX-minX)/D->lambda*D->divisionLambda));
   D->omega=2*pi*velocityC/D->lambda;
   D->dx=1.0/D->divisionLambda;
   D->dt=D->dx;
   if(D->boostOn==1)   D->minXSub=-D->nx;
   else 	       D->minXSub=0;

   //additional Boost parameters
   factor=D->gamma*(1+D->beta);
   D->minT=(int)(boostSaveStart/factor/factor); 	//boost frame iteration
   D->maxT=(int)((boostSaveStart+D->beta*D->nx)/(1+D->beta)-factor*D->gamma*D->minT*D->beta);	//boost frame iteration
   D->boostSaveStep=boostSaveStep;


   //additional external field parameters
   normalB=eMass*D->omega/(-eCharge);
   normalE=normalB*velocityC;
   Ext->Pr=(E2/normalE+B3/normalB)*0.5;
   Ext->Pl=(E2/normalE-B3/normalB)*0.5;
   Ext->Sr=(E3/normalE-B2/normalB)*0.5;
   Ext->Sl=(E3/normalE+B2/normalB)*0.5;

   //Probe parameter
   if(FindParameters("Domain",1,"probeNum",input,str)) D->probeNum=atoi(str);
   else { 
     printf("in [Domain], probeNum=?\n"); 
     exit(0); 
   } 
   
   if(D->probeNum>0)
   {
     D->probeX = (int *)malloc(D->probeNum*sizeof(int));
     for(i=0; i<D->probeNum; i++)
     {
       sprintf(name,"probeX%d",i);
       if(FindParameters("Probe",1,name,input,str)) 
         D->probeX[i]=((int)(atof(str)/D->lambda/D->dx));   
       else { 
         printf("in [Probe], probeX%d=?\n",i);  fail=1; 
       } 
     }
     D->probeIndex = (int *)malloc(D->probeNum*sizeof(int));
   }

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input)) 
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

   //Plasma parameter setting
   if(FindParameters("Domain",1,"crystal",input,str)) 
      D->crystal=atoi(str);
   else  
      D->crystal=0;
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;

   while(findLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;
  
   if(fail==1)
      exit(0);
}

int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   float position;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"polarity",input,str)) polarity=atoi(str);
   else  polarity=0;

   if(polarity)
   {
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
     {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1+D->beta);
     }
     else  L->lambda=D->lambda;
  
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  {
        printf("in [Laser], flat=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPosition",input,str)) position=atof(str);
     else  {
        printf("in [Laser], loadPosition=?  [m]\n");
        fail=1;
     }

     //additional laser parameters
     L->polarity=polarity;
     L->omega=2*pi*velocityC/L->lambda;
     L->loadPoint=((int)(position/D->lambda/D->dx));   
     if(fail==1)
        exit(0);
   }
   return polarity;
}

int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   int whatSpecies();
   double whatMass();
   float pointPosition;
   int whatCharge();
   char name[100], str[100];
   int i,species,fail=0;

   if(FindParameters("Plasma",rank,"species",input,name)) 
      species = whatSpecies(name);
   else  species = 0;

   if(species)
   {
      LL->species=species;

      if(FindParameters("Plasma",rank,"density",input,str)) 
      {
         LL->density=atof(str);
         LL->density*=D->gamma;
      }
      else  {
         printf("in [Plasma], density=? [m-3]\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"numberInCell",input,str)) 
         LL->numberInCell=atof(str);
      else  {
         printf("in [Plasma], numberInCell=? \n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"with_next_species",input,str)) 
         LL->withNextSpcs=atoi(str);
      else  {
         printf("in [Plasma], with_next_species=?  (0:no 1:yes)\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"with_prev_species",input,str)) 
         LL->withPrevSpcs=atoi(str);
      else  {
         printf("in [Plasma], with_prev_species=?  (0:no 1:yes)\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"startIndex",input,str)) 
         LL->index=atoi(str);
      else  {
         printf("in [Plasma], startIndex=? \n");
         printf("It may be 0 as default.\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"Lnodes",input,str)) LL->lnodes=atoi(str);
      else  {
         printf("in [Plasma], Lnodes=? \n");
         printf("Each nodes indicates plasma density changing.\n");
         fail=1;
      }
   
      LL->lpoint = (float *)malloc(LL->lnodes*sizeof(float));
      LL->ln = (float *)malloc(LL->lnodes*sizeof(float));   
      for(i=0; i<LL->lnodes; i++)
      {
         sprintf(name,"X%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) 
            LL->lpoint[i] = 2 + atof(str)/D->lambda/D->dx*(1+D->beta);
         else 
         { printf("in [Plasma], X%d=?  [m]\n",i);  fail=1; } 

         sprintf(name,"Ln%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) 
            LL->ln[i] = atof(str);
         else 
         { printf("in [Plasma], Ln%d=? [0 <= Ln%d <= 1] \n",i,i);  fail=1; } 
      }

      if(FindParameters("Plasma",rank,"pointPosition",input,str))  {
         pointPosition=atof(str);
         LL->pointPosition=((int)(pointPosition/D->lambda/D->dx));
      }
      else   LL->pointPosition=0;	//no points load.  

      if(FindParameters("Plasma",rank,"temperature",input,str))  
         LL->temperature=atof(str);
      else   LL->temperature=0.0;	

      LL->mass=whatMass(species);
      LL->charge=whatCharge(species);
      LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
      LL->superP=LL->density*D->lambda*D->dx/LL->numberInCell;

      LL->cnt=1;
      
   }	//end of if(species)
   
   if(fail==1)
      exit(0);

   return species;
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)	  	return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)	  	return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)	  	return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)	  	return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)	  	return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)	  	return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)	  	return (12.0111-6*eMassU)/eMassU;
   else {  printf("Species' mass not defined\n");  exit(0);  }
}

int whatCharge(int species)
{
   if(species == Electron) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0) 	 	return 0;
   else if(species == CPlus1) 	 	return 1;
   else if(species == CPlus2) 	 	return 2;
   else if(species == CPlus3) 	 	return 3;
   else if(species == CPlus4) 	 	return 4;
   else if(species == CPlus5) 	 	return 5;
   else if(species == CPlus6) 	 	return 6;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

