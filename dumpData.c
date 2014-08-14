#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

// Temporary routine to dump and restore data
void saveDump1D(Domain D,int iteration)
{
   FILE *out;
   char name[100];
   int i,j,istart,iend,n,s,cnt;
   Particle *particle;
   particle=D.particle;
   Probe **probe;
   probe=D.probe;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   out = fopen(name, "w");   

   // Save simulation Domain information
   fwrite(&(D.nxSub),sizeof(int),1,out);
   fwrite(&(D.minXSub),sizeof(int),1,out);
   fwrite(&(D.maxXSub),sizeof(int),1,out);
   fwrite(&(D.probeNum),sizeof(int),1,out);

   // Save informations of particles inside the domain
   for(s=0; s<D.nSpecies; s++)
     for (i=2; i<D.nxSub+2; i++) 
     {
       p = particle[i].head[s]->pt;
       cnt = 0;
       while(p)  { 
         p=p->next; 
         cnt++;    
       }
       fwrite(&cnt,sizeof(int),1,out);

       p = particle[i].head[s]->pt;
       while(p)  
       { 
         fwrite(&(p->x1),sizeof(float),1,out);   
         fwrite(&(p->p1),sizeof(float),1,out);   
         fwrite(&(p->p2),sizeof(float),1,out);   
         fwrite(&(p->p3),sizeof(float),1,out);   
         fwrite(&(p->index),sizeof(float),1,out);   

         p = p->next; 
       }            
     }		//End of for(i)

     FieldElement *field;
     field=D.field;

     for(i=0; i<D.nxSub+5; i++) 
     {
       fwrite(&(field[i].E1),sizeof(float),1,out);   
       fwrite(&(field[i].Pr),sizeof(float),1,out);   
       fwrite(&(field[i].Pl),sizeof(float),1,out);   
       fwrite(&(field[i].Sr),sizeof(float),1,out);   
       fwrite(&(field[i].Sl),sizeof(float),1,out);   
       fwrite(&(field[i].J1),sizeof(float),1,out);   
       fwrite(&(field[i].J2),sizeof(float),1,out);   
       fwrite(&(field[i].J3),sizeof(float),1,out);   
     }
  
   //save Probe data
   if(D.probeNum>0)
   {
     for(n=0; n<D.probeNum; n++)
       for(i=0; i<=D.maxStep; i++)
       { 
         fwrite(&(probe[n][i].Pr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Pl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].E1),sizeof(float),1,out);   
         fwrite(&(probe[n][i].B1),sizeof(float),1,out);   

         fwrite(&(probe[n][i].x),sizeof(float),1,out);   
         fwrite(&(probe[n][i].p1),sizeof(float),1,out);   
         fwrite(&(probe[n][i].p2),sizeof(float),1,out);   
         fwrite(&(probe[n][i].p3),sizeof(float),1,out);   

         fwrite(&(probe[n][i].J1),sizeof(float),1,out);   
         fwrite(&(probe[n][i].J2),sizeof(float),1,out);   
         fwrite(&(probe[n][i].J3),sizeof(float),1,out);   
       }
   }

   fclose(out);
}


// Temporary routine to dump and restore data
void restoreData1D(Domain *D, int iteration)
{
   FILE *in;
   char name[100];
   int i,j,n,s,cnt;
   float tmp;
   double tmp1;
   ptclList *p;
   Particle *particle;
   particle=D->particle;
   Probe **probe;
   probe=D->probe;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   in = fopen(name, "r");   

   // restore simulation Domain information
   fread(&(D->nxSub),sizeof(int),1,in);
   fread(&(D->minXSub),sizeof(int),1,in);
   fread(&(D->maxXSub),sizeof(int),1,in);
   fread(&(D->probeNum),sizeof(int),1,in);

   // restore informations of particles inside the domain
   for(s=0; s<D->nSpecies; s++)
     for (i=2; i<D->nxSub+2; i++) 
       {
         fread(&cnt,sizeof(int),1,in);

         for(n=0; n<cnt; n++)  
         { 
           p = (ptclList *)malloc(sizeof(ptclList)); 
           p->next = particle[i].head[s]->pt;
           particle[i].head[s]->pt = p;

           fread(&(p->x1),sizeof(float),1,in);   
           fread(&(p->p1),sizeof(float),1,in);   
           fread(&(p->p2),sizeof(float),1,in);   
           fread(&(p->p3),sizeof(float),1,in);   
           fread(&(p->index),sizeof(float),1,in); 
         }
       }

     FieldElement *field;
     field=D->field;

     for(i=0; i<D->nxSub+5; i++) 
     {
       fread(&(field[i].E1),sizeof(float),1,in);
       fread(&(field[i].Pr),sizeof(float),1,in);
       fread(&(field[i].Pl),sizeof(float),1,in);
       fread(&(field[i].Sr),sizeof(float),1,in);
       fread(&(field[i].Sl),sizeof(float),1,in);
       fread(&(field[i].J1),sizeof(float),1,in);
       fread(&(field[i].J2),sizeof(float),1,in);
       fread(&(field[i].J3),sizeof(float),1,in);
     }

   //restore Probe data
   if(D->probeNum>0)
   {
     for(n=0; n<D->probeNum; n++)
       for(i=0; i<=D->maxStep; i++)
       { 
         fread(&(probe[n][i].Pr),sizeof(float),1,in);   
         fread(&(probe[n][i].Pl),sizeof(float),1,in);   
         fread(&(probe[n][i].Sr),sizeof(float),1,in);   
         fread(&(probe[n][i].Sl),sizeof(float),1,in);   
         fread(&(probe[n][i].E1),sizeof(float),1,in);   
         fread(&(probe[n][i].B1),sizeof(float),1,in);   

         fread(&(probe[n][i].x),sizeof(float),1,in);   
         fread(&(probe[n][i].p1),sizeof(float),1,in);   
         fread(&(probe[n][i].p2),sizeof(float),1,in);   
         fread(&(probe[n][i].p3),sizeof(float),1,in);   

         fread(&(probe[n][i].J1),sizeof(float),1,in);   
         fread(&(probe[n][i].J2),sizeof(float),1,in);   
         fread(&(probe[n][i].J3),sizeof(float),1,in);   
       }
   }
   
   fclose(in);
}

