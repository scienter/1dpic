#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

/*
void MPI_TransferP_Moving(Domain *D)
{
    int i,numP=0,cnt,numberData,s;
    int myrank, nTasks;
    particle=D->particle;     
    float *frontP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            


    //Even - odd
    if(myrank%2==0 && myrank!=0)
    {
       numP=0;
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[1].head[s]->pt;
         while(p)   {
           p=p->next;
           numP++;
         }
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1 && myrank!=nTasks-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    numberData=numP*6;
    frontP=(float *)malloc(numberData*sizeof(float ));             

    if(myrank%2==0 && myrank!=0)
    {
       p=particle[1].head[s]->pt;
       frontP[0]=p->x1;     
       frontP[1]=p->oldX1;  
       frontP[2]=p->p1;     
       frontP[3]=p->p2;     
       frontP[4]=p->p3;     
       frontP[5]=p->index;  
       tmp=p->next;
       p->next=NULL;
       free(p);
       particle[1].head[s]->pt=tmp; 
       MPI_Send(frontP,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
    }
    
       else if(myrank%2==1 && myrank!=D->L-1) 
       {    
          MPI_Recv(frontP,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = particle[D->nxSub+1].head[s]->pt;
          particle[D->nxSub+1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
       
    //Odd - evem
    if(myrank%2==1)
    {
       numP=0;
       p=particle[1].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=D->L-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==1)
       {
          p=particle[1].head[s]->pt;
          frontP[0]=p->x1;     
          frontP[1]=p->oldX1;  
          frontP[2]=p->p1;     
          frontP[3]=p->p2;     
          frontP[4]=p->p3;     
          frontP[5]=p->q;      
          frontP[6]=p->m;      
          frontP[7]=p->index;  
          frontP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          particle[1].head[s]->pt=tmp; 

          MPI_Send(frontP,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==0 && myrank!=D->L-1)  
       {    
          MPI_Recv(frontP,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = particle[D->nxSub+1].head[s]->pt;
          particle[D->nxSub+1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }	//End of for(s)    

  free(frontP);
}
*/

void MPI_TransferP_Xplus(Domain *D)
{
    int i,s,n,numP=0,cnt,numberData;
    int myrank, nTasks;
    Particle *particle;
    particle=D->particle;     
    float *behindP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    //Even -> odd
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       numP=0;
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[D->nxSub+2].head[s]->pt;
         while(p)   {
           p=p->next;
           numP++;
         }
       }      
       MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1) 
       MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    numberData=numP*7;
    behindP=(float *)malloc(numberData*sizeof(float ));             
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       n=0;
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[D->nxSub+2].head[s]->pt;
         while(p)
         {
           behindP[7*n+0]=p->x1;     
           behindP[7*n+1]=p->oldX1-D->nxSub;  
           behindP[7*n+2]=p->p1;     
           behindP[7*n+3]=p->p2;     
           behindP[7*n+4]=p->p3;     
           behindP[7*n+5]=p->index;  
           behindP[7*n+6]=(float)s;  
           p=p->next;
           n++;
         }
       }

       MPI_Send(behindP,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD); 
    }
 
    else if(myrank%2==1) 
    {    
      MPI_Recv(behindP,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
         s=(int)(behindP[n*7+6]);
         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[2].head[s]->pt;
         particle[2].head[s]->pt = New;      

         New->x1=behindP[n*7+0];
         New->oldX1=behindP[n*7+1];
         New->p1=behindP[n*7+2];
         New->p2=behindP[n*7+3];
         New->p3=behindP[n*7+4];
         New->index=behindP[n*7+5];
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(behindP);   

    //Odd -> evem
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       numP=0;
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[D->nxSub+2].head[s]->pt;
         while(p)   {
           p=p->next;
           numP++;
         }
       }      
       MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=0) 
       MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    numberData=numP*7;
    behindP=(float *)malloc(numberData*sizeof(float ));             

    if(myrank%2==1 && myrank!=nTasks-1)
    {
       n=0;
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[D->nxSub+2].head[s]->pt;
         while(p)
         {
           behindP[n*7+0]=p->x1;     
           behindP[n*7+1]=p->oldX1-D->nxSub;  
           behindP[n*7+2]=p->p1;     
           behindP[n*7+3]=p->p2;     
           behindP[n*7+4]=p->p3;     
           behindP[n*7+5]=p->index;  
           behindP[n*7+6]=(float)s;  
           p=p->next;
           n++;
         }
       }

       MPI_Send(behindP,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD); 
    }
    
    else if(myrank%2==0 && myrank!=0) 
    {    
      MPI_Recv(behindP,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
         s=(int)(behindP[n*7+6]);
         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[2].head[s]->pt;
         particle[2].head[s]->pt = New;      

         New->x1=behindP[n*7+0];
         New->oldX1=behindP[n*7+1];
         New->p1=behindP[n*7+2];
         New->p2=behindP[n*7+3];
         New->p3=behindP[n*7+4];
         New->index=behindP[n*7+5];
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(behindP);   

}

void MPI_TransferP_Xminus(Domain *D)
{
    int i,n,numP=0,cnt,numberData,s;
    int myrank, nTasks;
    Particle *particle;
    particle=D->particle;     
    float *frontP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

    //Even - odd
    if(myrank%2==0 && myrank!=0)
    {
       numP=0;
       for(s=0; s<D->nSpecies; s++)
       {
         for(i=0; i<2; i++)
         {
           p=particle[i].head[s]->pt;
           while(p)   {
             p=p->next;
             numP++;
           }
         }
       }      
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1 && myrank!=nTasks-1) 
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    numberData=numP*8;
    frontP=(float *)malloc(numberData*sizeof(float ));             

    if(myrank%2==0 && myrank!=0)
    {
       n=0;
       for(s=0; s<D->nSpecies; s++)
       {
         for(i=0; i<2; i++)
         {
           p=particle[i].head[s]->pt;
           while(p)
           {
             frontP[8*n+0]=p->x1;     
             frontP[8*n+1]=p->oldX1;  
             frontP[8*n+2]=p->p1;     
             frontP[8*n+3]=p->p2;     
             frontP[8*n+4]=p->p3;     
             frontP[8*n+5]=p->index;  
             frontP[8*n+6]=(float)s;  
             frontP[8*n+7]=(float)i;  
             p=p->next;
             n++;
           }
         }
       }
       MPI_Send(frontP,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
    }
    
    else if(myrank%2==1 && myrank!=nTasks-1) 
    {    
      MPI_Recv(frontP,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
         s=(int)(frontP[8*n+6]);
         i=(int)(frontP[8*n+7]);
         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[D->nxSub+i].head[s]->pt;
         particle[D->nxSub+i].head[s]->pt = New;      

         New->x1=frontP[8*n+0];
         New->oldX1=frontP[8*n+1]+D->nxSub;
         New->p1=frontP[8*n+2];
         New->p2=frontP[8*n+3];
         New->p3=frontP[8*n+4];
         New->index=frontP[8*n+5];
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(frontP);
       
    //Odd - evem
    if(myrank%2==1)
    {
       numP=0;
       for(s=0; s<D->nSpecies; s++)
       {
         for(i=0; i<2; i++)
         {
           p=particle[i].head[s]->pt;
           while(p)   {
             p=p->next;
             numP++;
           }
         }
       }      
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=nTasks-1) 
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    numberData=numP*8;
    frontP=(float *)malloc(numberData*sizeof(float ));             

    if(myrank%2==1)
    {
       n=0;
       for(s=0; s<D->nSpecies; s++)
       {
         for(i=0; i<2; i++)
         {
           p=particle[i].head[s]->pt;
           while(p)
           {
             frontP[8*n+0]=p->x1;     
             frontP[8*n+1]=p->oldX1;  
             frontP[8*n+2]=p->p1;     
             frontP[8*n+3]=p->p2;     
             frontP[8*n+4]=p->p3;     
             frontP[8*n+5]=p->index;  
             frontP[8*n+6]=(float)s;  
             frontP[8*n+7]=(float)i;  
             p=p->next;
             n++;
           }
         }
       }
       MPI_Send(frontP,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
    }
    
    else if(myrank%2==0 && myrank!=nTasks-1)  
    {    
      MPI_Recv(frontP,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
         s=(int)(frontP[8*n+6]);
         i=(int)(frontP[8*n+7]);
         New = (ptclList *)malloc(sizeof(ptclList)); 
         New->next = particle[D->nxSub+i].head[s]->pt;
         particle[D->nxSub+i].head[s]->pt = New;      

         New->x1=frontP[8*n+0];
         New->oldX1=frontP[8*n+1]+D->nxSub;
         New->p1=frontP[8*n+2];
         New->p2=frontP[8*n+3];
         New->p3=frontP[8*n+4];
         New->index=frontP[8*n+5];
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(frontP);

}

