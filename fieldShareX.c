#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_TransferF_XplusFilter(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    behindF=(float *)malloc(numberData*sizeof(float )); 
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nxSub].E1;
    behindF[1]=D->field[D->nxSub].Pr;
    behindF[2]=D->field[D->nxSub].Pl;
    behindF[3]=D->field[D->nxSub].Sr;
    behindF[4]=D->field[D->nxSub].Sl;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nxSub].E1;
    behindF[1]=D->field[D->nxSub].Pr;
    behindF[2]=D->field[D->nxSub].Pl;
    behindF[3]=D->field[D->nxSub].Sr;
    behindF[4]=D->field[D->nxSub].Sl;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}


void MPI_TransferF_Xminus(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    float *frontF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=7;
    frontF=(float *)malloc(numberData*sizeof(float ));             
    //Transferring even ~ odd cores             
    frontF[0]=D->field[2].E1;
    frontF[1]=D->field[2].Pr;
    frontF[2]=D->field[2].Pl;
    frontF[3]=D->field[2].Sr;
    frontF[4]=D->field[2].Sl;
    frontF[5]=D->field[1].Pl;
    frontF[6]=D->field[1].Sl;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+2].E1=frontF[0];
       D->field[D->nxSub+2].Pr=frontF[1];
       D->field[D->nxSub+2].Pl=frontF[2];
       D->field[D->nxSub+2].Sr=frontF[3];
       D->field[D->nxSub+2].Sl=frontF[4];
       D->field[D->nxSub+1].Pl=frontF[5];
       D->field[D->nxSub+1].Sl=frontF[6];
    }
    else if(myrank%2==1)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    frontF[0]=D->field[2].E1;
    frontF[1]=D->field[2].Pr;
    frontF[2]=D->field[2].Pl;
    frontF[3]=D->field[2].Sr;
    frontF[4]=D->field[2].Sl;
    frontF[5]=D->field[1].Pl;
    frontF[6]=D->field[1].Sl;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+2].E1=frontF[0];
       D->field[D->nxSub+2].Pr=frontF[1];
       D->field[D->nxSub+2].Pl=frontF[2];
       D->field[D->nxSub+2].Sr=frontF[3];
       D->field[D->nxSub+2].Sl=frontF[4];
       D->field[D->nxSub+1].Pl=frontF[5];
       D->field[D->nxSub+1].Sl=frontF[6];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontF);
}

void MPI_TransferF_Xplus(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    float *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    behindF=(float *)malloc(numberData*sizeof(float )); 
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nxSub+1].E1;
    behindF[1]=D->field[D->nxSub+1].Pr;
    behindF[2]=D->field[D->nxSub+1].Pl;
    behindF[3]=D->field[D->nxSub+1].Sr;
    behindF[4]=D->field[D->nxSub+1].Sl;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       D->field[1].E1=behindF[0];
       D->field[1].Pr=behindF[1];
       D->field[1].Pl=behindF[2];
       D->field[1].Sr=behindF[3];
       D->field[1].Sl=behindF[4];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nxSub+1].E1;
    behindF[1]=D->field[D->nxSub+1].Pr;
    behindF[2]=D->field[D->nxSub+1].Pl;
    behindF[3]=D->field[D->nxSub+1].Sr;
    behindF[4]=D->field[D->nxSub+1].Sl;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       D->field[1].E1=behindF[0];
       D->field[1].Pr=behindF[1];
       D->field[1].Pl=behindF[2];
       D->field[1].Sr=behindF[3];
       D->field[1].Sl=behindF[4];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}

void MPI_TransferJ_Xminus(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    double *frontJ;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3*2;
    frontJ=(double *)malloc(numberData*sizeof(double ));             
    //Transferring even ~ odd cores          
    for(n=0; n<2; n++)
    {   
      frontJ[3*n+0]=D->field[n].J1;
      frontJ[3*n+1]=D->field[n].J2;
      frontJ[3*n+2]=D->field[n].J3;
    }
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       for(n=0; n<2; n++)
       {
         D->field[D->nxSub+n].J1+=frontJ[3*n+0];
         D->field[D->nxSub+n].J2+=frontJ[3*n+1];
         D->field[D->nxSub+n].J3+=frontJ[3*n+2];
       }
    }
    else if(myrank%2==1)
       MPI_Send(frontJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    for(n=0; n<2; n++)
    {   
      frontJ[3*n+0]=D->field[n].J1;
      frontJ[3*n+1]=D->field[n].J2;
      frontJ[3*n+2]=D->field[n].J3;
    }
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       for(n=0; n<2; n++)
       {
         D->field[D->nxSub+n].J1+=frontJ[3*n+0];
         D->field[D->nxSub+n].J2+=frontJ[3*n+1];
         D->field[D->nxSub+n].J3+=frontJ[3*n+2];
       }
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontJ);
}

void MPI_TransferJ_Xplus(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    double *behindJ;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3*3;
    behindJ=(double *)malloc(numberData*sizeof(double )); 
    
    //Transferring even ~ odd cores             
    for(n=0; n<3; n++)
    {   
      behindJ[3*n+0]=D->field[D->nxSub+2+n].J1;
      behindJ[3*n+1]=D->field[D->nxSub+2+n].J2;
      behindJ[3*n+2]=D->field[D->nxSub+2+n].J3;
    }
        
    if(myrank%2==1)
    {
       MPI_Recv(behindJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       for(n=0; n<3; n++)
       {
         D->field[2+n].J1+=behindJ[3*n+0];
         D->field[2+n].J2+=behindJ[3*n+1];
         D->field[2+n].J3+=behindJ[3*n+2];
       }
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    for(n=0; n<3; n++)
    {   
      behindJ[3*n+0]=D->field[D->nxSub+2+n].J1;
      behindJ[3*n+1]=D->field[D->nxSub+2+n].J2;
      behindJ[3*n+2]=D->field[D->nxSub+2+n].J3;
    }
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       for(n=0; n<3; n++)
       {
         D->field[2+n].J1+=behindJ[3*n+0];
         D->field[2+n].J2+=behindJ[3*n+1];
         D->field[2+n].J3+=behindJ[3*n+2];
       }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindJ);
}

void MPI_TransferJ_Moving(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    float *frontJ;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3;
    frontJ=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores             
    frontJ[0]=D->field[1].J1;
    frontJ[1]=D->field[1].J2;
    frontJ[2]=D->field[1].J3;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].J1=frontJ[0];
       D->field[D->nxSub+1].J2=frontJ[1];
       D->field[D->nxSub+1].J3=frontJ[2];
    }
    else if(myrank%2==1)
       MPI_Send(frontJ,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    frontJ[0]=D->field[1].J1;
    frontJ[1]=D->field[1].J2;
    frontJ[2]=D->field[1].J3;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].J1=frontJ[0];
       D->field[D->nxSub+1].J2=frontJ[1];
       D->field[D->nxSub+1].J3=frontJ[2];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontJ,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontJ);
}

void MPI_TransferF_Moving(Domain *D)
{
    int i,n,rank,numberData;
    int myrank, nTasks; 
    float *frontF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    frontF=(float *)malloc(numberData*sizeof(float ));             
    //Transferring even ~ odd cores             
    frontF[0]=D->field[2].E1;
    frontF[1]=D->field[2].Pr;
    frontF[2]=D->field[2].Pl;
    frontF[3]=D->field[2].Sr;
    frontF[4]=D->field[2].Sl;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+2].E1=frontF[0];
       D->field[D->nxSub+2].Pr=frontF[1];
       D->field[D->nxSub+2].Pl=frontF[2];
       D->field[D->nxSub+2].Sr=frontF[3];
       D->field[D->nxSub+2].Sl=frontF[4];
    }
    else if(myrank%2==1)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    frontF[0]=D->field[2].E1;
    frontF[1]=D->field[2].Pr;
    frontF[2]=D->field[2].Pl;
    frontF[3]=D->field[2].Sr;
    frontF[4]=D->field[2].Sl;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+2].E1=frontF[0];
       D->field[D->nxSub+2].Pr=frontF[1];
       D->field[D->nxSub+2].Pl=frontF[2];
       D->field[D->nxSub+2].Sr=frontF[3];
       D->field[D->nxSub+2].Sl=frontF[4];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontF);
}

