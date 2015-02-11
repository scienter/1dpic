#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "laser.h"
#include "constants.h"
#include "plasma.h"
#include "mpi.h"

int main(int argc, char *argv[])
{
    int i,j,k,n,iteration=0,filter,boost,filterStep,labSaveStep;
    float factor;
    double t;
    char name[100];
    FILE *out;
    Domain D;
    Probe P;
    LaserList *L;  
    External Ext;
    Ionization I;
    FieldElement *field;
    int myrank, nTasks;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    void parameterSetting();
    void boundaray();    
    void saveField();  
    void saveRaman();  
    void saveParticle();
    void saveDensity();
    void saveProbe();
//    void boostShot();
    void solveField();
//    void filterField();
    void loadLaser();
//    void boostLoadLaser();
    void loadPlasma();    
    void interpolation_1st();
    void interpolation_2nd();
    void particlePush();       
    void updateCurrent_1st(); 
    void updateCurrent_2nd(); 
    void updateCurrent_3rd(); 
    void removeEdge();
    void movingDomain();
    void loadMovingPlasma();
    void probe();
    void MPI_TransferF_Xminus();
    void MPI_TransferF_Xplus();
//    void MPI_TransferF_XplusFilter();
    void MPI_TransferP_Xplus();
    void MPI_TransferP_Xminus();
    void MPI_TransferJ_Xminus();
    void MPI_TransferJ_Xplus();
//    void MPI_TransferP_Moving();
//    void MPI_TransferJ_Moving();
    void clean1D();
 
    if(FindParameters("Domain",1,"filter",argv[1],name)) filter=atoi(name);
    else  filter=0;
    if(FindParameters("Domain",1,"filterStep",argv[1],name)) filterStep=atoi(name);
    else  filterStep=10;

    //parameter setting
    parameterSetting(&D,&Ext,&I,argv[1]);


    //create mesh
    boundary(&D,&Ext);

    //load plasma or load dump file
    if(argc >= 3)
    {
      iteration = atoi(argv[2]);
      restoreData1D(&D,iteration);
      t=D.dt*iteration;
    }
    else
    {
      if(D.crystal==1) {
         loadPlasma_crystal(&D);
      }
      else 
         loadPlasma_random(&D);
      t=0;
    }

    //rooping time 
   labSaveStep=D.boostSaveStep;
   factor=D.gamma*(1+D.beta);

//    if(D.boostOn==1)      boostLoadLaser(&D);  
    while(iteration<=D.maxStep)
    {   
       //probe data
//       probe(&D,iteration);
//       if(D.probeNum>0)
//         findProbeParticle(&D);

       
       if(filter==1 && iteration%filterStep==0)
       {
          MPI_TransferF_Xminus(&D);
          MPI_TransferF_XplusFilter(&D);
          filterField(&D);       
       }
       
       if(D.boostOn==1)
       {
          if(iteration>=D.minT && iteration<=D.maxT)
             boostShot(&D,iteration,labSaveStep);    
          if(iteration>=D.maxT)   
          {

             boostSaveField(&D,labSaveStep);
exit(0);
          }
       }

       //save File      
       if(iteration%D.saveStep==0 && iteration>=D.saveStart)   
       {
          if(D.fieldSave==1)  {
            saveField(&D,iteration);
            saveRaman(&D,iteration);
            if(myrank==0)
              printf("field%d is made. \n",iteration);
          }
          if(D.particleSave==1)  {
            saveParticle(&D,iteration);
            if(myrank==0)
              printf("particle%d is made. \n",iteration);
          }
          if(D.rhoSave==1)  {
            saveRho(&D,iteration);
            if(myrank==0)
              printf("rho%d is made. \n",iteration);
          }
          if(D.probeNum>0)  {
            saveProbe(&D,iteration);
            if(myrank==0)
              printf("probe%d is made. \n",iteration);
          }           
          if(D.currentSave==1)  {
            saveCurrent(&D,iteration);
            if(myrank==0)
              printf("current%d is made. \n",iteration);
          }           
          if(D.dumpSave==1 && iteration>=D.dumpStart) { 
            saveDump1D(D,iteration);
            if(myrank==0)
              printf("dump%d is made.\n",iteration);              
          }
       }

       //load Laser
       if(D.boostOn==0)   {
         L=D.laserList;
         while(L->next)  {
           loadLaser(&D,L,&Ext,t);
           L=L->next;
         }
       }  
 
       solveField(&D);
       MPI_TransferF_Xminus(&D);
       MPI_TransferF_Xplus(&D);

       if(D.interpolationType==1)
         interpolation_1st(&D,&Ext);
       else if(D.interpolationType==2)
         interpolation_2nd(&D,&Ext);

       //field Ionization
       fieldIonization(&D,&I);

       //particle push
       particlePush(&D);

       if(D.currentType==1)
          updateCurrent_1st(&D);
       else if(D.currentType==2)
          updateCurrent_2nd(&D);
       else if(D.currentType==3)
          updateCurrent_3rd(&D);
       MPI_TransferJ_Xplus(&D);
       MPI_TransferJ_Xminus(&D);

       if(iteration>=D.nx && D.moving==1 && D.boostOn==0) {
          movingDomain(&D);
          if(myrank==nTasks-1)
          {
             if(D.crystal==0)
                loadMovingPlasma_random(&D);
             else if(D.crystal==1)
                loadMovingPlasma_crystal(&D);
          }          
          MPI_TransferF_Moving(&D);
       }

       else if(D.boostOn==1) {
//          MPI_TransferF_Xminus(&D);
 //         MPI_TransferP_Moving(&D);
//          MPI_TransferJ_Moving(&D);
//          if(myrank==nTasks-1)   
//             loadMovingPlasma(&D,&Ext);          
//          movingDomain(&D);
       }

       rearrangeParticles(&D);
       MPI_TransferP_Xminus(&D);
       MPI_TransferP_Xplus(&D);
       removeEdge(&D);

       //time update
       t+=D.dt;  
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;

    }     //end of time roop                  

    clean1D(&D);

    MPI_Finalize();

    return 0;
}
