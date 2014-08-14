#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void filterField(Domain *D)
{
    int i,a;  
    float alpha[2];
    float nowE1,nowPr,nowPl,nowSr,nowSl;
    float beforeE1,beforePr,beforePl,beforeSr,beforeSl;
    FieldElement *field;
    field=D->field;
    int myrank,nTasks;  
 
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    alpha[0]=0.5;
    alpha[1]=1.5;

    if(myrank==0)
    {    
       i=0;
       beforeE1=field[0].E1-(field[1].E1-field[0].E1); 
       beforePr=field[0].Pr-(field[1].Pr-field[0].Pr); 
       beforePl=field[0].Pl-(field[1].Pl-field[0].Pl); 
       beforeSr=field[0].Sr-(field[1].Sr-field[0].Sr); 
       beforeSl=field[0].Sl-(field[1].Sl-field[0].Sl); 

       nowE1=field[i].E1;
       nowPr=field[i].Pr;
       nowPl=field[i].Pl;
       nowSr=field[i].Sr;
       nowSl=field[i].Pl;

       for(a=0; a<1; a++)
       {
//          field[i].E1=alpha[a]*nowE1+(1.0-alpha[a])*(beforeE1+field[i+1].E1)*0.5;       
          field[i].Pr=alpha[a]*nowPr+(1.0-alpha[a])*(beforePr+field[i+1].Pr)*0.5;       
          field[i].Pl=alpha[a]*nowPl+(1.0-alpha[a])*(beforePl+field[i+1].Pl)*0.5;       
          field[i].Sr=alpha[a]*nowSr+(1.0-alpha[a])*(beforeSr+field[i+1].Sr)*0.5;     
          field[i].Sl=alpha[a]*nowSl+(1.0-alpha[a])*(beforeSl+field[i+1].Sl)*0.5;      
       }
    }

    else
    {
       nowE1=field[0].E1;
       nowPr=field[0].Pr;
       nowPl=field[0].Pl;
       nowSr=field[0].Sr;
       nowSl=field[0].Pl;
    }
          
    for(a=0; a<2; a++)
    {
       for(i=1; i<=D->nxSub; i++)
       {              
          beforeE1=nowE1;
          beforePr=nowPr;
          beforePl=nowPl;
          beforeSr=nowSr;
          beforeSl=nowSl;
          nowE1=field[i].E1;
          nowPr=field[i].Pr;
          nowPl=field[i].Pl;
          nowSr=field[i].Sr;
          nowSl=field[i].Sl;
//          field[i].E1=alpha[a]*nowE1+(1.0-alpha[a])*(beforeE1+field[i+1].E1)*0.5;       
          field[i].Pr=alpha[a]*nowPr+(1.0-alpha[a])*(beforePr+field[i+1].Pr)*0.5;       
          field[i].Pl=alpha[a]*nowPl+(1.0-alpha[a])*(beforePl+field[i+1].Pl)*0.5;       
          field[i].Sr=alpha[a]*nowSr+(1.0-alpha[a])*(beforeSr+field[i+1].Sr)*0.5;     
          field[i].Sl=alpha[a]*nowSl+(1.0-alpha[a])*(beforeSl+field[i+1].Sl)*0.5;       
       }
   }
    
    

}


