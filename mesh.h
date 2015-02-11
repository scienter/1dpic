#include "particle.h"

typedef struct _Domain 
{
   int currentType;
   int interpolationType;

   int maxStep;
   int saveStep;
   int saveStart;
   int fieldSave;
   int particleSave;
   int rhoSave;
   int currentSave;
   int dumpSave;
   int dumpStart;

   int nx;             //Total domain
   int nxSub;          //Each core has sub domain        
   int minXSub;        //Each core has start mesh point in total domain
   int maxXSub;
//   int numberInCell;
   int moving;         //Moving domain option. 1:on

   float lambda;				// basic parameter for dx,dt
   float omega;
   float divisionLambda;   // dx=1/divisionLambda      
   float dt;
   float dx;

   //MPI parameter
   int L;         //x cores
   int beforeCore;
   int nextCore;
   
   struct _FieldElement *field;    
   struct _Particle *particle;
   struct _Boost *boost;    

   //Plasma load
   struct _LoadList *loadList;
   int nSpecies;
   int crystal;

   //Laser load
   struct _LaserList *laserList;
   int nLaser;

   //Boost
   int boostOn;
   float gamma;
   float beta;
   int minT;	//boost frame's step
   int maxT;	//boost frame's step
   int boostSaveStep;	//lab frame's step
   
   //Probe
   int probeNum;
   int *probeX;
   int *probeIndex;
   int isThereProbe;
   struct _Probe **probe;    

}  Domain; 

typedef struct _Boost
{
   float X1;
   float E1;
   float Pr;
   float Pl;
   float Sr;
   float Sl;   
}  Boost;

typedef struct _FieldElement 
{
   float E1;
   float Pr;
   float Pl;
   float Sr;
   float Sl;
   double J1;
   double J2;
   double J3;

}  FieldElement;

typedef struct _Particle
{
   float rho;
   // Particle List Header
   ptclHead **head;            

}  Particle;

typedef struct _Probe
{
   float E1;
   float Pr;
   float Pl;
   float B1;
   float Sr;
   float Sl;

   float x;
   float p1;
   float p2;
   float p3;

   float J1;
   float J2;
   float J3;
}  Probe;

typedef struct _External 
{
   float E1;
   float B1;
   float Pr;
   float Pl;
   float Sr;
   float Sl;

}  External;


typedef struct _Ionization
{
   int l;
   int m;
   int numberOfMaterials;
   int **material;

}  Ionization;
