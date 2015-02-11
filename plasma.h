#define Electron 	1
#define HPlus0 	100
#define HPlus1 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0	 	600
#define CPlus1	 	601
#define CPlus2	 	602
#define CPlus3	 	603
#define CPlus4	 	604
#define CPlus5	 	605
#define CPlus6	 	606
#define userDefined   	9999999

typedef struct _LoadList  {
   int species;
   float superP;
   float density;
   float numberInCell;
   float criticalDensity;
   int index;  
   float num;      //exceeded number of particle which is less than 1
   int cnt;
   int lnodes;     //longitudinal point number
   float *ln;      //longitudinal density (1 is P->density)
   float *lpoint;    //longitudinal point

   int pointPosition;   
   float p1;
   float p2;
   float p3;

   float beta;
   float gamma;
   float mass;
   int charge;
   
   float temperature;
   int withNextSpcs;
   int withPrevSpcs;
   
   float ioniEnergy;
  
   struct _LoadList *next;
} LoadList;
