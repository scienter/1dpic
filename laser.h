typedef struct _LaserList  {
   int polarity;
   float lambda;
   float omega;
   float amplitude;   //unit is a0.
   float rU;
   float flat;
   float rD;
   int loadPoint; 

   struct _LaserList *next;
} LaserList;
