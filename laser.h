typedef struct _LaserList  {
   int addition;	//ON : laser field is added on previous field value.
   int polarity;
   float lambda;
   float omega;
   float amplitude;   //unit is a0.
   float rU;
   float flat;
   float rD;
   int loadPoint; 
   float phs;  // 20150219 mshur, 

   struct _LaserList *next;
} LaserList;
