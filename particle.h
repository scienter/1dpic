

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    float x1; 
    float oldX1;   
    float p1;    //momentum  
    float p2;
    float p3;
    float E1;    
    float Pr;    
    float Pl;    
    float Sr;    
    float Sl;    
    int index; 
    struct _ptclList *next;
} ptclList;

