#include "ppin.h"

void   EvalueNode(rnd,scr,e,size,a,ma,funcs)
float  *scr,*a,*e;
int     *rnd,size,ma;
void    (*funcs)();
{
 int    i,j;
 float  x, y,s1,s2,eps,min,max,e_max,e_min, *dyda, *fvector(),qromo(),qtrap();
 float  (*funk)(),(*choose)(),midexp(),midpnt(),midinf(),ffgauss(),ffdgumbel();
 void   free_fvector();
 
 
 dyda=fvector(0,ma+1);
 printf("Assigned E-values on scores\n");
 if (rnd[6]==1||rnd[6]==3){
   max=min=scr[0];
   for (i=0;i<size;i++){
    x=scr[i];
    if (min>x)min=x;
    if (max<x)max=x;
    (*funcs)(x,a,&y,dyda,ma);
    e[i]=1-y;
    printf("Node[%5d] \tE-value[%f]=\t%e\n",i,x,e[i]);
    }
   x=min-1;
   (*funcs)(x,a,&y,dyda,ma);
   e_min=1-y;  
   x=max+1;
   (*funcs)(x,a,&y,dyda,ma);
   e_max=1-y;
   printf("Extreme E-values: Minimum[%f] (%e)  Maximum[%f] (%e)\n",(min-1),e_min,(max+1),e_max); 
 }
 if (rnd[6]==2||rnd[6]==4){
   if (rnd[6]==4) funk=&ffgauss;
   if (rnd[6]==2) funk=&ffdgumbel;
   //choose=&midexp;
   max=min=scr[0];
   choose=&midinf;
   for (i=0;i<size;i++){
     x=scr[i];
     if (min>x)min=x;
     if (max<x)max=x;
   }
   choose=&midinf;
   s2=qromo(funk,max,1.e+30,choose,a,ma);  
   for (i=0;i<size;i++){
    x=scr[i];
    choose=&midpnt;
    s1=qromo(funk,x,max,choose,a,ma);
    //printf("Test QROMO S1=%e\n",s1);
    //s1=qtrap(funk,x,max,a,ma);
    //printf("Test QTRAP S1=%e\n",s1);
    y=s1+s2;
    e[i]=y; 
    printf("Node[%5d] \tE-value[%f]=\t%e\n",i,x,e[i]);
   } 
   x=min-1;
   choose=&midpnt;
   s1=qromo(funk,x,max,choose,a,ma);
   //s1=qtrap(funk,x,max,a,ma);
   y=s1+s2;
   e_min=y;  
   x=max+1;
   choose=&midinf;
   e_max=qromo(funk,x,1.e+30,choose,a,ma);
   printf("Extreme E-values: Minimum[%f] (%e)  Maximum[%f] (%e)\n",(min-1),e_min,(max+1),e_max); 
 }
 free_fvector(dyda,0,ma+1);

}


