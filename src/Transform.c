#include "ppin.h"
float Transform(i,maxmin,size,d_scr)
int       i,*size;
float     maxmin[];
gradient **d_scr;
{
 float e,max,min;
 int   j;

 max=maxmin[0];
 min=maxmin[1];
 e=0.0;
 for (j=0;j<size[i];j++){
  if (d_scr[i][j].j!=i){ e+=d_scr[i][j].size;}
 }
 if (e>max) {e=max;}
 if (e<min) {e=0;}
 return e;

}
