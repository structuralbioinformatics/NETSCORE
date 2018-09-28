#include "ppin.h"
float  CmpMatrix(old,new,n,m)
gradient **old,**new;
int     *n,m;
{
 float c,d;
 int   i,j;

 c=0.0;
 for (i=0;i<m;i++){
 for (j=0;j<n[i];j++){
  if (fabs(old[i][j].size)>1.0e-10){
   d=fabs(new[i][j].size-old[i][j].size);
  }else{
   d=fabs(new[i][j].size);
  }
  if (c<d){c=d;}
 }}
 
 return c;

}
