#include "ppin.h"
gauss normal(x,n)
float  *x;
int     n;
{
 gauss g;
 float a,aa,s;
 int  i;
 
 aa=a=0;
 for (i=0;i<n;i++){
  a+=x[i];
  aa+=x[i]*x[i];
 }

 if (n>0){a=a/n;aa=aa/n;}else{a=0.0;aa=0.0;}
 g.average=a;
 s=aa-a*a;
 if (s>0) {g.rmsd=sqrt(s);}else{g.rmsd=0.0;}

 return g;
 
}
