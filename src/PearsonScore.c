#include "ppin.h"

float PearsonScore(a,b,dim)
float *a,*b;
int    dim;
{
 float c, ma, mb, da, db, dab ;
 int   i,j;
 c=ma=mb=da=db=dab=0.0;
 if (dim<=0) return c;
 for (i=0;i<dim;i++){
     ma+=a[i];
     mb+=b[i];
     }
 ma=ma/dim;
 mb=mb/dim;
 for (i=0;i<dim;i++){
    da+=(a[i]-ma)*(a[i]-ma);
    db+=(b[i]-mb)*(b[i]-mb);
    dab+=(a[i]-ma)*(b[i]-mb);
    }
 if(da*db>0)c=dab/sqrt(da*db);
 
 return c;
}
