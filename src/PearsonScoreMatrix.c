#include "ppin.h"

float PearsonScoreMatrix(a,dim,n,m)
float **a;
int    dim,n,m;
{
 float c, ma, mb, da, db, dab ;
 int   i,j;
 c=ma=mb=da=db=dab=0.0;
 if (dim<=0) return c;
 for (i=0;i<dim;i++){
     ma+=a[n][i];
     mb+=a[m][i];
     }
 ma=ma/dim;
 mb=mb/dim;
 for (i=0;i<dim;i++){
    da+=(a[n][i]-ma)*(a[n][i]-ma);
    db+=(a[m][i]-mb)*(a[m][i]-mb);
    dab+=(a[n][i]-ma)*(a[m][i]-mb);
    }
 if(da*db>0)c=dab/sqrt(da*db);
 
 return c;
}
