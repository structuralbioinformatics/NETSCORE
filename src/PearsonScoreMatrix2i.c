#include "ppin.h"

float PearsonScoreMatrix2i(a,b,f,m1,m2)
float **a;
int   **b,**f;
int    m1,m2;
{
 float c, ma, mb, da, db, dab ;
 int   i,j,dim;
 c=ma=mb=da=db=dab=0.0;
 if (m1<=0 && m2<=0) return c;
 dim=0;
 for (i=0;i<m1;i++){
 for (j=0;j<m2;j++){
 if (f[i][j]==0){
     ma+=a[i][j];
     mb+=b[i][j];
     dim++;
     }}}
 if (dim==0)return c;
 ma=ma/dim;
 mb=mb/dim;
 for (i=0;i<m1;i++){
 for (j=0;j<m2;j++){
 if (f[i][j]==0){
    da+=(a[i][j]-ma)*(a[i][j]-ma);
    db+=(b[i][j]-mb)*(b[i][j]-mb);
    dab+=(a[i][j]-ma)*(b[i][j]-mb);
    }}}
 if(da*db>0)c=dab/sqrt(da*db);
 
 return c;
}
