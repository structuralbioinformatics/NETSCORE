#include "ppin.h"
void   NormalizeASC(c,m)
float **c;
int     m;
{
  gauss  normal(),g1,g2;
  int    i,j,k,n;
  float   *x,**y,z;
  float  **matrix(),*fvector();
  void    free_matrix(),free_fvector();

  x=fvector(0,m);
  y=matrix(0,m,0,m);
  z=0.0;
  n=0;
  for (i=0;i<m;i++){for (j=0;j<m;j++) {z+=c[i][j];n++;}}
  if (n>0)z=z/n;
  if (z>0){
   for (i=0;i<m;i++){
   for (j=i+1;j<m;j++){
    for (k=0;k<m;k++) x[k]=c[i][k];
    g1=normal(x,m);
    for (k=0;k<m;k++) x[k]=c[j][k];
    g2=normal(x,m);
    y[i][j]= c[i][j] - g1.average - g2.average + z;
   }}
   for (i=0;i<m;i++)for (j=i+1;j<m;j++) c[i][j]=c[j][i]=y[i][j];
  }
  free_matrix(y,0,m,0,m);
  free_fvector(x,0,m);

}
