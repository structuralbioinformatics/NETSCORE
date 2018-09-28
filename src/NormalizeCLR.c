#include "ppin.h"
#define EPS 1.0e-6

void   NormalizeCLR(c,m,n)
float **c;
int     m,n;
{
  gauss  normal(),g1,g2;
  int    i,j,k;
  float   *x,**y,z1,z2;
  float  **matrix(),*fvector();
  void    free_matrix(),free_fvector();

  x=fvector(0,n);
  y=matrix(0,m,0,n);
  for (i=0;i<m;i++){
  for (j=i+1;j<m;j++){
   for (k=0;k<n;k++) x[k]=c[i][k];
   g1=normal(x,n);
   for (k=0;k<n;k++) x[k]=c[j][k];
   g2=normal(x,n);
   if (g1.rmsd>EPS){ z1=(c[i][j]-g1.average)/g1.rmsd; }else{z1=0.0;}
   if (g2.rmsd>EPS){ z2=(c[i][j]-g2.average)/g2.rmsd; }else{z2=0.0;}
   if (z1*z2>EPS){ y[i][j]=y[j][i]=sqrt(z1*z1+z2*z2); }else{y[i][j]=y[j][i]=0.0;}
   if (z1*z2<0) y[i][j]=y[j][i]= (-1)*sqrt(z1*z1+z2*z2);
  }}
  for (i=0;i<m;i++)for (j=i+1;j<m;j++) c[i][j]=c[j][i]=y[i][j];
  free_matrix(y,0,m,0,n);
  free_fvector(x,0,n);

}
