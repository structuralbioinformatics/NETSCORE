#include "ppin.h"

void WFloydWarshall(a,n)
float **a;
int n;
{

 int i, j,k;
 float infinity;

 if (MINSHORT>0){infinity=(float)n/MINSHORT;}
 else           {infinity=1.0e+10;}

 for( k = 0; k < n; k++){
  for( i = 0; i < n; i++){
    if (k!=i){
     //if (a[k][i]>=1./MINSHORT ) {a[k][i]=infinity;}
     if (a[k][i]==0 ) {a[k][i]=infinity;}
    }
  }
 }

 for( k = 0; k < n; k++){
  for( i = 0; i < n; i++){
   if  (a[k][i]<infinity && k!=i){
    for( j = 0; j < n; j++){
      if(a[i][j]>a[i][k]+a[k][j]){
         a[i][j]=a[i][k]+a[k][j];
      }
    }
   }
  }
  //printf(".%d.",k);
 }
 //printf("\n");
}
