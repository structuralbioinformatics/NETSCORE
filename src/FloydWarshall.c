#include <time.h>
#include        <string.h>
#include        <stdlib.h>
#include        <stdio.h>

void FloydWarshall(a,n)
int **a,n;
{
 int i, j,k;

 for( k = 0; k < n; k++){
  for( i = 0; i < n; i++){
    if (a[k][i]==0 && k!=i) {a[k][i]=n+1;}
  }
 }
 
 for( k = 0; k < n; k++){
  for( i = 0; i < n; i++){
   if  (a[k][i]<n+1 && k!=i ){
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
