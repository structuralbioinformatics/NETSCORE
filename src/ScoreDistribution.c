#include "ppin.h"

void ScoreDistribution(rnd,interval,x,scr,size,kmax,pdf)
int    *rnd,kmax,size;
float  interval,*x,*scr,*pdf;
{
  int   i, k ;
  float max,min,a,b,sum;

  max=scr[0];
  for (i=0;i<size;i++){if (max<scr[i]){max=scr[i];}}
  min=scr[0];
  for (i=0;i<size;i++){if (min>scr[i]){min=scr[i];}}


  for (k=1;k<=kmax;k++) pdf[k]=0;
  
  for (k=1;k<=kmax;k++){
     a=x[k]-interval/2;
     b=x[k]+interval/2;
     if (rnd[6]==2 || rnd[6]==4){for (i=0;i<size;i++){ if (scr[i]>a && scr[i]<=b){ pdf[k]+=1;}}}
     if (rnd[6]==1 || rnd[6]==3){for (i=0;i<size;i++){ if (scr[i]<x[k]){ pdf[k]+=1;}}}
  }

  if (rnd[6]==2 || rnd[6]==4){
   b=x[kmax]+interval/2;
   for (i=0;i<size;i++){ if (scr[i]<=max && scr[i]> b){ pdf[kmax] +=1;}}
   a=x[0]-interval/2;
   for (i=0;i<size;i++){ if (scr[i]>=min && scr[i]< a){ pdf[1] +=1;}}
  }
  if (rnd[6]==1 || rnd[6]==3){
   b=x[kmax];
   for (i=0;i<size;i++){ if (scr[i]<=max && scr[i]> b){ pdf[kmax] +=1;}}
  }

  if (rnd[6]==1||rnd[6]==3){for (i=1;i<=kmax;i++){pdf[i]=pdf[i]/size;}}
  if (rnd[6]==2||rnd[6]==4){
    sum=0;
    for (k=1;k<=kmax;k++) sum+=pdf[k];
    if (sum>0) {for (k=1;k<=kmax;k++) {pdf[k]=pdf[k]/(sum*interval);}}
   }
  //for (k=1;k<=kmax;k++){printf("X[%d]=%f PDF[%d]=%f\n",k,x[k],k,pdf[k]);}
}
