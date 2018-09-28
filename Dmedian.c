#include "ppin.h"

float Dmedian(data1,m1,data2,m2)
float data1[],data2[];
int m1,m2;
{
 float x,max1,min1,med1,max2,min2,med2,max,min;
 int i;
 max1=max2=min1=min2=med1=med2=x=0;
 for (i=1;i<=m1;i++){
     if (i==1) {max1=data1[i];min1=data1[i];}
     if (i>1 && max1 < data1[i])max1=data1[i];
     if (i>1 && min1 > data1[i])min1=data1[i];
    }
 med1=(max1-min1)/2;
 for (i=1;i<=m2;i++){
     if (i==1) {max2=data2[i];min2=data1[i];}
     if (i>1 && max2 < data2[i])max2=data1[i];
     if (i>1 && min2 > data2[i])min2=data1[i];
    }
 med2=(max2-min2)/2;
 if (min1<min2) min=min1; else min=min2;
 if (max1>max2) max=max1; else max=max2;

 if (max>min) x=fabs(med2-med1)/(max-min);

 return x;
}
