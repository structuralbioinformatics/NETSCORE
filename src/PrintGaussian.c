#include "ppin.h"

void   PrintGaussian(OUT,data1,n1,m1,data2,n2,m2,bins1,bins2,nbins)
FILE     *OUT;
float    *data1,*data2,*bins1,*bins2;
int       n1,n2,m1,m2,nbins;
{

  int    i,j;
  float  max1,max2,min1,min2,min,max,d,step1,step2,*bin1,*bin2;
  float  *fvector();
  void   free_fvector();

  max1=max2=min1=min2=0.0;
  bin1=fvector(0,MAXSTPG+2);
  bin2=fvector(0,MAXSTPG+2);

  for (i=n1;i<m1;i++){
    if (i==0 || data1[i]>max1) max1=data1[i];
    if (i==0 || data1[i]<min1) min1=data1[i];
    }
  for (i=n2;i<m2;i++){
    if (i==0 || data2[i]>max2) max2=data2[i];
    if (i==0 || data2[i]<min2) min2=data2[i];
    }
  if (min1>min2)min=min2;else min=min1;  
  if (max1>max2)max=max1;else max=max2; 
  d=(max-min)/MAXSTPG;
  for (j=-1;j<MAXSTPG+1;j++){
      step1=d*j + min;
      step2=d*(j+1) + min;
      bin1[j+1]=0;
      for (i=n1;i<m1;i++)
       if (data1[i]>=step1 && data1[i]<step2 ) bin1[j+1]++;
      if (m1>n1) bin1[j+1]=bin1[j+1]/(m1-n1);
      }
  
  for (j=-1;j<MAXSTPG+1;j++){
      step1=d*j + min;
      step2=d*(j+1) + min;
      bin2[j+1]=0;
      for (i=n2;i<m2;i++)
       if (data2[i]>=step1 && data2[i]<step2 ) bin2[j+1]++;
      if (m2>n2) bin2[j+1]=bin2[j+1]/(m2-n2);
      }    
   fprintf(OUT,"%10s\t%10s\t%10s\n","Score","Distr[1]","Distr[2]");
   nbins=0;
   for (j=-1;j<MAXSTPG+1;j++){
      step1=d*j + min + d/2.0;
      if (bin1[j+1]>0 ||  bin2[j+1]>0 || j==-1 || j==MAXSTPG){
         fprintf(OUT,"%10.5e\t%10.5f\t%10.5f\n",step1,bin1[j+1],bin2[j+1]);
         nbins++;
         bins1[nbins]=bin1[j+1];
         bins2[nbins]=bin2[j+1];
         }
      }  
   fprintf(OUT,"----------------------------------\n\n");

   free_fvector(bin1,0,MAXSTPG+2);
   free_fvector(bin2,0,MAXSTPG+2);
   
}
