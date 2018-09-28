#include "ppin.h"

void   BinChi(data1,m1,data2,m2,statistic,prob)
float *data1,*data2,*statistic,*prob;
int m1,m2;
{

  int    i,j,n1,n2,nbins,kmax,kbins;
  float  max1,max2,min1,min2,min,max,d,step1,step2,*bin1,*bin2,freedom;
  float  *bins1,*bins2;  
  float  *fvector();
  int    SelectNbins2();
  void   free_fvector(),chstwo();


  kmax=100;

  kbins=SelectNbins2(data1,1,m1,data2,1,m2,kmax);


  n1=n2=1;
  max1=max2=min1=min2=0.0;
  bin1=fvector(0,kbins+2);
  bin2=fvector(0,kbins+2);
  bins1=fvector(0,kbins+2);
  bins2=fvector(0,kbins+2);

  for (i=n1;i<m1;i++){
    if (i==n1 || data1[i]>max1) max1=data1[i];
    if (i==n1 || data1[i]<min1) min1=data1[i];
    }
  for (i=n2;i<m2;i++){
    if (i==n2 || data2[i]>max2) max2=data2[i];
    if (i==n2 || data2[i]<min2) min2=data2[i];
    }
  if (min1>min2)min=min2;else min=min1;  
  if (max1>max2)max=max1;else max=max2;

 
  d=(max-min)/kbins;
  for (j=-1;j<kbins+1;j++){
      step1=d*j + min;
      step2=d*(j+1) + min;
      bin1[j+1]=0;
      for (i=n1;i<m1;i++)
       if (data1[i]>=step1 && data1[i]<step2 ) bin1[j+1]++;
      }
  
  for (j=-1;j<kbins+1;j++){
      step1=d*j + min;
      step2=d*(j+1) + min;
      bin2[j+1]=0;
      for (i=n2;i<m2;i++)
       if (data2[i]>=step1 && data2[i]<step2 ) bin2[j+1]++;
      }    
   nbins=0;
   for (j=-1;j<kbins+1;j++){
      step1=d*j + min + d/2.0;
      if (bin1[j+1]>0 ||  bin2[j+1]>0 || j==-1 || j==kbins){
         nbins++;
         bins1[nbins]=bin1[j+1];
         bins2[nbins]=bin2[j+1];
         }
      }  
   chstwo(bins1,bins2,nbins,0,&freedom,statistic,prob);

   free_fvector(bin1,0,kbins+2);
   free_fvector(bin2,0,kbins+2);
   free_fvector(bins1,0,kbins+2);
   free_fvector(bins2,0,kbins+2);
   


}
