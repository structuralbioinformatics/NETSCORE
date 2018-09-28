#include "ppin.h"

void   Bin2Cramer(data1,m1,data2,m2,chisq,df,prob,cramrv,ccc)
float *data1,*data2;
float *chisq,*df,*prob,*cramrv,*ccc;
int m1,m2;
{


  int    i,j,n,k,ni,nj,**nn,kmax,kbins1,kbins2,kbins;
  float stepi1,stepi2,stepj1,stepj2;
  float max1,max2,min1,min2,h1,h2;
  int    **imatrix();
  int    SelectNbins();
  void   free_imatrix(),cntab1();

  if (m1!=m2) nrerror("Inconsistent Statistical test, condition M1=M2 is required");

  kmax=MAXHPDF;
  //printf("Data1\n");
  kbins1= SelectNbins(data1,1,m1,kmax);
  //printf("Data2\n");
  kbins2= SelectNbins(data2,1,m2,kmax);
  kbins1=(int)sqrt(kbins1);
  kbins2=(int)sqrt(kbins2);
  if (kbins1<2)kbins1=2; 
  if (kbins2<2)kbins2=2;
 
  for (i=1;i<m1;i++){
    if (i==1 || data1[i]>max1) max1=data1[i];
    if (i==1 || data1[i]<min1) min1=data1[i];
    }
  for (i=1;i<m2;i++){
    if (i==1 || data2[i]>max2) max2=data2[i];
    if (i==1 || data2[i]<min2) min2=data2[i];
    }
  h1=(max1-min1)/kbins1;
  h2=(max2-min2)/kbins2;


  nn=imatrix(0,kbins1+3,0,kbins2+3);
  

  for (i=-1;i<kbins1+1;i++){
  for (j=-1;j<kbins2+1;j++){
    stepi1=h1*i + min1;
    stepi2=h1*(i+1) + min1;
    stepj1=h2*j + min2;
    stepj2=h2*(j+1) + min2;
    nn[i+2][j+2]=0;
    for (k=1;k<m1;k++)
    {
       if (   data1[k]>=stepi1 && data1[k]<stepi2 
           && data2[k]>=stepj1 && data2[k]<stepj2) nn[i+2][j+2]++;
    }

  }}

  ni=kbins1+2;
  nj=kbins2+2;


  cntab1(nn,ni,nj,chisq,df,prob,cramrv,ccc); 

  free_imatrix(nn,0,kbins+3,0,3);

}

