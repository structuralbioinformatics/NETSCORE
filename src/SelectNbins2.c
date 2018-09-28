#include "ppin.h"

int SelectNbins2(data1,n1,m1,data2,n2,m2,max)
float *data1,*data2;
int   n1,n2,m1,m2,max;
{
 int    i,j,k,q,kmax,k_shi,n,nn,kmin;
 float  width,dmax,dmin,sigma,q1,q2,q3,ave,v,h,*bin;
 float  shi, shimazaki,sturges,log2n,square,step2,step1,scott,freedman_diaconis;
 float  *fvector();
 void   free_fvector();

 if (max<=0)max=MAXSTPG;

 bin=fvector(0,max+3); 

  dmin=dmax=0.0; 
  for (i=n1;i<m1;i++){
    if (i==n1 || data1[i]>dmax) dmax=data1[i];
    if (i==n1 || data1[i]<dmin) dmin=data1[i];
    }
  for (i=n2;i<m2;i++){
    if (data2[i]>dmax) dmax=data2[i];
    if (data2[i]<dmin) dmin=data2[i];
    }
 
 nn=(m2-n2) + (m1-n1);

 if ( nn>0) log2n= log(nn)/log(2);
 else log2n=1;
 sturges=log2n+1;
 k=(int)sturges;
 if (k>max)k=max;
 h=fabs(dmax-dmin)/k;
  for (j=0;j<k;j++){
    step2=dmin+h*(j+1);
    step1=dmin+h*j;
    bin[j]=0;
    for (i=n1;i<m1;i++) if (data1[i]<step2 && data1[i]>=step1 )bin[j]++;
    for (i=n2;i<m2;i++) if (data2[i]<step2 && data2[i]>=step1 )bin[j]++;
   }
  ave=v=0.0;
  for (j=0;j<k;j++){ave+=bin[j];v+=bin[j]*bin[j];}
  ave=ave/(k);
  v=v/(k);
  v=v-ave*ave;
  shi=(2*ave-v)/(nn*nn*h*h);
  shimazaki=shi;
  k_shi=k;
  kmin=k;
  kmax=k;
  if (kmax>=max)kmax=max;
 //printf("Shimazaki=%f K=%d Sturges_h=%f Sturges_k=%d\n",shimazaki,k_shi,h,k);

 if ( nn>0)square=sqrt(nn);
 else square=2;
 k=(int)square;
 if (k>max)k=max;
 h=fabs(dmax-dmin)/k;
  for (j=0;j<k;j++){
    step2=dmin+h*(j+1);
    step1=dmin+h*j;
    bin[j]=0;
    for (i=n1;i<m1;i++) if (data1[i]<step2 && data1[i]>=step1 )bin[j]++;
    for (i=n2;i<m2;i++) if (data2[i]<step2 && data2[i]>=step1 )bin[j]++;
   }
  ave=v=0.0;
  for (j=0;j<k;j++){ave+=bin[j];v+=bin[j]*bin[j];}
  ave=ave/(k);
  v=v/(k);
  v=v-ave*ave;
  shi=(2*ave-v)/(nn*nn*h*h);
  if (shi<=shimazaki) {shimazaki=shi;k_shi=k;}

  if (k<kmin)kmin=k;
  if (k>kmax)kmax=k;
  if (kmax>=max)kmax=max;
 //printf("Shimazaki=%f K=%d Square_h=%f Square_k=%d\n",shimazaki,k_shi,h,k);

 sigma=0.0;
 ave=0.0;
 for (i=n1;i<m1;i++){sigma+=data1[i]*data1[i];ave=data1[i];}
 for (i=n2;i<m2;i++){sigma+=data2[i]*data2[i];ave=data2[i];}
 if (nn>0) ave=ave/nn;
 if (nn>0) sigma=sigma/nn;
 if (nn>0) sigma=sigma-ave*ave; 
 if (nn>0) h=3.5*sigma/exp(log(nn)/3); 
 else      h=fabs(dmax-dmin)/2;
 scott=fabs(dmax-dmin)/h;
 k=(int)scott;
 if (k>max)k=max;
 h=fabs(dmax-dmin)/k;
  for (j=0;j<k;j++){
    step2=dmin+h*(j+1);
    step1=dmin+h*j;
    bin[j]=0;
    for (i=n1;i<m1;i++) if (data1[i]<step2 && data1[i]>=step1 )bin[j]++;
    for (i=n2;i<m2;i++) if (data2[i]<step2 && data2[i]>=step1 )bin[j]++;
   }
  ave=v=0.0;
  for (j=0;j<k;j++){ave+=bin[j];v+=bin[j]*bin[j];}
  ave=ave/(k);
  v=v/(k);
  v=v-ave*ave;
  shi=(2*ave-v)/(nn*nn*h*h);
  if (shi<=shimazaki) {shimazaki=shi;k_shi=k;}

  if (k<kmin)kmin=k;
  if (k>kmax)kmax=k;
  if (kmax>=max)kmax=max;
 //printf("Shimazaki=%f K=%d Scott_h=%f Scott_k=%d\n",shimazaki,k_shi,h,k);


 h=fabs(dmax-dmin)/max;
 q1=q2=q3=0.0;
 for (j=0;j<max;j++){
   step2=dmin+h*(j+1);
   q=0;
   for (i=n1;i<m1;i++) if (data1[i]<step2) q++;
   for (i=n2;i<m2;i++) if (data2[i]<step2) q++;
   if ((float)q/nn > 0.25 && q1==0) q1=step2;  
   if ((float)q/nn > 0.50 && q2==0) q2=step2;  
   if ((float)q/nn > 0.75 && q3==0) q3=step2;  
 }
   
 if (nn>0)  h= 2*(q3-q1)/sqrt(nn);
 else       h= fabs(dmax-dmin)/2;
 freedman_diaconis=fabs(dmax-dmin)/h;
 k=(int)freedman_diaconis;
 if (k>max)k=max;
 h=fabs(dmax-dmin)/k;
  for (j=0;j<k;j++){
    step2=dmin+h*(j+1);
    step1=dmin+h*j;
    bin[j]=0;
    for (i=n1;i<m1;i++) if (data1[i]<step2 && data1[i]>=step1 )bin[j]++;
    for (i=n2;i<m2;i++) if (data2[i]<step2 && data2[i]>=step1 )bin[j]++;
   }
  ave=v=0.0;
  for (j=0;j<k;j++){ave+=bin[j];v+=bin[j]*bin[j];}
  ave=ave/(k);
  v=v/(k);
  v=v-ave*ave;
  shi=(2*ave-v)/(nn*nn*h*h);
  if (shi<=shimazaki) {shimazaki=shi;k_shi=k;}

  if (k<kmin)kmin=k;
  if (k>kmax)kmax=k;
  if (kmax>=max)kmax=max;
 //printf("Shimazaki=%f K=%d freedman_h=%f freedman_k=%d\n",shimazaki,k_shi,h,k);

 if (kmin==1)kmin=2;
 for (k=kmin;k<kmax+2;k++ )
 {
  for (j=0;j<max;j++)bin[j]=0.0;
  h=fabs(dmax-dmin)/k;
  for (j=0;j<k;j++){
    step2=dmin+h*(j+1);
    step1=dmin+h*j;
    bin[j+1]=0;
    for (i=n1;i<m1;i++) if (data1[i]<step2 && data1[i]>=step1 )bin[j]++;
    for (i=n2;i<m2;i++) if (data2[i]<step2 && data2[i]>=step1 )bin[j]++;
   }
  ave=v=0.0;
  for (j=0;j<k;j++){ave+=bin[j];v+=bin[j]*bin[j];}
  ave=ave/k;
  v=v/k;
  v=v-ave*ave;
  shi=(2*ave-v)/(nn*nn*h*h);
  if (shi<=shimazaki) {shimazaki=shi;k_shi=k;}
 }

 //printf("Shimazaki(MAX to 2)=%f K=%d \n",shimazaki,k_shi);



 h=fabs(dmax-dmin)/k_shi; 

 free_fvector(bin,0,max+3);
 return  k_shi;
}

