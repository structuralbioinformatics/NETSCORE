#include "ppin.h"

void   BinMI(data1,m1,data2,m2,uxy)
float *data1,*data2;
float *uxy;
int m1,m2;
{


  int    i,j,n,k,ni,nj,**nn,kmax,kbins1,kbins2,kbins3,kbins4,kbins;
  float  h,max,min,*term,step1,step2;
  float  hh,hx,hy,hygx,hxgy,uygx,uxgy;
  int    **imatrix();
  int    SelectNbins();
  float  *fvector();
  void   free_imatrix(),free_fvector(),cntab2();

  n=m1+m2+1;
  term=fvector(0,n);
  for(i = 0; i < n; i++)term[i]=0.0;

  j=1;
  for (i = 1; i <= m1; i++) {term[j]= data1[i];j++;}
  for (i = 1; i <= m2; i++) {term[j]= data2[i];j++;}
  max=min=term[1];
  for (i = 1; i < j ; i++ ) { if (term[i]>max)max=term[i]; if (term[i]<min)min=term[i];}
  kmax=MAXHPDF;
  //printf("Data1\n");
  kbins1= SelectNbins(data1,1,(m1+1),kmax);
  //printf("Data2\n");
  kbins2= SelectNbins(data2,1,(m2+1),kmax);
  kbins3= (kbins1+kbins2);
  //printf("Data2+Data1\n");
  kbins4= SelectNbins(term,1,j,kmax);
  if (kbins3>kbins4)kbins=kbins3;else kbins=kbins4;
  //printf("Selected Kbins=%d\n",kbins);
  if (kbins>4) kbins=(int) kbins/2;
  h=(max-min)/kbins;
  //printf("Half Kbins=%d H=%f \n",kbins,h);
  nn=imatrix(0,kbins+3,0,3);
  for (i=-1;i<kbins+1;i++){
    step1=h*i + min;
    step2=h*(i+1) + min;
    nn[i+2][1]=0;
    for (k=1;k<=m1;k++) if (data1[k]>=step1 && data1[k]<step2) nn[i+2][1]++;
    nn[i+2][2]=0;
    for (k=1;k<=m2;k++) if (data2[k]>=step1 && data2[k]<step2) nn[i+2][2]++;
   }
  ni=kbins+2;
  nj=2;
  //for (i=1;i<kbins+1;i++)printf("NN[%d][1]=%d NN[%d][2]=%d\n",i,nn[i][1],i,nn[i][2]);

  cntab2(nn,ni,nj,&hh,&hx,&hy,&hygx,&hxgy,&uygx,&uxgy,uxy); 

  free_imatrix(nn,0,kbins+3,0,3);
  free_fvector(term,0,n);

}

