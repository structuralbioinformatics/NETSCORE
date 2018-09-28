#include "ppin.h"

void  XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
 void  XIndividualGene(),XClusterGene(),XIndividualGeneMax(),XClusterGeneMax(),XRankProd();

 printf("Execute XGenePhenotype (Method[0]= %d)\n",xmethod[0]);
 switch (xmethod[0]){
    case 0:
    XIndividualGene(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
    break;
    case 1:
    XClusterGene(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
    break;
    case 2:
    XIndividualGeneMax(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xPhenoRef,xdgroup,xmethod,number_phe,number_samples);
    break;
    case 3:
    XClusterGeneMax(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
    break;
    case 4:
    XRankProd(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
    break;
    default:
    XIndividualGene(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
    }
    

}

void XIndividualGene(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
 int   i,j,jj,k,kk,n,n1,n2,m1,m2,skip,sum;
 float *data1,*data2,prob,pb,statistic,magnitude,cramrv,ccc,df,chisq;
 float *fvector(),Dmedian();
 void  ttest(), tutest(), kstwo(),BinChi(),BinCramer(),BinMI(),Bin2Cramer(),Bin2MI(); 
 void  free_fvector();


/*
   for (j=0;j<number_phe;j++){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
     xPvalue[xgroup[i][n]][j]=2.0;
     xES[xgroup[i][n]][j]=0.0;
     }}}
*/

   for (j=0;j<number_phe;j++){
   if (j==xPhenoRef || xmethod[7]==1){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
    statistic=0.0;
    magnitude=-1.0;
    prob=-1.0;
    for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
       //printf("PROB %f Mag %f \n",prob,magnitude);
       n1=number_samples[j]+1;
       sum=0;
       for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
       n2=sum+1;
       data1=fvector(0,n1);
       data2=fvector(0,n2);
       m1=0;
       for (k=0;k<number_samples[j];k++){
       if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0.0){
          m1++;
          if (m1>n1)nrerror("M1 > N1 in XIndividualGene");
          data1[m1]=xZscore[xgroup[i][n]].fragment[kk].phenotype[j].sample[k];
          }}
       m2=0;
       for (jj=0;jj<number_phe;jj++){
       if  (j!=jj){
       for (k=0;k<number_samples[jj];k++){
       if (xprobe[xgroup[i][n]].fragment[kk].phenotype[jj].defined[k]>=0.0){
          m2++;
          if (m2>n2)nrerror("M2 > N2 in XIndividualGene");
          data2[m2]=xZscore[xgroup[i][n]].fragment[kk].phenotype[jj].sample[k];
          }}}}
   
       if (m1>0 && m2>0){
      //printf("Data1 <-c(");
      //for(k=1;k<m1;k++){printf("%10.5f ,",data1[k]);}
      //printf("%10.5f )\n",data1[m1]);
      //printf("Data2 <-c(");
      //for(k=1;k<m2;k++){printf("%10.5f ,",data2[k]);}
      //printf("%10.5f )\n",data2[m2]);
        switch(xmethod[1]){ 
         case 0:
         kstwo(data1,m1,data2,m2,&statistic,&pb);
         break;
         case 1:
         tutest(data1,m1,data2,m2,&statistic,&pb);
         break;
         case 2:
         ttest(data1,m1,data2,m2,&statistic,&pb);
         break;
         case 3:
         statistic=pb=1-Dmedian(data1,m1,data2,m2);
         break;
         case 4:
         BinChi(data1,m1,data2,m2,&statistic,&pb);
         break;
         case 5:
         BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         statistic=cramrv;
         break;
         case 6:
         BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         statistic=ccc;
         break;
         case 7:
         BinMI(data1,m1,data2,m2,&statistic);
         BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         break;
         case 8:
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         statistic=1-cramrv;
         break;
         case 9:
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         statistic=1-ccc;
         break;
         case 10:
         Bin2MI(data1,m1,data2,m2,&statistic);
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
         statistic=1-statistic;
         break;
         default:
         kstwo(data1,m1,data2,m2,&statistic,&pb);
         }
        switch(xmethod[4]){
         case 0:
         if (prob<0)prob=0;
         if (magnitude<0)magnitude=0;
         if (pb>0) prob+= -log(pb)/xprobe[xgroup[i][n]].split;
         if (statistic>0) magnitude+=statistic/xprobe[xgroup[i][n]].split;
         //printf("Check Magnitude %f %f Prob %f %f\n",magnitude,statistic,prob,pb);
         break;
         case 1:
         if (prob < pb) prob=pb;
         if (magnitude>statistic || magnitude<0 ) magnitude=statistic;
         //printf("Check Magnitude %f %f Prob %f %f\n",magnitude,statistic,prob,pb);
         break;
         case 2:
         if (prob > pb || prob<0) prob=pb;
         if (magnitude < statistic) magnitude=statistic;
         break;
         default:
         if (prob<0)prob=0;
	 if (magnitude<0)magnitude=0;
         if (pb>0) prob+= -log(pb)/xprobe[xgroup[i][n]].split;
         if (statistic>0) magnitude+=statistic/xprobe[xgroup[i][n]].split;
         break;
         }
      //printf("Method %d J %d I %d N %d KK %d Size %d Statistic %f PB %15.5e Magnitude %f Probability %15.5e \n",xmethod[4],j,i,n,kk,xprobe[xgroup[i][n]].split,statistic,pb,magnitude,prob);
        }
      free_fvector(data1,0,n1);
      free_fvector(data2,0,n2);
    }
    skip=0;
    for (k=0;k<number_samples[j];k++){
    for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
    if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
         skip++;
         }}}
    if (prob<0)prob=0.0;
    if (magnitude<0) magnitude=0.0;
    if (skip>0){ xPvalue[xgroup[i][n]][j]=exp(-prob);xES[xgroup[i][n]][j]=magnitude;}
    //if (skip>0)  printf("Statistic score %f Probability %15.5e \n",magnitude,prob);
   }}}} 

}

void XClusterGene(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
 int   i,j,ii,jj,k,kk,n,n1,n2,m1,m2,skip,sum;
 float **mdata1,**mdata2,*data1,*data2,prob,statistic,cramrv,ccc,df,chisq;
 int   *ivector();
 float *fvector(),**matrix(),Dmedian();
 void  ttest(), tutest(), kstwo(),BinChi(),BinCramer(),BinMI(),Bin2Cramer(),Bin2MI(); 
 void  free_matrix(),free_fvector(),free_ivector();

   statistic=0.0;
   prob=1.0;

/*
   for (j=0;j<number_phe;j++){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
     xPvalue[xgroup[i][n]][j]=2.0;
     xES[xgroup[i][n]][j]=0.0;
     }}}
*/

   for (j=0;j<number_phe;j++){
   if  (j==xPhenoRef || xmethod[7]==1){
   for (i=0;i<xngroup;i++){
      statistic=0.0;
      prob=1.0;
      n1=MAXEST*(xdgroup[i]*number_samples[j]+1);
      sum=0;
      for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
      n2=MAXEST*(xdgroup[i]*sum+1);
      data1=fvector(0,n1);
      data2=fvector(0,n2);
      m1=0;
      for (n=0;n<xdgroup[i];n++){
      for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
      for (k=0;k<number_samples[j];k++){
      if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
          m1++;
          if (m1>n1)nrerror("M1 > N1 in XClusterGene");
          data1[m1]=xZscore[xgroup[i][n]].fragment[kk].phenotype[j].sample[k];
          }}}}
      m2=0;
      for (n=0;n<xdgroup[i];n++){
      for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
      for (jj=0;jj<number_phe;jj++){
      if  (j!=jj){
      for (k=0;k<number_samples[jj];k++){
      if (xprobe[xgroup[i][n]].fragment[kk].phenotype[jj].defined[k]>0){
          m2++;
          if (m2>n2)nrerror("M2 > N2 in XClusterGene");
          data2[m2]=xZscore[xgroup[i][n]].fragment[kk].phenotype[jj].sample[k];
          }}}}}}

      if (m1>0 && m2>0){
      //printf("Data1 <-c(");
      //for(k=1;k<m1;k++){printf("%10.5f ,",data1[k]);}
      //printf("%10.5f )\n",data1[m1]);
      //printf("Data2 <-c(");
      //for(k=1;k<m2;k++){printf("%10.5f ,",data2[k]);}
      //printf("%10.5f )\n",data2[m2]);
        switch(xmethod[1]){ 
         case 0:
         kstwo(data1,m1,data2,m2,&statistic,&prob);
         break;
         case 1:
         tutest(data1,m1,data2,m2,&statistic,&prob);
         break;
         case 2:
         ttest(data1,m1,data2,m2,&statistic,&prob);
         break;
         case 3:
         statistic=prob=1-Dmedian(data1,m1,data2,m2);
         break;
         case 4:
         BinChi(data1,m1,data2,m2,&statistic,&prob);
         break;
         case 5:
         BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         statistic=cramrv;
         break;
         case 6:
         BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         statistic=ccc;
         break;
         case 7:
         BinMI(data1,m1,data2,m2,&statistic);
         BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         break;
         case 8:
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         statistic=1-cramrv;
         break;
         case 9:
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         statistic=1-ccc;
         break;
         case 10:
         Bin2MI(data1,m1,data2,m2,&statistic);
         Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
         statistic=1-statistic;
         break;
         default:
         kstwo(data1,m1,data2,m2,&statistic,&prob);
         }
        }
      //printf("Testing Group[%d] Statistic %f Probability %15.5e \n",i, statistic,prob);
      for (n=0;n<xdgroup[i];n++){
       skip=0;
       for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
         for (k=0;k<number_samples[j];k++){
         if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
            skip++;
            }}
         }
       if (skip>0){xPvalue[xgroup[i][n]][j]=prob;xES[xgroup[i][n]][j]=statistic; }
       //if (skip>0){ printf("Check xPvalue[%d][%d]= %f xES[%d][%d]= %f \n",xgroup[i][n],j,prob,xgroup[i][n],j,statistic);}
       }
      free_fvector(data1,0,n1);
      free_fvector(data2,0,n2);
      }}} 

}

void XIndividualGeneMax(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
 int   i,j,jj,k,kk,n,n1,n2,m1,m2,skip,sum;
 float *data1,*data2,prob,pb,statistic,max1,min1,max,min,magnitude,cramrv,ccc,df,chisq;
 float *fvector(),Dmedian();
 void  ttest(), tutest(), kstwo(),BinChi(),BinCramer(),BinMI(),Bin2Cramer(),Bin2MI(); 
 void  free_fvector();

   max=0;
   min=-1;
   statistic=0.0;
   magnitude=-1.0;
   prob=-1.0;
/*
   for (j=0;j<number_phe;j++){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
     xPvalue[xgroup[i][n]][j]=2.0;
     xES[xgroup[i][n]][j]=0.0;
     }}}
*/

   for (j=0;j<number_phe;j++){
   if  (j==xPhenoRef || xmethod[7]==1){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
    for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
      max=0;
      min=-1;
      magnitude=-1.0;
      prob=-1.0;
      n1=number_samples[j]+1;
      sum=0;
      for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
      n2=sum+1;
      data1=fvector(0,n1);
      data2=fvector(0,n2);
      m1=0;
      for (k=0;k<number_samples[j];k++){
      if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0.0){
          m1++;
          data1[m1]=xZscore[xgroup[i][n]].fragment[kk].phenotype[j].sample[k];
          }}
      max1=0;
      min1=0;
      for (jj=0;jj<number_phe;jj++){
      if  (j!=jj){
        prob=-1.0;
        magnitude=-1.0;
        m2=0;
        for (k=0;k<number_samples[jj];k++){
        if (xprobe[xgroup[i][n]].fragment[kk].phenotype[jj].defined[k]>=0.0){
          m2++;
          data2[m2]=xZscore[xgroup[i][n]].fragment[kk].phenotype[jj].sample[k];
          }}
   
        if (m1>0 && m2>0){
      //printf("Data1 <-c(");
      //for(k=1;k<m1;k++){printf("%10.5f ,",data1[k]);}
      //printf("%10.5f )\n",data1[m1]);
      //printf("Data2 <-c(");
      //for(k=1;k<m2;k++){printf("%10.5f ,",data2[k]);}
      //printf("%10.5f )\n",data2[m2]);
         switch(xmethod[1]){ 
          case 0:
          kstwo(data1,m1,data2,m2,&statistic,&pb);
          break;
          case 1:
          tutest(data1,m1,data2,m2,&statistic,&pb);
          break;
          case 2:
          ttest(data1,m1,data2,m2,&statistic,&pb);
          break;
          case 3:
          statistic=pb=1-Dmedian(data1,m1,data2,m2);
          break;
          case 4:
          BinChi(data1,m1,data2,m2,&statistic,&pb);
          break;
          case 5:
          BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          statistic=cramrv;
          break;
          case 6:
          BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          statistic=ccc;
          break;
          case 7:
          BinMI(data1,m1,data2,m2,&statistic);
          BinCramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          break;
          case 8:
          Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          statistic=1-cramrv;
          break;
          case 9:
          Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          statistic=1-ccc;
          break;
          case 10:
          Bin2MI(data1,m1,data2,m2,&statistic);
          Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&pb,&cramrv,&ccc);
          statistic=1-statistic;
          break;
          default:
          kstwo(data1,m1,data2,m2,&statistic,&pb);
          }
         switch(xmethod[4]){
          case 0:
          if (prob<0)prob=0;
          if (magnitude<0)magnitude=0;
          if (pb>0) prob+= -log(pb)/xprobe[xgroup[i][n]].split;
          if (statistic>0) magnitude+=statistic/xprobe[xgroup[i][n]].split;
          break;
          case 1:
          if (prob < pb) prob=pb;
          if (magnitude>statistic || magnitude<0 ) magnitude=statistic;
          break;
          case 2:
          if (prob > pb || prob<0) prob=pb;
          if (magnitude < statistic) magnitude=statistic;
          break;
          default:
          if (prob<0)prob=0;
          if (magnitude<0)magnitude=0;
          if (pb>0) prob+= -log(pb)/xprobe[xgroup[i][n]].split;
          if (statistic>0) magnitude+=statistic/xprobe[xgroup[i][n]].split;
          }
         //printf("Check Magnitude %f %f Prob %f %f\n",magnitude,statistic,prob,pb);
         }
        if  (exp(-prob)>max1 || max1==0){max1=exp(-prob);}
        if  (magnitude < min1 || min1==0) min1=magnitude;
      //printf("Method %d J %d I %d N %d KK %d JJ %d Size %d MIN1 %f MAX1 %15.5e Magnitude %f Probability %15.5e \n",xmethod[4],j,i,n,kk,jj,xprobe[xgroup[i][n]].split,min1,max1,magnitude,prob);
       }}
      if (max1>max || max==0){max=max1;}
      if (min1<min || min<0 ){min=min1;}
      //printf("Statistic %f Probability %15.5e \n",magnitude,prob);
      //printf("Check KK %d i %d n %d J %d Max=%f Min=%f Max1=%f Min1=%f\n",kk,i,n,j,max,min,max1,min1);
      free_fvector(data1,0,n1);
      free_fvector(data2,0,n2);
      }
      skip=0;
      for (k=0;k<number_samples[j];k++){
      for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
      if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
         skip++;
         }}}
      if (skip>0){ xPvalue[xgroup[i][n]][j]=max;xES[xgroup[i][n]][j]=min;}
      //if (skip>0){ printf("Check xPvalue[%d][%d]= %f xES[%d][%d]= %f \n",xgroup[i][n],j,max,xgroup[i][n],j,min);}
      }}}} 

}

void XClusterGeneMax(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
float    **xPvalue,**xES;
express  *xprobe;
express  *xZscore;
int      **xgroup,*xdgroup,xngroup,xPhenoRef,*xmethod,number_phe,number_samples[];

{
 int   i,j,jj,k,kk,n,n1,n2,m1,m2,skip,sum;
 float *data1,*data2,prob,statistic,max1,max,min1,min,cramrv,ccc,df,chisq;
 int   *ivector();
 float *fvector(),**matrix(),Dmedian();
 void  ttest(), tutest(), kstwo(),BinChi(),BinCramer(),BinMI(),Bin2Cramer(),Bin2MI(); 
 void  free_fvector(),free_ivector(),free_matrix();

   statistic=0.0;
   prob=1.0;
/*
   for (j=0;j<number_phe;j++){
   for (i=0;i<xngroup;i++){
   for (n=0;n<xdgroup[i];n++){
     xPvalue[xgroup[i][n]][j]=2.0;
     xES[xgroup[i][n]][j]=0.0;
     }}}
*/

   for (j=0;j<number_phe;j++){
   if  (j==xPhenoRef || xmethod[7]==1){
   for (i=0;i<xngroup;i++){
      statistic=0.0;
      prob=1.0;
      n1=MAXEST*(xdgroup[i]*number_samples[j]+1);
      sum=0;
      for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
      n2=MAXEST*(xdgroup[i]*sum+1);
      data1=fvector(0,n1);
      data2=fvector(0,n2);
      m1=0;
      for (n=0;n<xdgroup[i];n++){
      for (k=0;k<number_samples[j];k++){
      for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
      if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
          m1++;
          if (m1>n1)nrerror("M1 > N1 in XClusterGeneMax");
          data1[m1]=xZscore[xgroup[i][n]].fragment[kk].phenotype[j].sample[k];
          }}}}
      max=0;
      min=-1;
      for (jj=0;jj<number_phe;jj++){
      if  (j!=jj){
       m2=0;
       for (n=0;n<xdgroup[i];n++){
       for (k=0;k<number_samples[jj];k++){
       for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
       if (xprobe[xgroup[i][n]].fragment[kk].phenotype[jj].defined[k]>0){
          m2++;
          if (m2>n2)nrerror("M2 > N2 in XClusterGeneMax");
          data2[m2]=xZscore[xgroup[i][n]].fragment[kk].phenotype[jj].sample[k];
          }}}}
       if (m1>0 && m2>0){
      //printf("Data1 <-c(");
      //for(k=1;k<m1;k++){printf("%10.5f ,",data1[k]);}
      //printf("%10.5f )\n",data1[m1]);
      //printf("Data2 <-c(");
      //for(k=1;k<m2;k++){printf("%10.5f ,",data2[k]);}
      //printf("%10.5f )\n",data2[m2]);
          switch(xmethod[1]){ 
           case 0:
           kstwo(data1,m1,data2,m2,&statistic,&prob);
           break;
           case 1:
           tutest(data1,m1,data2,m2,&statistic,&prob);
           break;
           case 2:
           ttest(data1,m1,data2,m2,&statistic,&prob);
           break;
           case 3:
           statistic=prob=1-Dmedian(data1,m1,data2,m2);
           break;
           case 4:
           BinChi(data1,m1,data2,m2,&statistic,&prob);
           break;
           case 5:
           BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           statistic=cramrv;
           break;
           case 6:
           BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           statistic=ccc;
           break;
           case 7:
           BinMI(data1,m1,data2,m2,&statistic);
           BinCramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           break;
           case 8:
           Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           statistic=1-cramrv;
           break;
           case 9:
           Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           statistic=1-ccc;
           break;
           case 10:
           Bin2MI(data1,m1,data2,m2,&statistic);
           Bin2Cramer(data1,m1,data2,m2,&chisq,&df,&prob,&cramrv,&ccc);
           statistic=1-statistic;
           break;
           default:
           kstwo(data1,m1,data2,m2,&statistic,&prob);
           }
      // printf("Check J %d JJ %d Group[%d][%d] Prob= %f Statistic=%f\n",j,jj,i,n,prob,statistic);
          }
      //if (m1>0 && m2>0){printf("Testing Group[%d] Statistic %f Probability %15.5e \n",i, statistic,prob);}
       if (prob>max || max==0){ max=prob; }
       if (statistic<min || min<0){ min=statistic; }
      //printf("Check J %d JJ %d Group[%d][%d] Max= %f Min=%f\n",j,jj,i,n,max,min);
       }}
      for (n=0;n<xdgroup[i];n++){
       skip=0;
       for (kk=0;kk<xprobe[xgroup[i][n]].split;kk++){
         for (k=0;k<number_samples[j];k++){
         if (xprobe[xgroup[i][n]].fragment[kk].phenotype[j].defined[k]>0){
            skip++;
            }}
         }
       if (skip>0){ xPvalue[xgroup[i][n]][j]=max; xES[xgroup[i][n]][j]=min;}
       //if (skip>0){ printf("Check xPvalue[%d][%d]= %f xES[%d][%d]= %f \n",xgroup[i][n],j,max,xgroup[i][n],j,min);}
       }
      free_fvector(data1,0,n1);
      free_fvector(data2,0,n2);
      }}} 

}


