#include "ppin.h"

void  XAssignPhenoScore(protein,interaction,xprobe,xsimilar,xPscore,size,edges,number_phe,number_samples,xmethod,nref)
node      *protein;
edge      *interaction;
express   *xprobe;
float    **xPscore,**xsimilar;
int        size,edges,number_phe,number_samples[],nref,*xmethod;
{

 int    i,j,k,kk,ii,jj,sum,n1,n2,m1,m2;
 float *data1,*data2,prob,pb,statistic;
 float *fvector();
 void   free_fvector(),kstwo(),tutest(),ttest();


 printf("PHENOTYPE SCORES ASSIGNMENT\n");
 printf("--------------------------------------------------------------------------------------------\n");
 printf("NODES\n");
 switch (xmethod[17]){
  case 0:
   for (i=0;i<size;i++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",i,protein[i].name1,protein[i].copy.score,(xPscore[i][nref]*protein[i].copy.score));
     protein[i].copy.score=xPscore[i][nref]*protein[i].copy.score;
     }
  break;
  case 1:
   for (i=0;i<size;i++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",i,protein[i].name1,protein[i].copy.score,((xPscore[i][nref]+protein[i].copy.score)));
     protein[i].copy.score=(xPscore[i][nref]+protein[i].copy.score);
     }
  break;
  case 2:
   for (i=0;i<size;i++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",i,protein[i].name1,protein[i].copy.score,((1.0+xPscore[i][nref])*protein[i].copy.score));
     protein[i].copy.score=(1.0+xPscore[i][nref])*protein[i].copy.score;
     }
  break;
  default:
   for (i=0;i<size;i++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",i,protein[i].name1,protein[i].copy.score,(xPscore[i][nref]*protein[i].copy.score));
     protein[i].copy.score=xPscore[i][nref]*protein[i].copy.score;
     }

  }

 switch(xmethod[3]){
 case 0:
       break;
 case 1:
       printf("EDGES\n");
       sum=0;
       for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
       n1=n2=sum+1;
       data1=fvector(0,n1);
       data2=fvector(0,n2);
       for (ii=0;ii<edges;ii++){
        prob=0;
        for (k=0;k<xprobe[interaction[ii].i].split;k++){
        for (kk=0;kk<xprobe[interaction[ii].j].split;kk++){
           m1=m2=0;
           statistic=0;
           for (i=0;i<number_phe;i++){
           for (j=0;j<number_samples[i];j++){
              if (xprobe[interaction[ii].i].fragment[k].phenotype[i].defined[j]>0){
                 m1++;
                 data1[m1]=xprobe[interaction[ii].i].fragment[k].phenotype[i].sample[j];
                 }
              if (xprobe[interaction[ii].j].fragment[kk].phenotype[i].defined[j]>0){
                 m2++;
                 data2[m2]=xprobe[interaction[ii].j].fragment[kk].phenotype[i].sample[j];
                 }
              }}
           if (m1>0 && m2>0){
            kstwo(data1,m1,data2,m2,&statistic,&pb);
            switch(xmethod[4]){
             case 0:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
               break;
             case 1:
               if (prob < pb) prob=pb;
               break;
             case 2:
               if (prob > pb) prob=pb;
               break;
             default:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
             }
            }
          }}
	if (xmethod[4]!=1 && xmethod[4]!=2)prob=exp(-prob);
        prob=1-prob;
        if (prob<0)prob=1.0;
        printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
        interaction[ii].association=prob*interaction[ii].association;
        }
       free_fvector(data1,0,n1);
       free_fvector(data2,0,n2);
       break;
 case 2:
       printf("EDGES\n");
       sum=0;
       for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
       n1=n2=sum+1;
       data1=fvector(0,n1);
       data2=fvector(0,n2);
       for (ii=0;ii<edges;ii++){
        prob=0.0;
        for (k=0;k<xprobe[interaction[ii].i].split;k++){
        for (kk=0;kk<xprobe[interaction[ii].j].split;kk++){
           m1=m2=0;
           statistic=0;
           for (i=0;i<number_phe;i++){
           for (j=0;j<number_samples[i];j++){
              if (xprobe[interaction[ii].i].fragment[k].phenotype[i].defined[j]>0){
                 m1++;
                 data1[m1]=xprobe[interaction[ii].i].fragment[k].phenotype[i].sample[j];
                 }
              if (xprobe[interaction[ii].j].fragment[kk].phenotype[i].defined[j]>0){
                 m2++;
                 data2[m2]=xprobe[interaction[ii].j].fragment[kk].phenotype[i].sample[j];
                 }
              }}
           if (m1>0 && m2>0){
            tutest(data1,m1,data2,m2,&statistic,&pb);
            switch(xmethod[4]){
             case 0:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
               break;
             case 1:
               if (prob < pb) prob=pb;
               break;
             case 2:
               if (prob > pb) prob=pb;
               break;
             default:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
             }
            }
          }}
	if (xmethod[4]!=1 && xmethod[4]!=2)prob=exp(-prob);
        prob=1-prob;
        if (prob<0)prob=1.0;
        printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
        interaction[ii].association=prob*interaction[ii].association;
        }
       free_fvector(data1,0,n1);
       free_fvector(data2,0,n2);
       break;
 case 3:
       printf("EDGES\n");
       sum=0;
       for (jj=0;jj<number_phe;jj++){if (j!=jj){sum+=number_samples[jj];}}
       n1=n2=sum+1;
       data1=fvector(0,n1);
       data2=fvector(0,n2);
       for (ii=0;ii<edges;ii++){
        prob=0.0;
        for (k=0;k<xprobe[interaction[ii].i].split;k++){
        for (kk=0;kk<xprobe[interaction[ii].j].split;kk++){
           m1=m2=0;
           statistic=0;
           for (i=0;i<number_phe;i++){
           for (j=0;j<number_samples[i];j++){
              if (xprobe[interaction[ii].i].fragment[k].phenotype[i].defined[j]>0){
                 m1++;
                 data1[m1]=xprobe[interaction[ii].i].fragment[k].phenotype[i].sample[j];
                 }
              if (xprobe[interaction[ii].j].fragment[kk].phenotype[i].defined[j]>0){
                 m2++;
                 data2[m2]=xprobe[interaction[ii].j].fragment[kk].phenotype[i].sample[j];
                 }
              }}
           if (m1>0 && m2>0){
            ttest(data1,m1,data2,m2,&statistic,&prob);
            switch(xmethod[4]){
             case 0:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
               break;
             case 1:
               if (prob < pb) prob=pb;
               break;
             case 2:
               if (prob > pb) prob=pb;
               break;
             default:
               if (pb>0)prob+= -log(pb)/(xprobe[interaction[ii].i].split*xprobe[interaction[ii].j].split);
             }
            }
           }}
	if (xmethod[4]!=1 && xmethod[4]!=2)prob=exp(-prob);
        prob=1-prob;
        if (prob<0)prob=1.0;
        printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
        interaction[ii].association=prob*interaction[ii].association;
        }
       free_fvector(data1,0,n1);
       free_fvector(data2,0,n2);
       break;
 case 4:
       printf("EDGES\n");
       for (ii=0;ii<edges;ii++){
           prob=fabs(xsimilar[interaction[ii].i][interaction[ii].j]);
           printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
           interaction[ii].association=prob*interaction[ii].association;
           }
       break;
 case 5:
       printf("EDGES\n");
       for (ii=0;ii<edges;ii++){
           prob=0.5*(1 + xsimilar[interaction[ii].i][interaction[ii].j]);
           printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
           interaction[ii].association=prob*interaction[ii].association;
           }
       break;
 case 6:
       printf("EDGES\n");
       // xsimilar changes sign in XSimilarityMatrix
       for (ii=0;ii<edges;ii++){
           prob=0.5*(1 + xsimilar[interaction[ii].i][interaction[ii].j]);
           printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",ii,interaction[ii].a.name1,interaction[ii].b.name1,interaction[ii].association,(prob*interaction[ii].association));
           interaction[ii].association=prob*interaction[ii].association;
           }
       break;
 default:
       break;
 }
 printf("--------------------------------------------------------------------------------------------\n");

}
