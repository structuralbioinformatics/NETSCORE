#include "ppin.h"

void  XGeneZscore(xgZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
express  *xprobe;
express  *xgZscore;
int      **xgroup,*xdgroup,*xmethod,xngroup,xPhenoRef,number_phe,*number_samples;
{
 int   i,j,k,kk,m,n,sign,positive,negative;
 float meanProbe,sigmaProbe,min;

  printf("Execute XGeneZscore\n");
    for (i=0;i<xngroup;i++){
     for (k=0;k<xdgroup[i];k++){
     for (j=0;j<number_phe;j++){
     for (n=0;n<number_samples[j];n++){
     for (kk=0;kk<xprobe[xgroup[i][k]].split;kk++){
            xgZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=0.0;}}}}
     for (k=0;k<xdgroup[i];k++){
     for (kk=0;kk<xprobe[xgroup[i][k]].split;kk++){
        meanProbe=0.0;
        sigmaProbe=0.0;
        m=0;
        for (j=0;j<number_phe;j++){
        for (n=0;n<number_samples[j];n++){
        if (xprobe[xgroup[i][k]].fragment[kk].phenotype[j].defined[n]==1){
             meanProbe+=xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n];
             sigmaProbe+=xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]*xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]; 
             m++;  
             }}}
        if (m>0){
             meanProbe=meanProbe/m;
             sigmaProbe=sigmaProbe/m - meanProbe*meanProbe;
             }
        sign=1;
        if (xmethod[2]==2){
           positive=0;
           negative=0;
           for (n=0;n<number_samples[xPhenoRef];n++){
           if ((xprobe[xgroup[i][k]].fragment[kk].phenotype[xPhenoRef].sample[n]-meanProbe)>0){
              positive++;
              }else{
              negative++;
              }}
           if (positive!=negative)sign=(positive-negative)/fabs(positive-negative);
           }
        for (j=0;j<number_phe;j++){
        for (n=0;n<number_samples[j];n++){
             strcpy(xgZscore[xgroup[i][k]].fragment[kk].name,xprobe[xgroup[i][k]].fragment[kk].name);
        if (sigmaProbe>0 && xprobe[xgroup[i][k]].fragment[kk].phenotype[j].defined[n]==1){
             xgZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]+=sign*(xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]-meanProbe)/sqrt(sigmaProbe);
             }else{
             xgZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=0.0;
             }}}
        }}
     for (k=0;k<xdgroup[i];k++){
     for (kk=0;kk<xprobe[xgroup[i][k]].split;kk++){
        for (j=0;j<number_phe;j++){
        for (n=0;n<number_samples[j];n++){
           if (xdgroup[i]>0)
           xgZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=xgZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]/sqrt(xdgroup[i]);
        }}}}
    }

}
