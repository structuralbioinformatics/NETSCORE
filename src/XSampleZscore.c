#include "ppin.h"

void  XSampleZscore(xsZscore,xprobe,xgroup,xngroup,xdgroup,number_phe,number_samples)
express  *xprobe;
express  *xsZscore;
int      **xgroup,*xdgroup,xngroup,number_phe,number_samples[];

{
 int   i,j,k,kk,m,n;
 float meanProbe,sigmaProbe;

    for (j=0;j<number_phe;j++){
    for (n=0;n<number_samples[j];n++){
    for (i=0;i<xngroup;i++){
        meanProbe=0.0;
        sigmaProbe=0.0;
        m=0;
        for (k=0;k<xdgroup[i];k++){
        for (kk=0;kk<xprobe[xgroup[i][k]].split;kk++){
             meanProbe+=xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n];
             sigmaProbe+=xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]*xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]; 
             m++;  
             }
        if (m>0){
             meanProbe=meanProbe/m;
             sigmaProbe=sigmaProbe/m - meanProbe*meanProbe;
             }
        for (k=0;k<xdgroup[i];k++){
             strcpy(xsZscore[xgroup[i][k]].fragment[kk].name,xprobe[xgroup[i][k]].fragment[kk].name);
        if (sigmaProbe>0){
             xsZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=(xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]-meanProbe)/sqrt(sigmaProbe);
             }else if(meanProbe>0){
             xsZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=0.0;
             }else{
             xsZscore[xgroup[i][k]].fragment[kk].phenotype[j].sample[n]=xprobe[xgroup[i][k]].fragment[kk].phenotype[j].sample[n];
             }}
        }}}}

}
