#include "ppin.h"

float Dirac(a)
float a;
{
 float x; 
 if (a<0.05 && a>=0)  x=1.0; 
 else  x= 0.0; 
 return x;
}

void  XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size)
float    **xES;
float    **xPvalue;
float    **xPscore;
express  *xprobe;
int      *xmethod,xPhenoRef, number_phe,size;
int      *number_samples;
{

 int   i,j,k,kk,skip;
 float Dirac();

 printf("Execute XGenePhenoScore\n");
 for (i=0;i<size;i++){
 for (j=0;j<number_phe;j++){
    xPscore[i][j]=0.0;
    }}
 for (i=0;i<size;i++){
 for (j=0;j<number_phe;j++){
 if (j==xPhenoRef || xmethod[7]==1){ 
    skip=0;
    for (kk=0;kk<xprobe[i].split;kk++){
    for (k=0;k<number_samples[j];k++){
    if  (xprobe[i].fragment[kk].phenotype[j].defined[k]>0){
        skip++;
        }}}
    if (skip>0 && xPvalue[i][j]<=1.0){
         switch(xmethod[13]){
         case 0:
              xPscore[i][j]=1.0-xPvalue[i][j];
              break;
         case 1:
              xPscore[i][j]=xES[i][j];
              break;
         case 2:
              xPscore[i][j]=(1.0-xPvalue[i][j])*xES[i][j];
              break;
         case 3:
              xPscore[i][j]= Dirac(xPvalue[i][j])*xES[i][j];
              break;
         case 4:
              xPscore[i][j]= Dirac(xPvalue[i][j]);
              break;
         default:
              xPscore[i][j]=1.0-xPvalue[i][j];
         }
        }
    }}}
}

