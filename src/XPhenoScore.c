#include "ppin.h"

void  XPhenoScore(protein,interaction,xPhenoRef,fphe,xsimilar,xfilter,xES,xPvalue,xZscore,xPscore,xprobe,xgroup,xngroup,xdgroup,xmethod,number_phe,number_samples,size,edges,clusparam,fclus)
node      *protein;
edge      *interaction;
char       fphe[MAXPHE][MAXS];
float    **xsimilar;
int      **xfilter,xPhenoRef;
float    **xPscore;
float    **xPvalue;
float    **xES;
express  *xprobe;
express  *xZscore;
int      size,edges,**xgroup,*xdgroup,*xngroup,*xmethod,number_phe,*number_samples;
float      clusparam[];
char       fclus[];
{
 int        Clustering(),xngrp;
 void       XAssignPhenoScore(),XGeneZscore(),XSampleZscore(),XGenePhenotype(),XGenePhenoScore();
 int        i,j,k,kk,skip;

/*
     printf("     \t");
     for (k=0;k<size;k++){printf("%15s\t",protein[k].name1);}
     printf("\n");
     for (i=0;i<size;i++){
         printf("%5d\t",i);
         for (k=0;k<size;k++){printf("%10.3e (%2d)\t",xsimilar[i][k],xfilter[i][k]);}
         printf("\n");
         }
*/



     if (xmethod[11]<5) xngrp=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
     if (xmethod[11]==5) xngrp=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
     
     printf("\nClusters\n");
     for (i=0;i<xngrp;i++){
        printf("Group[%d] (%d):\t",i,xdgroup[i]);
        for (k=0;k<xdgroup[i];k++){
             printf("%15s ",protein[xgroup[i][k]].name1);
             }
         printf("\n");
         }
     fflush(stdout);
     //XSampleZscore(xZscore,xprobe,xgroup,xngrp,xdgroup,number_phe,number_samples);
     XGeneZscore(xZscore,xprobe,xgroup,xngrp,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
     for (j=0;j<number_phe;j++){
         if (j==xPhenoRef || xmethod[7]==1){
         printf("\nZscore Phenotype[%d]: %s\n",j,fphe[j]);
         for (k=0;k<number_samples[j];k++){
              printf(" => Zscoring Sample[%d]: \n",k+1);
              for (i=0;i<size;i++){
              for (kk=0;kk<xprobe[i].split;kk++){
                 if (xprobe[i].fragment[kk].phenotype[j].sample[k]!=0.0 && 
                     xZscore[i].fragment[kk].phenotype[j].sample[k]!=0.0) printf("\tNode[%d]Probe[%d][%d][%d]: \t%s  \tExpression: %10.3e\tZscore: %10.3e\n",i,kk,j,k,protein[i].name1,xprobe[i].fragment[kk].phenotype[j].sample[k],xZscore[i].fragment[kk].phenotype[j].sample[k]);
              }}}
         fflush(stdout);
         }}
      
     for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
     XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngrp,xdgroup,xPhenoRef,xmethod,number_phe,number_samples); 
     fflush(stdout);
     XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size); 
     fflush(stdout);
     if (xmethod[7]==1){
        printf("\n P-values:\n");
        printf("               ");
        for (j=0;j<number_phe;j++){ printf("\t%15s",fphe[j]);}
        printf("\n");
        for (i=0;i<size;i++){
            skip=0;
            for (j=0;j<number_phe;j++) if (xPvalue[i][j]!=0.0 && xPvalue[i][j]<=1.0) skip++;
            if (skip>0){
            printf("%15s\t",protein[i].name1);
            for (j=0;j<number_phe;j++){printf("%15.6e\t",xPvalue[i][j]);}
            printf("\n");
            } 
        }
        printf("\n Statistic Enrichment Score:\n");
        printf("               ");
        for (j=0;j<number_phe;j++){ printf("\t%15s",fphe[j]);}
        printf("\n");
        for (i=0;i<size;i++){
            skip=0;
            for (j=0;j<number_phe;j++) if (xES[i][j]!=0.0) skip++;
            if (skip>0){
            printf("%15s\t",protein[i].name1);
            for (j=0;j<number_phe;j++){printf("%15.6e\t",xES[i][j]);}
            printf("\n");
            }
        }
        fflush(stdout);
        printf("\n");
        printf("\n Phenotype-Scores:\n");
        for (j=0;j<number_phe;j++)printf(" # Phenotype[%d] is %s\n",j,fphe[j]);
        printf("\n               ");
        for (j=0;j<number_phe;j++){ printf("\t%15d",j);}
        printf("\n");
        for (i=0;i<size;i++){
            skip=0;
            for (j=0;j<number_phe;j++) if (xPscore[i][j]!=0.0) skip++;
            if (skip>0){
            printf("%15s\t",protein[i].name1);
            for (j=0;j<number_phe;j++){printf("%15.6e\t",xPscore[i][j]);}
            printf("\n");
            }
        }
        printf("\n");
        fflush(stdout);
     }
     XAssignPhenoScore(protein,interaction,xprobe,xsimilar,xPscore,size,edges,number_phe,number_samples,xmethod,xPhenoRef);
     fflush(stdout);
     printf("Phenotype reference: %s\n",fphe[xPhenoRef]);
     printf("%15s\t%15s\t%15s\t%15s\n","Node","Score","xES","P-value");
     for (i=0;i<size;i++){
         //if (protein[i].copy.score>0) printf("%15s\t%15.6e\t%15.6e\t%15.6e\n",protein[i].name1,protein[i].copy.score,xES[i][xPhenoRef],xPvalue[i][xPhenoRef]);
         printf("%15s\t%15.6e\t%15.6e\t%15.6e\n",protein[i].name1,protein[i].copy.score,xES[i][xPhenoRef],xPvalue[i][xPhenoRef]);
         }
     fflush(stdout);
     *xngroup=xngrp;

}
