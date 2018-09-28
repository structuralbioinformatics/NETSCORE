#include "ppin.h"

float  XModifyPhenoScore(iter,protein,interaction,xPhenoRef,fphe,fclus,xsimilar,xfilter,xES,xPvalue,xZscore,xPscore,xprobe,xgroup,xngroup,xdgroup,xgroup2,xngroup2,xdgroup2,xmethod,number_phe,number_samples,size,edges,clusparam,average_nodescr,rms_nodescr,average_edgescr,rms_edgescr)
int      iter;
node     *protein;
edge     *interaction;
char     fphe[MAXPHE][MAXS],fclus[MAXS];
float    **xsimilar;
int      **xfilter,xPhenoRef;
float    **xPscore;
float    **xPvalue;
float    **xES;
express  *xprobe;
express  *xZscore;
int      size,edges,**xgroup,*xdgroup,*xngroup,**xgroup2,*xdgroup2,*xngroup2,*xmethod,number_phe,number_samples[];
float    clusparam[];
float    *average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
{
 FILE  *OUT;
 FILE  *LOG;
 int   i,j,k,kk,xngrp,xngrp2,syscall,inflation,status,skip;
 pid_t pid;
 float phi,score_i,score_j,xerror,qu,qu2,squ,xflation,maximum,center,r1,r2,delta,accu;
 float *xscr;
 float distr[51];
 int   ln=1023; 
 char* mouve=malloc((ln)*sizeof(char));
 char* output=malloc((ln)*sizeof(char));
 char* cluster=malloc((ln)*sizeof(char));
 char* logfile=malloc((ln)*sizeof(char));
 char* output2=malloc((ln)*sizeof(char));
 char* cluster2=malloc((ln)*sizeof(char));
 char* logfile2=malloc((ln)*sizeof(char));
 
 int        Clustering();
 float      *fvector();
 void       XGeneZscore(),XSampleZscore(),XGenePhenotype(),XGenePhenoScore(),XAnalyzeCluster(),XAssignModifyPhenoScore(),nrerror(),free_fvector();
 
   xscr=fvector(0,size);
   xngrp=*xngroup; 
   maximum=0.0;
   qu2=qu=squ=0.0;
   for (i=0;i<size;i++)qu2+=average_nodescr[i]*average_nodescr[i];
   for (i=0;i<size;i++)qu+=average_nodescr[i];
   for (i=0;i<size;i++)if (maximum<=average_nodescr[i])maximum=average_nodescr[i];
   qu2=qu2/size;
   qu=qu/size;
   squ=qu2-qu*qu;
   if (squ>1) squ=sqrt(squ); else squ=1;
   delta=(maximum-qu)/50;
   for (i=0;i<51;i++)distr[i]=0.0;
   for (i=0;i<size;i++){
     for (j=0;j<51;j++){
       r1=qu+j*delta;
       r2=qu+(j+1)*delta;
       if (average_nodescr[i]>=r1 && average_nodescr[i]<r2)distr[j]+=1.;
       }
   }
   //printf("Distribution Scores\n");
   //for (j=0;j<51;j++) printf("DISTR[%d][%f]=%f\n",j,(j*delta+qu),distr[j]);
   accu=0.0;
   center=0.0;
   for (j=50;j>=0;j--){
      if (accu >= clusparam[13]*size/100 && center==0.0) {center=j*delta;break;}
      accu+=distr[j];
   }
   center+=qu;
   //printf("Center(accu=%f max=%f \%=%f) = QU(%f) + J(%d) X Delta(%f) = %f\n",accu, maximum, (clusparam[13]*size/100), qu, j, delta, center);
   for (i=0;i<size;i++)  xscr[i]= (average_nodescr[i] - center)/squ;
   //printf("Normalize NODES\n");
   //for (i=0;i<size;i++) printf("XSCR[%d] = ( AVE(%f) - Center(%f)) / SQU(%f) = %f\n",i,average_nodescr[i],center,squ,xscr[i]);
   for (i=0;i<size;i++){
   for (j=0;j<size;j++){
     if (xmethod[2]>0){
       if (xscr[i] > 0.0 && xscr[j] >0.0){
         phi=fabs(xsimilar[i][j]);
         }else{
         phi=xsimilar[i][j];
         }
       }else{
       if (xsimilar[i][j]!=0.0){
         phi=xsimilar[i][j];
         }else if (xscr[i] > 0.0 && xscr[j] >0.0){
         phi=MINXSIM;
         }else{
         phi=0.0;
         }
       }
     //score_i=protein[i].copy.score + xscr[i];
     //score_j=protein[j].copy.score + xscr[j];
     if (xscr[i]>0) score_i=xscr[i];else score_i=0.0;
     if (xscr[j]>0) score_j=xscr[j];else score_j=0.0;
     xsimilar[i][j]=phi*(1+sqrt(score_i*score_j));
     //if (xsimilar[i][j]>  1.0){xsimilar[i][j]=  1.0;}
     //if (xsimilar[i][j]< -1.0){xsimilar[i][j]= -1.0;}
     }}

   if (xmethod[11]<5){
     if (xmethod[11]==1){
      xflation=clusparam[6];
      inflation= (int)(100*xflation);
      sprintf(output,".network_of_similarity.I%3d.out",inflation);
      sprintf(cluster,".network_of_similarity.I%3d.clus",inflation);
      sprintf(logfile,".network_of_similarity.I%3d.log",inflation);
      sprintf(output2,".network_of_similarity.I%3d.out_%d",inflation,iter);
      sprintf(cluster2,".network_of_similarity.I%3d.clus_%d",inflation,iter);
      sprintf(logfile2,".network_of_similarity.I%3d.log_%d",inflation,iter);
      sprintf(mouve,"\\mv .network_of_similarity.mci .network_of_similarity.mci_%d",iter);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv .errors_in_MCL.out .errors_in_MCL.out_%d",iter);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv .errors_in_MCLabc.out .errors_in_MCLabc.out_%d",iter);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv .errors_in_MCXLOAD.out .errors_in_MCXLOAD.out_%d",iter);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv %s %s",output,output2);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv %s %s",cluster,cluster2);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      sprintf(mouve,"\\mv %s %s",logfile,logfile2);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      }
     if (xmethod[7]==1 || xmethod[11]==1) {
      sprintf(mouve,"\\mv .network_of_similarity.dat .network_of_similarity.dat_%d",iter);
      printf("Execute %s\n",mouve);
      fflush(stdout);
       //pid=fork();
       //if (pid<0) {fprintf(stderr,"Can't fork, error(%d) %s \n",errno,strerror(errno));nrerror("FORK() is not possible\n");}
       //if (pid==0){ syscall=system(mouve);exit(syscall);}
       //else       { wait(&status); printf("Status %d\n",WEXITSTATUS(status));}
       //if (WEXITSTATUS(status)!=0)printf("Error status of cleaning execution\n");status=0;
      syscall=system(mouve);
      printf("Done %d \n",syscall);
      printf("Matrix of Similarity\n");
      OUT=fopen(".network_of_similarity.dat","w");
      for (i=0;i<size;i++){
      for (j=i+1;j<size;j++){
        if (xmethod[7]==1 && xsimilar[i][j]>0) printf("%s %s %f ",protein[i].name1,protein[j].name1,xsimilar[i][j]);
        if (xmethod[7]==1 && xsimilar[i][j]>clusparam[12] ) printf(" *** PASS THRESHOLD  ***\n");else if (xmethod[7]==1 && xsimilar[i][j]>0) printf(" \n");
        if (clusparam[12]==0 || xsimilar[i][j] > clusparam[12]) fprintf(OUT,"%s %s %f\n",protein[i].name1,protein[j].name1,xsimilar[i][j]);
      }}
      fflush(OUT);
      fclose(OUT);
      }
     if (xmethod[11]<5) xngrp2=Clustering(protein,size,xsimilar,xfilter,size,xgroup2,xdgroup2,clusparam,xmethod[11]);
   }
   if (xmethod[11]==5) xngrp2=ReadClusterFile(protein,size,fclus,xgroup2,xdgroup2);
   if (xmethod[7]==1 || xmethod[18]==1) XAnalyzeCluster(xprobe,size,protein,xgroup2,xngrp2,xdgroup2,xgroup,xngrp,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
   XGeneZscore(xZscore,xprobe,xgroup2,xngrp2,xdgroup2,xPhenoRef,xmethod,number_phe,number_samples);
  if (xmethod[7]==1){
   for (j=0;j<number_phe;j++){
         printf("\nZscore Phenotype[%d]: %s\n",j,fphe[j]);
         for (k=0;k<number_samples[j];k++){
              printf(" => Zscoring Sample[%d]: \n",k+1);
              for (i=0;i<size;i++){
              for (kk=0;kk<xprobe[i].split;kk++){
                  if (xprobe[i].fragment[kk].phenotype[j].sample[k]!=0.0 &&
                      xZscore[i].fragment[kk].phenotype[j].sample[k]!=0.0 ) printf("\tNode[%d]Probe[%d][%d][%d]: \t%s  \tExpression: %10.3e\tZscore: %10.3e\n",i,kk,j,k,protein[i].name1,xprobe[i].fragment[kk].phenotype[j].sample[k],xZscore[i].fragment[kk].phenotype[j].sample[k]);
              }}}
         }
  }
   XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup2,xngrp2,xdgroup2,xPhenoRef,xmethod,number_phe,number_samples); 
   XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size); 

  if (xmethod[7]==1){
   printf("\n P-values:\n");
   for (j=0;j<number_phe;j++)printf(" # Phenotype[%d] is %s\n",j,fphe[j]);
   printf("\n");
   for (j=0;j<number_phe;j++){ printf("\t%15d",j);}
   printf("\n");
   for (i=0;i<size;i++){
         skip=0;
         for (j=0;j<number_phe;j++) if (xPvalue[i][j]!=0.0 && xPvalue[i][j]<=1) skip++;
         if (skip>0){
          printf("%15s\t",protein[i].name1);
          for (j=0;j<number_phe;j++){printf("%15.6e\t",xPvalue[i][j]);}
          printf("\n");
          }
         }
   printf("\n");
   printf("\n Statistic Enrichment Score:\n");
   for (j=0;j<number_phe;j++)printf(" # Phenotype[%d] is %s\n",j,fphe[j]);
   printf("\n");
   for (j=0;j<number_phe;j++){ printf("\t%15d",j);}
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
  }

   xerror=0;
   for (i=0;i<size;i++){
     if (xPvalue[i][xPhenoRef]<=1){
        xerror+=(protein[i].copy.score-xPscore[i][xPhenoRef])*(protein[i].copy.score-xPscore[i][xPhenoRef]);
        }
     }
   xerror=sqrt(xerror)/size;
   XAssignModifyPhenoScore(protein,interaction,xprobe,xsimilar,xPscore,average_nodescr,size,edges,number_phe,number_samples,xmethod,xPhenoRef);
   printf("Phenotype reference: %s\n",fphe[xPhenoRef]);
   printf("%15s\t%15s\t%15s\t%15s\n","Node","Score","xES","P-value");
   for (i=0;i<size;i++){
	 printf("%15s\t%15.6e\t%15.6e\t%15.6e\n",protein[i].name1,protein[i].copy.score,xES[i][xPhenoRef],xPvalue[i][xPhenoRef]);
         }
   fflush(stdout);


   for (i=0;i<MAXXG;i++){for (j=0;j<MAXXSG;j++){xgroup[i][j]=0;}}
   *xngroup=xngrp2;
   for (i=0;i<xngrp2;i++){
       xdgroup[i]=xdgroup2[i];
       for (j=0;j<xdgroup[i];j++){
           xgroup[i][j]=xgroup2[i][j];
           }}


   free_fvector(xscr,0,size);
   return xerror;

}
