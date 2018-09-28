#include "ppin.h"

void   XAnalyzePhenoNet(OUT,protein,interaction,size,edges,threshold,xprobe,number_phe,number_samples,xPhenoRef,shortest_path,clusparam,xmethod,fclus)
FILE     *OUT;
node     *protein;
edge     *interaction;
int      size,edges,number_phe,*number_samples,*xmethod,xPhenoRef;
express  *xprobe;
int     **shortest_path;
float    *clusparam,threshold;
char      fclus[MAXS];
{
 express  *xZscore;
 float    **xsimilar,**xPscore,**xPvalue,**xES,**ScoreTest,**ScoreEdge,**SimilarExpressionKS,**SimilarExpressionTU,**SimilarExpressionTT,**CorrelationABS,**CorrelationPOS,**CorrelationNEG,**PvalueKS,**PvalueTU,**PvalueTT,**PvalueDM,**PvalueCH,*data1,*data2,*bins1,*bins2,freedom,statistic,prob,split;
 int      **xfilter,**xgroup,*xdgroup,xngroup,*xmethod2;
 int        i,j,jj,k,ii,n,sum,n1,n2,m1,m2,kk,skip,skip2,nbins;
 char       dist;
 express   *xvector();
 float    **matrix(),*fvector();
 int      **imatrix(),*ivector(),Clustering(),**XFilterMatrix();
 float    **XSimilarityMatrix(),PearsonScoreMatrix(),PearsonScoreMatrix2i(),PearsonScore(),Dmedian();
 void       XGeneZscore(),XSampleZscore(),XGenePhenotype(),XGenePhenoScore(),PrintGaussian(),kstwo(),tutest(),ttest(),chstwo();
 void       free_imatrix(),free_ivector(),free_matrix(),free_fvector(),free_xvector(),nrerror();

 xZscore=xvector(0,size);
 xsimilar=matrix(0,size,0,size);
 xfilter=imatrix(0,size,0,size);
 xgroup=imatrix(0,MAXXG,0,MAXXSG);
 xdgroup=ivector(0,MAXXG);
 xPvalue=matrix(0,size,0,MAXPHE);
 xES=matrix(0,size,0,MAXPHE);
 xPscore=matrix(0,size,0,MAXPHE);
 xmethod2=ivector(0,MAXMETHOD);
 SimilarExpressionKS=matrix(0,size,0,size);
 SimilarExpressionTU=matrix(0,size,0,size);
 SimilarExpressionTT=matrix(0,size,0,size);
 CorrelationABS=matrix(0,size,0,size);
 CorrelationPOS=matrix(0,size,0,size);
 CorrelationNEG=matrix(0,size,0,size);
 ScoreEdge=matrix(0,size,0,size);
 ScoreTest=matrix(0,5,0,size);

  switch(xmethod[5])
  { case 7: dist='e';
    case 8: dist='b';
    case 2: dist='c';
    case 4: dist='a';
    case 1: dist='u';
    case 3: dist='x';
    case 5: dist='s';
    case 6: dist='k';
    default: dist='c';
  }


 if (size > 100 ) skip=(int)size/MAXSTPS; else skip=1;
 fprintf(OUT,"\nANALYSES OF NETWORK AND EXPRESSION \n");
 fprintf(OUT,"--------------------------------------------------------------------------------------------\n");

 // Compare shortest_path and expression of nodes
       for (ii=0;ii<size;ii++)for (jj=0;jj<size;jj++)xfilter[ii][jj]=0;
       sum=0;
       for (jj=0;jj<number_phe;jj++){sum+=number_samples[jj];}
       n1=n2=sum+1;
       data1=fvector(0,n1);
       data2=fvector(0,n2);
       for (ii=0;ii<size;ii++){
           for (jj=0;jj<size;jj++){
               SimilarExpressionKS[ii][jj]=SimilarExpressionTU[ii][jj]=SimilarExpressionTT[ii][jj]=0.0;}}
       for (ii=0;ii<size;ii++){
       for (jj=ii;jj<size;jj++){
        if (xprobe[ii].split<1 || xprobe[jj].split<1)xfilter[jj][ii]=xfilter[ii][jj]=1;
       for (k=0;k<xprobe[ii].split;k++){
       for (kk=0;kk<xprobe[jj].split;kk++){
           m1=m2=0;
           statistic=0;
           prob=2.0;
           for (i=0;i<number_phe;i++){
           for (j=0;j<number_samples[i];j++){
              if (xprobe[ii].fragment[k].phenotype[i].defined[j]>0){
                 m1++;
                 if (m1>n1)nrerror("Inconsistent N1<M1 in XAnalysePhenoNet\n");
                 data1[m1]=xprobe[ii].fragment[k].phenotype[i].sample[j];
                 }
              if (xprobe[jj].fragment[kk].phenotype[i].defined[j]>0){
                 m2++;
                 if (m2>n2)nrerror("Inconsistent N2<M2 in XAnalysePhenoNet\n");
                 data2[m2]=xprobe[jj].fragment[kk].phenotype[i].sample[j];
                 }
              }}

           if (m1>0 && m2>0){kstwo(data1,m1,data2,m2,&statistic,&prob);}
           switch(xmethod[4]){
           case 0:
               if (prob>0) SimilarExpressionKS[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0) SimilarExpressionKS[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               break;
           case 1:
               if (SimilarExpressionKS[jj][ii] < prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionKS[ii][jj] < prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           case 2:
               if (SimilarExpressionKS[jj][ii] > prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionKS[ii][jj] > prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           default:
               if (prob>0) SimilarExpressionKS[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0) SimilarExpressionKS[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
           }

           if (m1>0 && m2>0){tutest(data1,m1,data2,m2,&statistic,&prob);}
           switch(xmethod[4]){
           case 0:
               if (prob>0)SimilarExpressionTU[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0)SimilarExpressionTU[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               break;
           case 1:
               if (SimilarExpressionTU[jj][ii] < prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionTU[ii][jj] < prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           case 2:
               if (SimilarExpressionTU[jj][ii] > prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionTU[ii][jj] > prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           default:
               if (prob>0)SimilarExpressionTU[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0)SimilarExpressionTU[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
           }

           if (m1>0 && m2>0){ttest(data1,m1,data2,m2,&statistic,&prob);}
           switch(xmethod[4]){
           case 0:
               if (prob>0)SimilarExpressionTT[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0)SimilarExpressionTT[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               break;
           case 1:
               if (SimilarExpressionTT[jj][ii] < prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionTT[ii][jj] < prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           case 2:
               if (SimilarExpressionTT[jj][ii] > prob) SimilarExpressionKS[jj][ii]=prob;
               if (SimilarExpressionTT[ii][jj] > prob) SimilarExpressionKS[ii][jj]=prob;
               break;
           default:
               if (prob>0)SimilarExpressionTT[jj][ii]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
               if (prob>0)SimilarExpressionTT[ii][jj]+= -log(prob)/(xprobe[ii].split*xprobe[jj].split);
           }

           if (m1<1 || m2<1)xfilter[jj][ii]=xfilter[ii][jj]=1;

           }}}}

       if (xmethod[4]!=1 && xmethod[4]!=2){
       for (ii=0;ii<size;ii++){
       for (jj=ii;jj<size;jj++){       
		SimilarExpressionKS[jj][ii]=exp(-SimilarExpressionKS[jj][ii]);
		SimilarExpressionKS[ii][jj]=exp(-SimilarExpressionKS[ii][jj]);
		SimilarExpressionTU[jj][ii]=exp(-SimilarExpressionTU[jj][ii]);
		SimilarExpressionTU[ii][jj]=exp(-SimilarExpressionTU[ii][jj]);
		SimilarExpressionTT[jj][ii]=exp(-SimilarExpressionTT[jj][ii]);
		SimilarExpressionTT[ii][jj]=exp(-SimilarExpressionTT[ii][jj]);
	   }}}

       free_fvector(data1,0,n1);
       free_fvector(data2,0,n2);

       xmethod2[0]=0;
       xmethod2[1]=0;
       xmethod2[2]=0;
       xmethod2[3]=0;
       for (ii=4;ii<MAXMETHOD;ii++){xmethod2[ii]=xmethod[ii];}
       xsimilar=XSimilarityMatrix(protein,xprobe,size,number_phe,number_samples,xmethod2,clusparam,shortest_path,dist);  
       for (ii=0;ii<size;ii++){
       for (jj=ii;jj<size;jj++){
           CorrelationABS[jj][ii]=CorrelationABS[ii][jj]=fabs(xsimilar[ii][jj]);
           CorrelationPOS[jj][ii]=CorrelationPOS[ii][jj]=0.5*(1 + xsimilar[ii][jj]);
           CorrelationNEG[jj][ii]=CorrelationNEG[ii][jj]=0.5*(1 - xsimilar[ii][jj]);
           }}
       fprintf(OUT,"\nComparison of Shortest_path and Expression Correlations\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Node1","Node2","Short_Path","Sim_KS","Sim_TU","Sim_TT","PCE_abs","PCE_pos","PCE_neg");
       kk=0;
       for (ii=0;ii<size;ii++){
       for (jj=ii+1;jj<size;jj++){
           if ( kk < MAXSTPS  ){
             if ( xfilter[ii][jj]==0) 
               fprintf(OUT,"%10s\t%10s\t%10d\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\n",protein[ii].name1,protein[jj].name1,shortest_path[ii][jj],SimilarExpressionKS[ii][jj], SimilarExpressionTU[ii][jj],SimilarExpressionTT[ii][jj],CorrelationABS[ii][jj],CorrelationPOS[ii][jj],CorrelationNEG[ii][jj]);
           }else{ 
             if ( kk == skip*((int)kk/skip) && xfilter[ii][jj]==0 )
               fprintf(OUT,"%10s\t%10s\t%10d\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\n",protein[ii].name1,protein[jj].name1,shortest_path[ii][jj],SimilarExpressionKS[ii][jj], SimilarExpressionTU[ii][jj],SimilarExpressionTT[ii][jj],CorrelationABS[ii][jj],CorrelationPOS[ii][jj],CorrelationNEG[ii][jj]);
           }
           kk++;
           }}
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"Correlation Shortest Path vs Sim_KS: \t%15.7e\n",PearsonScoreMatrix2i(SimilarExpressionKS,shortest_path,xfilter,size,size));    
       fprintf(OUT,"Correlation Shortest Path vs Sim_TU: \t%15.7e\n",PearsonScoreMatrix2i(SimilarExpressionTU,shortest_path,xfilter,size,size));    
       fprintf(OUT,"Correlation Shortest Path vs Sim_TT: \t%15.7e\n",PearsonScoreMatrix2i(SimilarExpressionTT,shortest_path,xfilter,size,size));    
       fprintf(OUT,"Correlation Shortest Path vs PCE_abs:\t%15.7e\n",PearsonScoreMatrix2i(CorrelationABS,shortest_path,xfilter,size,size));    
       fprintf(OUT,"Correlation Shortest Path vs PCE_pos:\t%15.7e\n",PearsonScoreMatrix2i(CorrelationPOS,shortest_path,xfilter,size,size));    
       fprintf(OUT,"Correlation Shortest Path vs PCE_neg:\t%15.7e\n",PearsonScoreMatrix2i(CorrelationNEG,shortest_path,xfilter,size,size));    
       fflush(OUT);
 // Calculate PhenoScores with/without PPI and/or Similar Expression levels


       //Without Interactions
       n=0;
       for (i=0;i<size;i++){for (j=0;j<size;j++){ xfilter[i][j]=xfilter[j][i]=1;}}
       if (xmethod[11]<5) xngroup=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
       if (xmethod[11]==5) xngroup=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
       XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
       XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);
       for (i=0;i<size;i++)ScoreTest[n][i]=xPscore[i][xPhenoRef];

       fprintf(OUT,"\nClusters Without Interactions\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for (i=0;i<xngroup;i++){
        fprintf(OUT,"Group[%d]:\t",i);
        for (k=0;k<xdgroup[i];k++){
             fprintf(OUT,"%15s ",protein[xgroup[i][k]].name1);
             }
         fprintf(OUT,"\n");
         }
       fflush(OUT);

       //With Interactions
       n=1;
       xfilter=XFilterMatrix(protein,interaction,size,threshold);
       if (xmethod[11]<5) xngroup=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
       if (xmethod[11]==5) xngroup=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
       XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
       XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);
       for (i=0;i<size;i++)ScoreTest[n][i]=xPscore[i][xPhenoRef];

       fprintf(OUT,"\nClusters With Interactions\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for (i=0;i<xngroup;i++){
        fprintf(OUT,"Group[%d]:\t",i);
        for (k=0;k<xdgroup[i];k++){
             fprintf(OUT,"%15s ",protein[xgroup[i][k]].name1);
             }
         fprintf(OUT,"\n");
         }
       fflush(OUT);

       //With Interactions and Similar Expression KS
       n=2;
       xfilter=XFilterMatrix(protein,interaction,size,threshold);
       for (ii=0;ii<edges;ii++)
           if (SimilarExpressionKS[interaction[ii].i][interaction[ii].j]<0.5)xfilter[interaction[ii].i][interaction[ii].j]=xfilter[interaction[ii].j][interaction[ii].i]=0;
       if (xmethod[11]<5) xngroup=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
       if (xmethod[11]==5) xngroup=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
       XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
       XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);
       for (i=0;i<size;i++)ScoreTest[n][i]=xPscore[i][xPhenoRef];

       fprintf(OUT,"\nClusters With Interactions and Similar Expression according to Kolmogorov Smirnov (KS)\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for (i=0;i<xngroup;i++){
        fprintf(OUT,"Group[%d]:\t",i);
        for (k=0;k<xdgroup[i];k++){
             fprintf(OUT,"%15s ",protein[xgroup[i][k]].name1);
             }
         fprintf(OUT,"\n");
         }
       fflush(OUT);
       //With Interactions and Similar Expression TU
       n=3;
       xfilter=XFilterMatrix(protein,interaction,size,threshold);
       for (ii=0;ii<edges;ii++)
           if (SimilarExpressionTU[interaction[ii].i][interaction[ii].j]<0.5)xfilter[interaction[ii].i][interaction[ii].j]=xfilter[interaction[ii].j][interaction[ii].i]=0;
       if (xmethod[11]<5) xngroup=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
       if (xmethod[11]==5) xngroup=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
       XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
       XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);
       for (i=0;i<size;i++)ScoreTest[n][i]=xPscore[i][xPhenoRef];

       fprintf(OUT,"\nClusters With Interactions and Similar Expression according to T-test unequal variance (TU)\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for (i=0;i<xngroup;i++){
        fprintf(OUT,"Group[%d]:\t",i);
        for (k=0;k<xdgroup[i];k++){
             fprintf(OUT,"%15s ",protein[xgroup[i][k]].name1);
             }
         fprintf(OUT,"\n");
         }
       fflush(OUT);
       //With Interactions and Similar Expression TT
       n=4;
       xfilter=XFilterMatrix(protein,interaction,size,threshold);
       for (ii=0;ii<edges;ii++)
           if (SimilarExpressionTT[interaction[ii].i][interaction[ii].j]<0.5)xfilter[interaction[ii].i][interaction[ii].j]=xfilter[interaction[ii].j][interaction[ii].i]=0;
       if (xmethod[11]<5) xngroup=Clustering(protein,size,xsimilar,xfilter,size,xgroup,xdgroup,clusparam,xmethod[11]);
       if (xmethod[11]==5) xngroup=ReadClusterFile(protein,size,fclus,xgroup,xdgroup);
       XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       for (i=0;i<size;i++)for (j=0;j<number_phe;j++){xPvalue[i][j]=2; xES[i][j]=0;}
       XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
       XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);
       for (i=0;i<size;i++)ScoreTest[n][i]=xPscore[i][xPhenoRef];

       fprintf(OUT,"\nClusters With Interactions and Similar Expression according to Student T-test (TT)\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for (i=0;i<xngroup;i++){
        fprintf(OUT,"Group[%d]:\t",i);
        for (k=0;k<xdgroup[i];k++){
             fprintf(OUT,"%15s ",protein[xgroup[i][k]].name1);
             }
         fprintf(OUT,"\n");
         }

       fflush(OUT);

       fprintf(OUT,"\nComparison of PhenoScores\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Node","Name","No_PPI [0]","+PPI [1]","PPI+KS [2]","PPI+TU [3]","PPI+TT [4]");
       for (i=0;i<size;i++){
          fprintf(OUT,"%10d\t%10s\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\n",i,protein[i].name1,ScoreTest[0][i],ScoreTest[1][i],ScoreTest[2][i],ScoreTest[3][i],ScoreTest[4][i]);
          }
       fprintf(OUT,"\nCorrelation Matrix between PhenoScores\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10d\t%10d\t%10d\t%10d\t%10d\n"," ",0,1,2,3,4);
       for (i=0;i<5;i++){
           fprintf(OUT,"%10d\t",i);
           for (j=0;j<4;j++){
               fprintf(OUT,"%10.5f\t",PearsonScoreMatrix(ScoreTest,size,i,j));
               }
           fprintf(OUT,"%10.5f\n",PearsonScoreMatrix(ScoreTest,size,i,4));
           }
       fprintf(OUT,"\nComparison of Quadratic Mean of PhenoScores and Shortest Path\n");
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Node1","Node2","Short_Path","No_PPI [0]","+PPI [1]","PPI+KS [2]","PPI+TU [3]","PPI+TT [4]");
       for (ii=0;ii<size;ii++){
         for(jj=ii+1;jj<size;jj++){
           skip2=0;
           for(k=0;k<5;k++){skip2+=ScoreTest[k][ii]*ScoreTest[k][jj];}
           if (skip2==0){
             xfilter[ii][jj]=1;
             xfilter[jj][ii]=1;
             }else{  
             xfilter[ii][jj]=0;
             xfilter[jj][ii]=0;
             }
           }}
       kk=0;
       for (ii=0;ii<size;ii++){
       for (jj=ii+1;jj<size;jj++){
       if ( kk == skip*((int)kk/skip) && xfilter[ii][jj]==0 ){
           fprintf(OUT,"%10s\t%10s\t%10d\t",protein[ii].name1,protein[jj].name1,shortest_path[ii][jj]);
           for(k=0;k<4;k++){
               fprintf(OUT,"%10.7f\t",sqrt(ScoreTest[k][ii]*ScoreTest[k][jj]));
               }
           fprintf(OUT,"%10.7f\n",sqrt(ScoreTest[4][ii]*ScoreTest[4][jj]));
           }
           kk++;
           }}
       fprintf(OUT,  "--------------------------------------------------------------------------------------------\n");
       for(k=0;k<5;k++){
          for (ii=0;ii<size;ii++){
             for (jj=ii;jj<size;jj++){
                 if (ScoreTest[k][ii]*ScoreTest[k][jj]==0){
                    xfilter[ii][jj]=1;
                    }else{  
                    xfilter[ii][jj]=0;
                    }}}
          for (ii=0;ii<size;ii++)
          for (jj=ii;jj<size;jj++)
              ScoreEdge[ii][jj]=ScoreEdge[jj][ii]=sqrt(ScoreTest[k][ii]*ScoreTest[k][jj]);
          fprintf(OUT,"Correlation Shortest Path vs Quadratic Mean of PhenoScore[%d]:\t%15.7e\n",k,PearsonScoreMatrix2i(ScoreEdge,shortest_path,xfilter,size,size));
          }

 // Distribution of Shortest Path distances among PhenoScores under a threshold

       fflush(OUT);

/*
       PvalueKS=matrix(0,5,0,MAXSTPS);
       PvalueTT=matrix(0,5,0,MAXSTPS);
       PvalueTU=matrix(0,5,0,MAXSTPS);
       PvalueDM=matrix(0,5,0,MAXSTPS);
       PvalueCH=matrix(0,5,0,MAXSTPS);

       for (j=0;j<5;j++){for (k=0;k<MAXSTPS;k++){PvalueKS[j][k]=PvalueTT[j][k]=PvalueTU[j][k]=PvalueDM[j][k]=PvalueCH[j][k]=0.0;}}
       for (j=0;j<5;j++){
       for (k=0;k<MAXSTPS;k++){
        split=(float)k/(float)MAXSTPS;
        n1=n2=0;
        for (ii=0;ii<size;ii++){
        for (jj=ii+1;jj<size;jj++){
        if (shortest_path[ii][jj]<size){
            if (  ScoreTest[j][ii]<=split && ScoreTest[j][jj]<=split ) n1++; else n2++; 
            }}}
        m1=n1+1;
        m2=n2+1;
        data1=fvector(0,m1);
        data2=fvector(0,m2);
        bins1=fvector(0,MAXSTPG+1);
        bins2=fvector(0,MAXSTPG+1);
        n1=n2=0;
        for (ii=0;ii<size;ii++){
        for (jj=ii+1;jj<size;jj++){
        if (shortest_path[ii][jj]<size){
            if (  ScoreTest[j][ii]<=split && ScoreTest[j][jj]<=split ){
               n1++;
               data1[n1]=(float)shortest_path[ii][jj];
               }else{
               n2++;
               data2[n2]=(float)shortest_path[ii][jj];
               }
            }}}
        fprintf(OUT,"\nDistribution of Shortest Path among PhenoScores[%d] higher than %10.5e (N1=%5d) (N2=%5d) \n",j,split,n1,n2);
        fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
        if (n1>1 && n2>1){   
         PrintGaussian(OUT,data1,1,m1,data2,1,m2,bins1,bins2,nbins);
         kstwo(data1,n1,data2,n2,&statistic,&prob);
         PvalueKS[j][k]=prob;
         ttest(data1,n1,data2,n2,&statistic,&prob);
         PvalueTT[j][k]=prob;
         tutest(data1,n1,data2,n2,&statistic,&prob);
         PvalueTU[j][k]=prob;
         statistic=prob=1-Dmedian(data1,n1,data2,n2);
         PvalueDM[j][k]=prob;
         chstwo(bins1,bins2,nbins,0,&freedom,&statistic,&prob);
         PvalueCH[j][k]=prob;
         }
         fflush(OUT);
         free_fvector(data1,0,m1);
         free_fvector(data2,0,m2);
         free_fvector(bins1,0,MAXSTPG+1);
         free_fvector(bins2,0,MAXSTPG+1);
       }}

       fprintf(OUT,"\nKolmogorov-Smirnov P-values between Shortest Path distribution low and high PhenoScores\n"); 
       fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","threshold","No_PPI","+PPI","PPI+SimKS","PPI+SimTU","PPI+SimTT");
       for (k=0;k<MAXSTPS;k++){
           fprintf(OUT,"%10.5f\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",(float)k/(float)MAXSTPS,PvalueKS[0][k],PvalueKS[1][k],PvalueKS[2][k],PvalueKS[3][k],PvalueKS[4][k]);
           }
       fprintf(OUT,"Unequal variance T-test P-values between Shortest Path distribution low and high PhenoScores \n"); 
       fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","threshold","No_PPI","+PPI","PPI+SimKS","PPI+SimTU","PPI+SimTT");
       for (k=0;k<MAXSTPS;k++){
           fprintf(OUT,"%10.5f\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",(float)k/(float)MAXSTPS,PvalueTU[0][k],PvalueTU[1][k],PvalueTU[2][k],PvalueTU[3][k],PvalueTU[4][k]);
           }
       fprintf(OUT,"Student T-test P-values between Shortest Path distribution low and high PhenoScores\n"); 
       fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","threshold","No_PPI","+PPI","PPI+SimKS","PPI+SimTU","PPI+SimTT");
       for (k=0;k<MAXSTPS;k++){
           fprintf(OUT,"%10.5f\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",(float)k/(float)MAXSTPS,PvalueTT[0][k],PvalueTT[1][k],PvalueTT[2][k],PvalueTT[3][k],PvalueTT[4][k]);
           }
       fprintf(OUT,"P-values by Ratio of Median Difference between Shortest Path distribution low and high PhenoScores\n"); 
       fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","threshold","No_PPI","+PPI","PPI+SimKS","PPI+SimTU","PPI+SimTT");
       for (k=0;k<MAXSTPS;k++){
           fprintf(OUT,"%10.5f\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",(float)k/(float)MAXSTPS,PvalueDM[0][k],PvalueDM[1][k],PvalueDM[2][k],PvalueDM[3][k],PvalueDM[4][k]);
           }
       fprintf(OUT,"P-values by Chi-Square between Shortest Path distribution low and high PhenoScores\n"); 
       fprintf(OUT,"------------------------------------------------------------------------------------------------\n");
       fprintf(OUT,"%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","threshold","No_PPI","+PPI","PPI+SimKS","PPI+SimTU","PPI+SimTT");
       for (k=0;k<MAXSTPS;k++){
           fprintf(OUT,"%10.5f\t%10.4e\t%10.4e\t%10.4e\t%10.4e\t%10.4e\n",(float)k/(float)MAXSTPS,PvalueCH[0][k],PvalueCH[1][k],PvalueCH[2][k],PvalueCH[3][k],PvalueCH[4][k]);
           }
       free_matrix(PvalueKS,0,5,0,MAXSTPS);
       free_matrix(PvalueTT,0,5,0,MAXSTPS);
       free_matrix(PvalueTU,0,5,0,MAXSTPS);
       free_matrix(PvalueDM,0,5,0,MAXSTPS);
       free_matrix(PvalueCH,0,5,0,MAXSTPS);

*/
       fflush(OUT);

 free_matrix(xsimilar,0,size,0,size);
 free_imatrix(xfilter,0,size,0,size);
 free_imatrix(xgroup,0,MAXXG,0,MAXXSG);
 free_ivector(xdgroup,0,MAXXG);
 free_matrix(xPvalue,0,size,0,MAXPHE);
 free_matrix(xES,0,size,0,MAXPHE);
 free_matrix(xPscore,0,size,0,MAXPHE);
 free_ivector(xmethod2,0,MAXMETHOD);
 free_matrix(SimilarExpressionKS,0,size,0,size);
 free_matrix(SimilarExpressionTU,0,size,0,size);
 free_matrix(SimilarExpressionTT,0,size,0,size);
 free_matrix(CorrelationABS,0,size,0,size);
 free_matrix(CorrelationPOS,0,size,0,size);
 free_matrix(CorrelationNEG,0,size,0,size);
 free_matrix(ScoreEdge,0,size,0,size);
 free_matrix(ScoreTest,0,5,0,size);
 free_xvector(xZscore,0,size);
 
}
