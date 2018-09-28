#include "ppin.h"

void   XAnalyzeCluster(xprobe,size,protein,xgroup2,xngroup2,xdgroup2,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples)
express  *xprobe;
node     *protein;
int      size,**xgroup,*xdgroup,xngroup,**xgroup2,*xdgroup2,xngroup2,xPhenoRef,*xmethod,number_phe,*number_samples;
{
 express   *xZscore,*xZscore2;
 int        i,j,k,n,m,rank,*select, **skip,*pair1,*pair2; 
 float     **xES,**xES2,**xPvalue,**xPscore,**xPvalue2,**xPscore2;
 float     *common,phscr1,phscr2,min1,min2,max1,max2,increase,decrease;

 express   *xvector();
 int       *ivector(),**imatrix();
 float     *fvector(),Common2(),**matrix();
 void       XGeneZscore(),XGenePhenotype(),XGenePhenoScore(),SortCommon();
 void       free_xvector(),free_ivector(),free_fvector();

 xZscore=xvector(0,size); 
 xZscore2=xvector(0,size); 
 xES=matrix(0,size,0,MAXPHE);
 xES2=matrix(0,size,0,MAXPHE);
 xPvalue=matrix(0,size,0,MAXPHE);
 xPvalue2=matrix(0,size,0,MAXPHE);
 xPscore=matrix(0,size,0,MAXPHE);
 xPscore2=matrix(0,size,0,MAXPHE);
 skip=imatrix(0,xngroup,0,xngroup2);
 pair1=ivector(0,xngroup2*xngroup+1);
 pair2=ivector(0,xngroup2*xngroup+1);
 select=ivector(0,xngroup2*xngroup+1);
 common=fvector(0,xngroup2*xngroup+1);

 printf("\nANALYSIS OF THE NEW CLUSTERS\n");
 printf("--------------------------------------------------------------------------------\n");
 printf("Old Clusters\n");
 for (i=0;i<xngroup;i++){
        printf("Group[%d] (%d):\t",i,xdgroup[i]);
        for (k=0;k<xdgroup[i];k++){
             if (xgroup[i][k]<size && xgroup[i][k]>=0){printf("%15s ",protein[xgroup[i][k]].name1);}
             else                                    {printf("Error XGROUP[%d][%d]= %d SIZE %d\n",i,k,xgroup[i][k],size);}
             }
         printf("\n");
         }
 printf("\nNew Clusters\n");
 for (i=0;i<xngroup2;i++){
        printf("Group[%d] (%d):\t",i,xdgroup2[i]);
        for (k=0;k<xdgroup2[i];k++){
             if (xgroup2[i][k]<size && xgroup2[i][k]>=0){printf("%15s ",protein[xgroup2[i][k]].name1);}
             else                                      {printf("Error XGROUP2[%d][%d]= %d SIZE %d\n",i,k,xgroup2[i][k],size);}
             }
         printf("\n");
         }
 printf("Assign Zscore with new clusters\n");
 XGeneZscore(xZscore2,xprobe,xgroup2,xngroup2,xdgroup2,xPhenoRef,xmethod,number_phe,number_samples);
 printf("Assign Phenotype with new clusters\n");
 XGenePhenotype(xES2,xPvalue2,xZscore2,xprobe,xgroup2,xngroup2,xdgroup2,xPhenoRef,xmethod,number_phe,number_samples); 
 printf("Calculate Phenoscore with new clusters\n");
 XGenePhenoScore(xES2,xPvalue2,xPscore2,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);

 printf("Assign Zscore with old clusters\n");
 XGeneZscore(xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples);
 printf("Assign Phenotype with old clusters\n");
 XGenePhenotype(xES,xPvalue,xZscore,xprobe,xgroup,xngroup,xdgroup,xPhenoRef,xmethod,number_phe,number_samples); 
 printf("Calculate Phenoscore with old clusters\n");
 XGenePhenoScore(xES,xPvalue,xPscore,xprobe,xPhenoRef,xmethod,number_phe,number_samples,size);

 printf("Check Common clusters\n");
 n=0;
 rank=0;
 for (i=0;i<xngroup;i++){
 for (j=0;j<xngroup2;j++){
  skip[i][j]=1;
  if (xdgroup[i]>1 || xdgroup2[j]>1){
    n++;
    common[n]=Common2(xgroup,xgroup2,xdgroup,xdgroup2,i,j);
    pair1[n]=i;
    pair2[n]=j;
    if (common[n]>0) skip[i][j]=0;
    //printf("common[%d]=%f IJ[%d][%d]\n",n,common[n],pair1[n],pair2[n]);
    }}}
 if (n>1){
  SortCommon(n,common,pair1,pair2);  
  for (i=n;i>0;i--){
   //printf("Reorder common[%d]=%f IJ[%d][%d]\n",i,common[i],pair1[i],pair2[i]);
   if (skip[pair1[i]][pair2[i]]==0){
       for (j=0;j<xngroup2;j++){skip[pair1[i]][j]=1;}
       for (j=0;j<xngroup;j++) {skip[j][pair2[i]]=1;}
       skip[pair1[i]][pair2[i]]=0;
      }
   } 

  for (i=n;i>0;i--){
  if (skip[pair1[i]][pair2[i]]==0){
   rank++;
   select[rank]=i;
   }}
 }
 printf("\nImprovement of the PhenoScores between similar clusters (with common nodes)\n"); 
 printf("--------------------------------------------------------------------------------\n");
 printf("\n%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n","Rank","#Common","Cluster[0]","Cluster[1]","<PhScr[0]>","<PhScr[1]>","Dif<PhScr>" ,"minPval[0]","minPval[1]" ,"maxPval[0]","maxPval[1]");
 increase=0.0;
 decrease=0.0;
 for (i=1;i<=rank;i++){
   k=select[i];
   phscr1=phscr2=min1=max1=min2=max2=0;
   //printf("\nGroup_old[%d] => ",pair1[k]);
   for (j=0;j<xdgroup[pair1[k]];j++){
      //printf(" Add[%s] %10.5f ",protein[xgroup[pair1[k]][j]].name1,xPscore[xgroup[pair1[k]][j]][xPhenoRef]);
      phscr1+=xPscore[xgroup[pair1[k]][j]][xPhenoRef];
      if (j==0 || min1>xPvalue[xgroup[pair1[k]][j]][xPhenoRef])min1=xPvalue[xgroup[pair1[k]][j]][xPhenoRef];
      if (j==0 || max1<xPvalue[xgroup[pair1[k]][j]][xPhenoRef])max1=xPvalue[xgroup[pair1[k]][j]][xPhenoRef];
      }
   //printf("\nGroup_new[%d] => ",pair2[k]);
   for (j=0;j<xdgroup2[pair2[k]];j++){
      //printf(" Add[%s] %10.5f ",protein[xgroup2[pair2[k]][j]].name1,xPscore2[xgroup2[pair2[k]][j]][xPhenoRef]);
      phscr2+=xPscore2[xgroup2[pair2[k]][j]][xPhenoRef];
      if (j==0 || min2>xPvalue2[xgroup2[pair2[k]][j]][xPhenoRef])min2=xPvalue2[xgroup2[pair2[k]][j]][xPhenoRef];
      if (j==0 || max2<xPvalue2[xgroup2[pair2[k]][j]][xPhenoRef])max2=xPvalue2[xgroup2[pair2[k]][j]][xPhenoRef];
      }
   //printf("\n");
   phscr1=phscr1/xdgroup[pair1[k]];
   phscr2=phscr2/xdgroup2[pair2[k]];
   printf("%10d\t%10.5f\t%10d\t%10d\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n",i,common[k],pair1[k],pair2[k],phscr1,phscr2,phscr2-phscr1,min1,min2,max1,max2);
   if (phscr2>=phscr1) increase+=phscr2-phscr1;else decrease+=phscr2-phscr1;
   } 
 printf("--------------------------------------------------------------------------------\n");
 printf("Total increase of <PhenoScores>: %10.5f\n",increase);
 printf("Total decrease of <PhenoScores>: %10.5f\n",decrease);


 free_imatrix(skip,0,xngroup,0,xngroup2);
 free_ivector(pair1,0,xngroup2*xngroup+1);
 free_ivector(pair2,0,xngroup2*xngroup+1);
 free_ivector(select,0,xngroup2*xngroup+1);
 free_xvector(xZscore,0,size);
 free_xvector(xZscore2,0,size);
 free_matrix(xES,0,size,0,MAXPHE);
 free_matrix(xES2,0,size,0,MAXPHE);
 free_matrix(xPvalue,0,size,0,MAXPHE);
 free_matrix(xPvalue2,0,size,0,MAXPHE);
 free_matrix(xPscore,0,size,0,MAXPHE); 
 free_matrix(xPscore2,0,size,0,MAXPHE); 

}


float  Common2(g1,g2,d1,d2,i1,i2)
int   **g1,**g2;
int     i1,i2;
int    *d1,*d2;
{
 int   j,k;
 float x;
 x=0.0;
 //printf("Check Common[%d][%d]\n",i1,i2);
 for (k=0;k<d1[i1];k++){
 for (j=0;j<d2[i2];j++){
   if (g1[i1][k]==g2[i2][j])x++; 
 }}
 return  x;
}

void SortCommon(n,ra,ir1,ir2)
int n,*ir1,*ir2;
float ra[];
{
        int l,j,ir,i,ira,irb;
        float rra;

        //printf("Sort Common %d clusters\n",n);

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1){
                        rra=ra[--l];
                        ira=ir1[l];
                        irb=ir2[l];
                } else {
                        rra=ra[ir];
                        ira=ir1[ir];
                        irb=ir2[ir];
                        ra[ir]=ra[1];
                        ir1[ir]=ir1[1];
                        ir2[ir]=ir2[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                ir1[1]=ira;
                                ir2[1]=irb;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                ir1[i]=ir1[j];
                                ir2[i]=ir2[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
                ir1[i]=ira;
                ir2[i]=irb;
        }
}

