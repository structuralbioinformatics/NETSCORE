#include "ppin.h"
#define EPS 1.0e-6

float  *NetPrioritization(edges,size2,interaction,protein,nodepval,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter)
int       **distribution,*ndstr;
int       *size2;
int       edges;
node      *protein;
edge      *interaction;
int       *n_edgescr,*n_nodescr;
float     tolerance,threshold,thero;
int       *iteration,diter;
float     maxmin[];
int       rnd[],linker[];
float    *nodepval,*nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
{
 int        i,j,kk,size,intact,iter,step,iteration_global, iteration_netscore, iteration_netzscore, iteration_netshort;      
 int        kmax,nseries,ma,mfit,ndata,kbins,k;
 int       *lista,*n_noderanscr,*n_edgeranscr;
 float     errmr, *error,*noderanscr,*average_noderanscr,*edgeranscr,*average_edgeranscr;
 float     *rms_noderanscr,*rms_edgeranscr,*param,chisq2,chisq1,dchi;
 float     *pdf,*rpdf,*mpdf,*xpdf,*testpval,**covar,**alpha,chisq,alamda, interval, max, min, origin;
 float     *NetScoreMotor(),*NetShortMotor(),*NetZScoreMotorS(),*NetZScoreMotorD(),*NetComboMotorD(),*NetComboMotorS();
 float     *fvector(),**matrix(),(*funk)(),ffgauss(),ffdgumbel(),ffgumbel(), mrqmin();
 int       SelectNbins(),*ivector();
 void      (*funcs)(), EvalueNode(),fgauss(), fgumbel(),fdgumbel(),free_matrix(),free_fvector(),free_ivector(),ScoreDistribution(),RandomScores();
 int       jmax, InitParameters();

   size=size2[0];
   intact=size2[1];

   error=fvector(0,MAXERROR);
   switch(rnd[4]){
    case 0:
     iteration_netshort=iteration[0];
     if (iteration_netshort<=0) iteration_netshort=1;
     error=NetShortMotor(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netshort,diter);
     break;
    case 1:
     iteration_netscore=iteration[0];
     if (iteration_netscore<=0) iteration_netscore=1;
     error=NetScoreMotor(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netscore,diter);
     break;
    case 2:
     iteration_netzscore=iteration[0];
     if (iteration_netzscore<=0) iteration_netzscore=1;
     error=NetZScoreMotorS(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netzscore,diter);
     break;
    case 3:
     iteration_netzscore=iteration[0];
     if (iteration_netzscore<=0) iteration_netzscore=1;
     error=NetZScoreMotorD(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netzscore,diter);
     break;
    case 4:
     error=NetComboMotorS(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
     break;
    case 5:
     error=NetComboMotorD(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
     break;
    default:
     iteration_netshort=iteration[0];
     if (iteration_netshort<=0) iteration_netshort=1;
     error=NetShortMotor(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netshort,diter);
   }
   
   if (rnd[6] <=0){                   /* P-values of nodes are irrelevant */

    for (i=0;i<size;i++){nodepval[i]=1.0;}

   }else{                  /* Random Nodes for p-value calculation */

    printf("\nRUNNING E-VALUE CALCULATION\n");
    if (rnd[6]==1) {funcs=&fgumbel;kmax=ma=mfit=ndata=3*NGUMB;}
    if (rnd[6]==2) {funcs=&fdgumbel;kmax=ma=mfit=ndata=3*NGUMB;}
    if (rnd[6]==3 || rnd[6]==4) {funcs=&fgauss;kmax=ma=mfit=ndata=3*NGAUSS;}
    //printf("RND[6] %d RND[7] %d INIT_MA %d INIT_MFIT %d \n",rnd[6],rnd[7],ma,mfit);
    if (rnd[7]>0) {mfit =3*rnd[7];}
    if (mfit>ma) mfit=ma; 

    nseries=rnd[3];
    if (MAXHPDF>kmax)kmax=MAXHPDF;
    if (kmax<MINHPDF)kmax=MINHPDF;

    param=fvector(0,ma+1);
    covar=matrix(0,ma+1,0,ma+1);
    alpha=matrix(0,ma+1,0,ma+1);
    lista=ivector(0,ma+1);
    noderanscr=fvector(0,size);
    edgeranscr=fvector(0,intact);
    average_noderanscr=fvector(0,size);
    average_edgeranscr=fvector(0,intact);
    rms_noderanscr=fvector(0,size);
    rms_edgeranscr=fvector(0,intact);
    n_noderanscr=ivector(0,size);
    n_edgeranscr=ivector(0,intact);
    for (i=0;i<size;i++){noderanscr[i]=0.0;average_noderanscr[i]=0.0;rms_noderanscr[i]=0.0;n_noderanscr[i]=0;}
    for (i=0;i<edges;i++){edgeranscr[i]=0.0;average_edgeranscr[i]=0.0;rms_edgeranscr[i]=0.0;n_edgeranscr[i]=0;}
 
    printf("\nRANDOM SCORES FOR E-VALUE CALCULATION\n");
    RandomScores(edges,size2,interaction,protein,noderanscr,average_noderanscr,rms_noderanscr,n_noderanscr,edgeranscr,average_edgeranscr,rms_edgeranscr,n_edgeranscr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
// Initialize number of bins and interval size with the first random set of scores. The same will be used on next randoms
    kbins=SelectNbins(noderanscr,0,size,kmax);
    max=min=noderanscr[0];
    for (j=0;j<size;j++) {if (max<noderanscr[j])max=noderanscr[j];if (min>noderanscr[j])min=noderanscr[j];}

    if (kbins < MINHPDF) kbins=MINHPDF;
    if (kbins >= kmax ) kmax=kbins+20;
    ndata=kmax;

    pdf=fvector(0,ndata+1);
    xpdf=fvector(0,ndata+1);
    mpdf=fvector(0,ndata+1);
    rpdf=fvector(0,ndata+1);
    testpval=fvector(0,ndata+1);

    interval=(max-min)/kbins;
    kk     = (kmax-kbins)/2;
    origin = min - kk*interval; 
    xpdf[0]= origin;
    for (k=1;k<=kmax;k++) xpdf[k]= origin + (k-1)*interval + interval/2; 
//***********************************************
    ScoreDistribution(rnd,interval,xpdf,noderanscr,size,kmax,pdf);
    for (k=0;k<=kmax;k++){mpdf[k]=0.0;rpdf[k]=0.0;}
    for (k=1;k<=kmax;k++){mpdf[k]+=pdf[k];rpdf[k]+=pdf[k]*pdf[k];}

    for (j=0;j<nseries-1;j++){
      for (i=0;i<size;i++){noderanscr[i]=0.0;average_noderanscr[i]=0.0;rms_noderanscr[i]=0.0;n_noderanscr[i]=0;}
      for (i=0;i<edges;i++){edgeranscr[i]=0.0;average_edgeranscr[i]=0.0;rms_edgeranscr[i]=0.0;n_edgeranscr[i]=0;}
      printf("\nRANDOM SCORES FOR E-VALUE CALCULATION (Iteration=%d <%d)\n",j,nseries);
      RandomScores(edges,size2,interaction,protein,noderanscr,average_noderanscr,rms_noderanscr,n_noderanscr,edgeranscr,average_edgeranscr,rms_edgeranscr,n_edgeranscr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
      ScoreDistribution(rnd,interval,xpdf,noderanscr,size,kmax,pdf);
      for (k=1;k<=kmax;k++){mpdf[k]+=pdf[k];rpdf[k]+=pdf[k]*pdf[k];}
    }

    if (nseries>1){
      for (k=1;k<=kmax;k++){
//         printf("START XPDF=%f MPDF[%d]=%f RPDF[%d]=%f\n",xpdf[k],k,mpdf[k],k,rpdf[k]);
         mpdf[k]=mpdf[k]/nseries;
         rpdf[k]=rpdf[k]/nseries -mpdf[k]*mpdf[k] ;
         if (rpdf[k]>0.0) {rpdf[k]=sqrt(rpdf[k]);}else{rpdf[k]=0.0001;}
      }
    }else{
      for (k=1;k<=kmax;k++){rpdf[k]=0.0001;}
    }
    //printf("\nRandom Scores Distribution\n");
    //for (k=1;k<=kmax;k++){printf("Bin[%5d]\t < F(x>%f) > =%e (%e) \n",k,xpdf[k],mpdf[k],rpdf[k]);}
    iter=0;
    step=0;
    dchi=1.0;
    chisq2=0;
    alamda = -1.0;
    //for (k=0;k<=ma;k++)lista[k]=k;
    //param[0]=0.0;
    //for (k=1;k<=mfit;k+=3){ param[k]=1.0; param[k+1]=0.0;param[k+2]=1.0;}
    //for (k=mfit+1;k<=ma;k+=3){ param[k]=0.0; param[k+1]=0.0; param[k+2]=1.0;}
    jmax=InitParameters(rnd,xpdf,mpdf,ndata,param,ma,lista,mfit);
    printf("Fitting of FUNCTION by Levenberg Marquard \n");
    errmr=0;
    while ( iter<100 ){
     chisq1=chisq2;
     errmr=mrqmin(xpdf,mpdf,rpdf,ndata,param,ma,lista,mfit,covar,alpha,&chisq,funcs,&alamda);
     chisq2=chisq;
     dchi=chisq2-chisq1;
     if (dchi!=0){step=0;}else{step++;}
     printf("%5d> Lambda: %e Dchisq=%e CHISQ %e (Check Zero STEP %d) \n",iter,alamda,dchi,chisq2,step);
     if (alamda<1.0e-20) alamda=-1;
     if (errmr == 0){ 
      if (sqrt(dchi*dchi)< EPS*fabs(chisq2) && dchi < 0 ) {break;}
      if (dchi >=0 && step> 4){break;}
      iter++;
     }else{
      if (mfit+3<=ma && mfit+3<=3*jmax && mfit+3<=3*rnd[8]){ 
       mfit += 3; 
       printf("Increase the number of funtions MFIT=%d\n",mfit);
       jmax=InitParameters(rnd,xpdf,mpdf,ndata,param,ma,lista,mfit);
      }else{
       printf("Iteration to fit the function has reached the limit MFIT=%d JMAX=%d NFPV2=%d \n",mfit,jmax,rnd[8]);
       break;
      }
      iter=0; 
      step=0;
      alamda = -1.0;
     }
    }
    alamda=0.0;
    errmr=mrqmin(xpdf,mpdf,rpdf,ndata,param,ma,lista,mfit,covar,alpha,&chisq,funcs,&alamda);
    if (errmr == 1) printf("\nFitting function failed\n");
    printf("\nFitting function parameters \n");
    for (k=0;k<=mfit;k++){printf("Parameter[%5d] =\t%e\n",k,param[k]);}

    printf("\nTEST OF FUNCTION PARAMETERS \n");
    if (rnd[6]==4 || rnd[6]==3) funk=&ffgauss;
    if (rnd[6]==2 ) funk=&ffdgumbel;
    if (rnd[6]==1 ) funk=&ffgumbel;
    printf("\nComparison between (CDF/PDF) binned scores and its fitted function\n");
    for (k=1;k<=kmax;k++){
      printf("Bin[%5d]\tF[%f]=%e (+/- %e) \tFitting= %e \n",k,xpdf[k],mpdf[k],rpdf[k],funk(xpdf[k],param,mfit));
     }

    for (i=0;i<=kmax;i++){testpval[i]=1.0;}
    printf("\nTest of function parameters to calculate E-values\n");
    EvalueNode(rnd,xpdf,testpval,kmax+1,param,mfit,funcs); 
    printf("\nE-values of theoretical binned Random Scores\n");
    for (k=1;k<=kmax;k++){printf("Bin[%5d]\tx>%f\t CDF/PDF= %e (+/- %e)   \tE-value= %e \n",k,xpdf[k],mpdf[k],rpdf[k],testpval[k]);}

    for (i=0;i<size;i++){nodepval[i]=1.0;}
    EvalueNode(rnd,nodescr,nodepval,size,param,mfit,funcs);
     
    free_fvector(noderanscr,0,size);
    free_fvector(edgeranscr,0,intact);
    free_fvector(average_noderanscr,0,size);
    free_fvector(average_edgeranscr,0,intact);
    free_fvector(rms_noderanscr,0,size);
    free_fvector(rms_edgeranscr,0,intact);
    free_ivector(n_noderanscr,0,size);
    free_ivector(n_edgeranscr,0,intact);
    free_fvector(param,0,ma+1);
    free_ivector(lista,0,ma+1);
    free_matrix(covar,0,ma+1,0,ma+1);
    free_matrix(alpha,0,ma+1,0,ma+1);
    free_fvector(pdf,0,ndata+1);
    free_fvector(mpdf,0,ndata+1);
    free_fvector(rpdf,0,ndata+1);
    free_fvector(xpdf,0,ndata+1);
    free_fvector(testpval,0,ndata+1);
    
   } /* End  p-value calculation with random nodes */


  return error;
}

#undef EPS
