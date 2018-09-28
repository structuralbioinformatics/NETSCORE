#include "ppin.h"

float  *NetZScoreMotorS(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter)
int       **distribution,*ndstr;
int       *size2;
int       edges;
node      *protein;
edge      *interaction;
int       *n_edgescr,*n_nodescr;
float     tolerance,threshold,thero;
int       iteration,diter;
float     maxmin[];
int       rnd[],linker[];
float    *nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
{
 int       size,random_iterations;
 int       **edistribution;
 int       *nedstr;
 int       **swap,**swap_edges;
 gradient  **d_nodescr,**d_edgescr;
 node      *dnode;
 edge      *dedge;
 float     *error;
 float     *dnodescr;
 float     rho,mind;
 int       i,j,jj,k,n,m,maxd,nstep,dsize,dedges;
 int       sizei2[2];
 float     *ModifyZScoreS();
 void      DualNodeSCR(),RandomSwapNodes();
 int       DualNode(),DualEdge(),mod();

 float     *fvector();
 int       *ivector();
 int       **imatrix();
 edge      *evector();
 node      *avector();
 float     *fvector();
 gradient  *gvector();
 gradient  **gmatrix();
 void      free_fvector();
 void      free_gvector();
 void      free_ivector();
 void      free_evector();
 void      free_avector();
 void      free_xvector();
 void      free_gmatrix();
 void      free_matrix();
 void      free_imatrix();

  size=size2[0];
  
  random_iterations=rnd[0];
  error=fvector(0,MAXERROR);
  swap=imatrix(0,size,0,random_iterations);
  d_nodescr=gmatrix(0,size,0,MAXWN);
  d_edgescr=gmatrix(0,edges,0,MAXWE);

  printf("\n*******************************\n");
  printf("RUNING NetZScore Static Random \n");
  printf("*******************************\n");
 
  for (i=0;i<size;i++){for (j=0;j<random_iterations;j++){swap[i][j]=0;}} 
  maxd=MAXID;
  RandomSwapNodes(size2,interaction,protein,distribution,ndstr,maxd,rnd,swap);
  printf("Swap random nodes. Done\n");
  if (maxmin[3]>0){
    edistribution=imatrix(0,MAXD,0,MAXIDE);
    nedstr=ivector(0,MAXD);
    for (i=0;i<MAXD;i++){for (j=0;j<MAXIDE;j++){ edistribution[i][j]=0;}}
    for (i=0;i<MAXD;i++){nedstr[i]=1;}
  }


  for (i=0;i<size;i++){nodescr[i]=0.0;average_nodescr[i]=0.0;rms_nodescr[i]=0.0;n_nodescr[i]=0;}
  for (i=0;i<edges;i++){edgescr[i]=0.0;average_edgescr[i]=0.0;rms_edgescr[i]=0.0;n_edgescr[i]=0;}

  rho=2*tolerance;
  n=0;
  mind=maxmin[5];
  nstep=0;

  if (maxmin[3]>0){
    printf("\nNodes on Dual space (Dimension %d )\n",edges);
    dnode=avector(0,edges);
    dsize=DualNode(dnode,interaction,protein,edges,mind,nodescr);
    printf("\nEdges on Dual space (Dimension %d )\n",dsize);
    dedge=evector(0,dsize);
    dedges=DualEdge(dedge,dnode,interaction,protein,dsize,edges,mind,nodescr);
    printf("\nIncrease Scoring on Dual-Dual Space (Dimension %d ) \n",dedges);
    dnodescr=fvector(0,dedges);
    DualNodeSCR(dnodescr,nodescr,dedges,dedge,interaction);
    printf("\nCalculate Edge Degree Distribution\n");
    EdgeDistribution(dnode,edges,edistribution,nedstr,maxmin);
    sizei2[0]=edges;
    sizei2[1]=dsize;
    printf("Swap random edges. Make space\n");
    swap_edges=imatrix(0,edges,0,random_iterations);
    printf("Swap random edges. Do\n");
    maxd=MAXIDE;
    RandomSwapNodes(sizei2,dedge,dnode,edistribution,nedstr,maxd,rnd,swap_edges);
    printf("Swap random edges. Done\n");
    nstep++;
    for(k=0;k<size;k++){
       if (nodescr[k]>0){
           average_nodescr[k]+=nodescr[k];
           rms_nodescr[k]+=nodescr[k]*nodescr[k];
           }
       }
    for(k=0;k<edges;k++){
       if (edgescr[k]>0){
           average_edgescr[k]+=edgescr[k];
           rms_edgescr[k]+=edgescr[k]*edgescr[k];
           }
       }
    }    


  while ( rho >= tolerance  && (n < iteration || iteration==0)){

   printf("\n\nIteration %d Starting Error: %f Tolerance: %e \n",n,rho,tolerance);

   if (mod(n,diter) && maxmin[3]>0 && n>0 ){
        free_evector(dedge,0,dsize);
        free_fvector(dnodescr,0,dedges);
        //printf("free swap_edges\n");
        free_imatrix(swap_edges,0,edges,0,random_iterations);
        printf("\nNodes on Dual space (Dimension %d Iteration %d)\n",edges,n);
        dsize=DualNode(dnode,interaction,protein,edges,mind,nodescr);
        printf("\nEdges on Dual space (Dimension %d Iteration %d)\n",dsize,n);
        dedge=evector(0,dsize);
        dedges=DualEdge(dedge,dnode,interaction,protein,dsize,edges,mind,nodescr);
        printf("\nIncrease Scoring on Dual-Dual Space (Dimension %d ) \n",dedges);
        dnodescr=fvector(0,dedges);
        DualNodeSCR(dnodescr,nodescr,dedges,dedge,interaction);
        printf("\nRestart Degree Distribution\n");
        for (i=0;i<MAXD;i++){for (j=0;j<MAXIDE;j++){ edistribution[i][j]=0;}}
        for (i=0;i<MAXD;i++){nedstr[i]=1;}
        printf("\nRecalculate Edge Degree Distribution\n");
        EdgeDistribution(dnode,edges,edistribution,nedstr,maxmin);
        sizei2[0]=edges;
        sizei2[1]=dsize;
        swap_edges=imatrix(0,edges,0,random_iterations);
        printf("\nRecalculate Swaping\n");
	maxd=MAXIDE;
        RandomSwapNodes(sizei2,dedge,dnode,edistribution,nedstr,maxd,rnd,swap_edges);
        printf("Swap random edges. Done\n");
        nstep++;
        for(k=0;k<size;k++){
         if (nodescr[k]>0){
           average_nodescr[k]+=nodescr[k];
           rms_nodescr[k]+=nodescr[k]*nodescr[k];
           }
         }
        for(k=0;k<edges;k++){
         if (edgescr[k]>0){
           average_edgescr[k]+=edgescr[k];
           rms_edgescr[k]+=edgescr[k]*edgescr[k];
           }
         }
        }

   if (iteration>0) error=ModifyZScoreS(n,diter,interaction,protein,dnode,dedge,edges,size,dedges,dsize,linker,edgescr,nodescr,dnodescr,d_edgescr,n_edgescr,d_nodescr,n_nodescr,swap,swap_edges,rnd,maxmin);

   rho=0;
   for (j=0;j<3;j++){if (error[j]>0 && rho<error[j]){rho=error[j];}}
   printf("\nError: %10.5e %10.5e %10.5e %10.5e Iteration: %d \n",rho,error[0],error[1],error[2],n);
   printf("NODES\n");
   for(k=0;k<size;k++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(nodescr[k]+protein[k].copy.score));
   }
   if (maxmin[3]>0 ){ 
   printf("EDGES\n");
   for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(edgescr[k]+interaction[k].association));
   }}
   nstep++;
   for(k=0;k<size;k++){
      if (nodescr[k]>0){
        average_nodescr[k]+=nodescr[k];
        rms_nodescr[k]+=nodescr[k]*nodescr[k];
       }
     }
   for(k=0;k<edges;k++){
      if (edgescr[k]>0){
        average_edgescr[k]+=edgescr[k];
        rms_edgescr[k]+=edgescr[k]*edgescr[k];
       }
     }
   n++;
  }
  printf("\nAVERAGE OF SCORES AND FLUCTUATION\n");
  printf("NODES\n");
  for(k=0;k<size;k++){
    average_nodescr[k]=average_nodescr[k]/nstep;
    if ((rms_nodescr[k]/nstep-average_nodescr[k]*average_nodescr[k])>1.0e-5){rms_nodescr[k]=sqrt(rms_nodescr[k]/nstep-average_nodescr[k]*average_nodescr[k]);}else{rms_nodescr[k]=0.0;}
    printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(average_nodescr[k]+protein[k].copy.score),rms_nodescr[k]);
  }
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     average_edgescr[k]=average_edgescr[k]/nstep;
     if ((rms_edgescr[k]/nstep-average_edgescr[k]*average_edgescr[k])>1.0e-5){rms_edgescr[k]=sqrt(rms_edgescr[k]/nstep-average_edgescr[k]*average_edgescr[k]);}else{rms_edgescr[k]=0.0;}
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(average_edgescr[k]+interaction[k].association),rms_edgescr[k]);
  }}



  if (maxmin[3]>0){
  for (i=0;i<MAXD;i++){
  for (j=0;j<MAXIDE;j++){
   printf("Edistr[%d][%d]=%d nedstr[%d]=%d\n",i,j,edistribution[i][j],i,nedstr[i]);
  }}
  printf("free edistribution\n");
    free_imatrix(edistribution,0,MAXD,0,MAXIDE);
  printf("free dnode\n");
    free_avector(dnode,0,edges);
  printf("free dedge\n");
    free_evector(dedge,0,dsize);
  printf("free dnodescr\n");
    free_fvector(dnodescr,0,dedges);
  printf("free swap_edges\n");
    free_imatrix(swap_edges,0,edges,0,random_iterations);
  printf("free nedstr\n");
    free_ivector(nedstr,0,MAXD);
  }
  printf("free d_edgescr\n");
  free_gmatrix(d_edgescr,0,edges,0,MAXWE);
  printf("free d_nodescr\n");
  free_gmatrix(d_nodescr,0,size,0,MAXWN);
  printf("free swap\n");
  free_imatrix(swap,0,size,0,random_iterations);
  printf("Done NetZscore\n");

  return error;

}

