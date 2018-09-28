#include "ppin.h"

float  *NetZScoreMotorD(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter)
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
 int       size;
 gradient  **d_nodescr,**d_edgescr;
 node      *dnode;
 edge      *dedge;
 float     *error;
 float     *dnodescr;
 float     rho,mind;
 int       i,j,jj,k,n,m,nstep,dsize,dedges;
 float     *ModifyZScoreD();
 void      DualNodeSCR(),EdgeDistribution();
 int       DualNode(),DualEdge(),mod();
 int       **edistribution;
 int       *nedstr;

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

  error=fvector(0,MAXERROR);
  d_nodescr=gmatrix(0,size,0,MAXWN);
  d_edgescr=gmatrix(0,edges,0,MAXWE);



  printf("\n*******************************\n");
  printf("RUNING NetZScore Dynamic \n");
  printf("*******************************\n");


  for (i=0;i<size;i++){nodescr[i]=0.0;average_nodescr[i]=0.0;rms_nodescr[i]=0.0;n_nodescr[i]=0;}
  for (i=0;i<edges;i++){edgescr[i]=0.0;average_edgescr[i]=0.0;rms_edgescr[i]=0.0;n_edgescr[i]=0;}

 
  rho=2*tolerance;
  n=0;
  mind=maxmin[5];
  nstep=0;




  while ( rho >= tolerance  && (n < iteration || iteration==0)){

   printf("\n\nIteration %d Starting Error: %f Tolerance: %e \n",n,rho,tolerance);


   if (iteration>0) error=ModifyZScoreD(n,diter,interaction,protein,dnode,dedge,edges,size,dedges,dsize,linker,edgescr,nodescr,dnodescr,d_edgescr,n_edgescr,d_nodescr,n_nodescr,distribution,ndstr,edistribution,nedstr,rnd,maxmin);

   rho=0;
   for (j=0;j<3;j++){if (error[j]>0 && rho<error[j]){rho=error[j];}}
   printf("\nError: %10.5e %10.5e %10.5e %10.5e Iteration: %d \n",rho,error[0],error[1],error[2],n);
   /*
   printf("NODES\n");
   for(k=0;k<size;k++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(nodescr[k]+protein[k].copy.score));
   }
   */
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
    if ((rms_nodescr[k]/nstep-average_nodescr[k]*average_nodescr[k])>0){rms_nodescr[k]=sqrt(rms_nodescr[k]/nstep-average_nodescr[k]*average_nodescr[k]);}else{rms_nodescr[k]=0.0;}
    if (protein[k].copy.score>0){printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(average_nodescr[k]+protein[k].copy.score),rms_nodescr[k]);}
  }




  free_gmatrix(d_edgescr,0,edges,0,MAXWE);
  free_gmatrix(d_nodescr,0,size,0,MAXWN);

  return error;

}

