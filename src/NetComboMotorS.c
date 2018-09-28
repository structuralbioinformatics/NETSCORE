#include "ppin.h"

float  *NetComboMotorS(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter)
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
float    *nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
{
 int        i,k,size,iteration_global, iteration_netscore, iteration_netzscore, iteration_netshort;
 int        n_size,a_size,n_edges,a_edges;
 int       *n_edgescr_dummy,*n_nodescr_dummy; 
 float     *error, *error_dummy;
 float     mean,rms,mean_avg,rms_avg;
 float     mean_edge,rms_edge,mean_avg_edge,rms_avg_edge;
 float     *nodescr_dummy,*edgescr_dummy,*average_nodescr_dummy,*average_edgescr_dummy,*rms_nodescr_dummy,*rms_edgescr_dummy;
 float     *NetScoreMotor(),*NetShortMotor(),*NetZScoreMotorS();
 float     *fvector();
 int       *ivector();
 void      free_fvector(),free_ivector();
 clock_t   launch, done;
 double     diff;

 size=size2[0];
 error=fvector(0,MAXERROR);
 error_dummy=fvector(0,MAXERROR);
 nodescr_dummy=fvector(0,size);
 edgescr_dummy=fvector(0,edges);
 average_nodescr_dummy=fvector(0,size);
 average_edgescr_dummy=fvector(0,edges);
 rms_nodescr_dummy=fvector(0,size);
 rms_edgescr_dummy=fvector(0,edges);
 n_nodescr_dummy=ivector(0,size);
 n_edgescr_dummy=ivector(0,edges);

 iteration_netscore=iteration[1];
 iteration_netzscore=iteration[2];
 iteration_netshort=iteration[3];

  if (iteration_netscore<=0 && iteration_netzscore<=0 && iteration_netshort<=0){
    if (iteration[0]>0){
     iteration_netscore=iteration[0];
     iteration_netzscore=iteration[0];
     iteration_netshort=iteration[0];
    }else{
     iteration_netscore=1;
     iteration_netzscore=1;
     iteration_netshort=1;
    }
  }


  printf("\n*******************************\n");
  printf("RUNING NetCombo Static Random \n");
  printf("*******************************\n");


 for (i=0;i<size;i++){nodescr[i]=0.0;average_nodescr[i]=0.0;rms_nodescr[i]=0.0;n_nodescr[i]=0;}
 for (i=0;i<edges;i++){edgescr[i]=0.0;average_edgescr[i]=0.0;rms_edgescr[i]=0.0;n_edgescr[i]=0;}
 for(i=0;i<MAXERROR;i++){error[i]=error_dummy[i]=0.0;}

// ******************************************************************************** 
// Calculate NetScore

 printf("NetScore (iteration %d)\n",iteration_netscore); 
 if (iteration_netscore>0){

 launch = clock();
 error_dummy= NetScoreMotor(edges,size2,interaction,protein,nodescr_dummy,average_nodescr_dummy,rms_nodescr_dummy,n_nodescr_dummy,edgescr_dummy,average_edgescr_dummy,rms_edgescr_dummy,n_edgescr_dummy,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netscore,diter);
 for(i=0;i<MAXERROR;i++)if (error_dummy[i]>0 && error[i]<error_dummy[i]>0)error[i]=error_dummy[i];

// Transform into Zscores the increasing scores NODESCR
 mean=rms=mean_avg=rms_avg=mean_edge=rms_edge=mean_avg_edge=rms_avg_edge=0;

 n_size=a_size=0;
 for (i=0;i<size;i++){
   mean+=nodescr_dummy[i];
   rms+=nodescr_dummy[i]*nodescr_dummy[i];
   if (nodescr_dummy[i]!=0 && rnd[5]==1){ n_size++;}else{n_size++;}
   mean_avg+=average_nodescr_dummy[i];
   if (average_nodescr_dummy[i]!=0 && rnd[5]==1) {a_size++;}else{a_size++;}
   rms_avg+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   rms_nodescr[i]+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   }
 if (n_size>0){
   mean=mean/n_size;
   rms=rms/n_size-mean*mean;
   if (rms>0) {rms=sqrt(rms);}else{rms=1.0;}
   }
 if (a_size>0){ mean_avg=mean_avg/a_size; rms_avg=sqrt(rms_avg/a_size); }else{mean_avg=0.0;rms_avg=1.0;}
 if (rms==0.0)rms=1.0;
 if (rms_avg==0.0)rms_avg=1.0;
 
 for (i=0;i<size;i++){
    nodescr[i] += (nodescr_dummy[i] - mean)/rms;
    average_nodescr[i] += (average_nodescr_dummy[i] - mean_avg)/ rms_avg;
   }
 
 
// Transform into Zscores the increasing scores EDGESCR
 if (maxmin[3]>0 ){
   n_edges=a_edges=0;
   for (i=0;i<edges;i++){
     mean_edge+=edgescr_dummy[i];
     rms_edge+=edgescr_dummy[i]*edgescr_dummy[i];
     if (edgescr_dummy[i]!=0) n_edges++;
     mean_avg_edge+=average_edgescr_dummy[i];
     if (average_edgescr_dummy[i]!=0) a_edges++;
     rms_avg_edge+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     rms_edgescr[i]+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     }
   if (n_edges>0){
     mean_edge=mean_edge/n_edges;
     rms_edge=rms_edge/n_edges-mean_edge*mean_edge;
     if (rms_edge>0) {rms_edge=sqrt(rms_edge);}else{rms_edge=1.0;}
     }
   if (a_edges>0){ mean_avg_edge=mean_avg_edge/a_edges; rms_avg_edge=sqrt(rms_avg_edge/a_edges);}else{mean_avg_edge=0.0;rms_avg_edge=1.0;}
   if (rms_edge==0.0)rms_edge=1.0;
   if (rms_avg_edge==0.0)rms_avg_edge=1.0;

   for (i=0;i<edges;i++){
     edgescr[i] += (edgescr_dummy[i] - mean_edge)/rms_edge ;
     average_edgescr[i] += (average_edgescr_dummy[i] - mean_avg_edge)/ rms_avg_edge;
     }
   }else{
     for (i=0;i<edges;i++){
     edgescr[i] =0.0;
     average_edgescr[i] =0.0;
     rms_edgescr[i] =0.0;
     }
     rms_edge=1.0;
     rms_avg_edge=1.0;
   }

// Select maximum n_nodescr
   for (i=0;i<size;i++) if (n_nodescr_dummy[i]>n_nodescr[i]) n_nodescr[i]=n_nodescr_dummy[i];
   if (maxmin[3]>0 ) for (i=0;i<edges;i++) if (n_edgescr_dummy[i]>n_edgescr[i]) n_edgescr[i]=n_edgescr_dummy[i];

// Print Partial Modification by NetScore
  printf("\nNETSCORE Z-SCORE BY NETSCORE\n");   
  /*
  printf("NODES\n");
  for(k=0;k<size;k++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((nodescr_dummy[k] - mean)/rms +protein[k].copy.score));
   }
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (edgescr_dummy[k] - mean_edge)/rms_edge +interaction[k].association));
   }}
   */
  printf("\nAVERAGE OF SCORES AND FLUCTUATION\n");
  printf("NODES\n");
  for(k=0;k<size;k++){
    if (protein[k].copy.score>0){
    printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((average_nodescr_dummy[k] - mean_avg)/ rms_avg +protein[k].copy.score),rms_nodescr_dummy[k]);
  }}
  /*
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (average_edgescr_dummy[k] - mean_avg_edge)/ rms_avg_edge +interaction[k].association),rms_edgescr_dummy[k]);
  }}
  */
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("End calculation by NetScore. Time %e ms\n",diff);

  } // End Calculation with NetScore

// ******************************************************************************** 
// Calculate NetZScore (static) 

 printf("NetZScore Static (iteration %d)\n",iteration_netzscore); 
 if (iteration_netzscore>0){
 launch = clock();
 error_dummy= NetZScoreMotorS(edges,size2,interaction,protein,nodescr_dummy,average_nodescr_dummy,rms_nodescr_dummy,n_nodescr_dummy,edgescr_dummy,average_edgescr_dummy,rms_edgescr_dummy,n_edgescr_dummy,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netzscore,diter);


  for(i=0;i<MAXERROR;i++)if (error_dummy[i]>0 && error[i]<error_dummy[i]>0)error[i]=error_dummy[i];
 

// Transform into Zscores the increasing scores NODESCR
 mean=rms=mean_avg=rms_avg=mean_edge=rms_edge=mean_avg_edge=rms_avg_edge=0;
 n_size=a_size=0;
 for (i=0;i<size;i++){
   mean+=nodescr_dummy[i];
   rms+=nodescr_dummy[i]*nodescr_dummy[i];
   if (nodescr_dummy[i]!=0 && rnd[5]==1){ n_size++;}else{n_size++;}
   mean_avg+=average_nodescr_dummy[i];
   if (average_nodescr_dummy[i]!=0 && rnd[5]==1) {a_size++;}else{a_size++;}
   rms_avg+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   rms_nodescr[i]+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   }
 if (n_size>0){
   mean=mean/n_size;
   rms=rms/n_size-mean*mean;
   if (rms>0) {rms=sqrt(rms);}else{rms=1.0;}
   }
 if (a_size>0){ mean_avg=mean_avg/a_size; rms_avg=sqrt(rms_avg/a_size); }else{mean_avg=0.0;rms_avg=1.0;}
 if (rms==0.0)rms=1.0;
 if (rms_avg==0.0)rms_avg=1.0;
 
 for (i=0;i<size;i++){
    nodescr[i] += (nodescr_dummy[i] - mean)/rms;
    average_nodescr[i]  += (average_nodescr_dummy[i] - mean_avg)/ rms_avg;
   }
 
 
 printf("Transform into Zscores the increasing scores EDGESCR\n");
// Transform into Zscores the increasing scores EDGESCR
 if (maxmin[3]>0 ){
   n_edges=a_edges=0;
   for (i=0;i<edges;i++){
     mean_edge+=edgescr_dummy[i];
     rms_edge+=edgescr_dummy[i]*edgescr_dummy[i];
     if (edgescr_dummy[i]!=0) n_edges++;
     mean_avg_edge+=average_edgescr_dummy[i];
     if (average_edgescr_dummy[i]!=0) a_edges++;
     rms_avg_edge+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     rms_edgescr[i]+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     }
   if (n_edges>0){
     mean_edge=mean_edge/n_edges;
     rms_edge=rms_edge/n_edges-mean_edge*mean_edge;
     if (rms_edge>0) {rms_edge=sqrt(rms_edge);}else{rms_edge=1.0;}
     }
   if (a_edges>0){ mean_avg_edge=mean_avg_edge/a_edges; rms_avg_edge=sqrt(rms_avg_edge/a_edges);}else{mean_avg_edge=0.0;rms_avg_edge=1.0;}
   if (rms_edge==0.0)rms_edge=1.0;
   if (rms_avg_edge==0.0)rms_avg_edge=1.0;

   for (i=0;i<edges;i++){
     edgescr[i] += (edgescr_dummy[i] - mean_edge)/rms_edge ;
     average_edgescr[i] += (average_edgescr_dummy[i] - mean_avg_edge)/ rms_avg_edge;
     }
   }else{
     for (i=0;i<edges;i++){
     edgescr[i] =0.0;
     average_edgescr[i] =0.0;
     rms_edgescr[i] =0.0;
     }
     rms_edge=1.0;
     rms_avg_edge=1.0;
   }

// Select maximum n_nodescr
   for (i=0;i<size;i++) if (n_nodescr_dummy[i]>n_nodescr[i]) n_nodescr[i]=n_nodescr_dummy[i];
   if (maxmin[3]>0 ) for (i=0;i<edges;i++) if (n_edgescr_dummy[i]>n_edgescr[i]) n_edgescr[i]=n_edgescr_dummy[i];

// Print Partial Modification by NetZScore
  printf("\nNETSCORE Z-SCORE BY NETZSCORE\n");   
  /*
  printf("NODES\n");
  for(k=0;k<size;k++){
    if (protein[k].copy.score>0){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((nodescr_dummy[k] - mean)/rms +protein[k].copy.score));
   }}
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (edgescr_dummy[k] - mean_edge)/rms_edge +interaction[k].association));
   }}
   */
  printf("\nAVERAGE OF SCORES AND FLUCTUATION\n");
  printf("NODES\n");
  for(k=0;k<size;k++){
    if (protein[k].copy.score>0){
    printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((average_nodescr_dummy[k] - mean_avg)/ rms_avg +protein[k].copy.score),rms_nodescr_dummy[k]);
  }}
  /*
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (average_edgescr_dummy[k] - mean_avg_edge)/ rms_avg_edge +interaction[k].association),rms_edgescr_dummy[k]);
  }}
  */

  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("End calculation by NetZScore. Time %e ms\n",diff);
  } // End Calculation with NetZScore


// ******************************************************************************** 
// Calculate NetShort 

 printf("Netshort (iteration %d)\n",iteration_netshort);
 if (iteration_netshort>0){

 
 launch = clock();
 error_dummy= NetShortMotor(edges,size2,interaction,protein,nodescr_dummy,average_nodescr_dummy,rms_nodescr_dummy,n_nodescr_dummy,edgescr_dummy,average_edgescr_dummy,rms_edgescr_dummy,n_edgescr_dummy,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netshort,diter);
 for(i=0;i<MAXERROR;i++)if (error_dummy[i]>0 && error[i]<error_dummy[i]>0)error[i]=error_dummy[i];



// Transform into Zscores the increasing scores NODESCR
 mean=rms=mean_avg=rms_avg=mean_edge=rms_edge=mean_avg_edge=rms_avg_edge=0;
 n_size=a_size=0;
 for (i=0;i<size;i++){
   mean+=nodescr_dummy[i];
   rms+=nodescr_dummy[i]*nodescr_dummy[i];
   if (nodescr_dummy[i]!=0 && rnd[5]==1){ n_size++;}else{n_size++;}
   mean_avg+=average_nodescr_dummy[i];
   if (average_nodescr_dummy[i]!=0 && rnd[5]==1) {a_size++;}else{a_size++;}
   rms_avg+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   rms_nodescr[i]+=rms_nodescr_dummy[i]*rms_nodescr_dummy[i];
   }
 if (n_size>0){
   mean=mean/n_size;
   rms=rms/n_size-mean*mean;
   if (rms>0) {rms=sqrt(rms);}else{rms=1.0;}
   }
 if (a_size>0){ mean_avg=mean_avg/a_size; rms_avg=sqrt(rms_avg/a_size); }else{mean_avg=0.0;rms_avg=1.0;}
 if (rms==0.0)rms=1.0;
 if (rms_avg==0.0)rms_avg=1.0;
 
 for (i=0;i<size;i++){
    nodescr[i] += (nodescr_dummy[i] - mean)/rms;
    average_nodescr[i]  += (average_nodescr_dummy[i] - mean_avg)/ rms_avg;
   }
 
 
// Transform into Zscores the increasing scores EDGESCR
 if (maxmin[3]>0 ){
   n_edges=a_edges=0;
   for (i=0;i<edges;i++){
     mean_edge+=edgescr_dummy[i];
     rms_edge+=edgescr_dummy[i]*edgescr_dummy[i];
     if (edgescr_dummy[i]!=0) n_edges++;
     mean_avg_edge+=average_edgescr_dummy[i];
     if (average_edgescr_dummy[i]!=0) a_edges++;
     rms_avg_edge+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     rms_edgescr[i]+=rms_edgescr_dummy[i]*rms_edgescr_dummy[i];
     }
   if (n_edges>0){
     mean_edge=mean_edge/n_edges;
     rms_edge=rms_edge/n_edges-mean_edge*mean_edge;
     if (rms_edge>0) {rms_edge=sqrt(rms_edge);}else{rms_edge=1.0;}
     }
   if (a_edges>0){ mean_avg_edge=mean_avg_edge/a_edges; rms_avg_edge=sqrt(rms_avg_edge/a_edges);}else{mean_avg_edge=0.0;rms_avg_edge=1.0;}
   if (rms_edge==0.0)rms_edge=1.0;
   if (rms_avg_edge==0.0)rms_avg_edge=1.0;

   for (i=0;i<edges;i++){
     edgescr[i] += (edgescr_dummy[i] - mean_edge)/rms_edge ;
     average_edgescr[i] += (average_edgescr_dummy[i] - mean_avg_edge)/ rms_avg_edge;
     }
   }else{
     for (i=0;i<edges;i++){
     edgescr[i] =0.0;
     average_edgescr[i] =0.0;
     rms_edgescr[i] =0.0;
     }
     rms_edge=1.0;
     rms_avg_edge=1.0;
   }

 
// Select maximum n_nodescr
   for (i=0;i<size;i++) if (n_nodescr_dummy[i]>n_nodescr[i]) n_nodescr[i]=n_nodescr_dummy[i];
   if (maxmin[3]>0 ) for (i=0;i<edges;i++) if (n_edgescr_dummy[i]>n_edgescr[i]) n_edgescr[i]=n_edgescr_dummy[i];

// Print Partial Modification by NetSshor
  printf("\nNETSCORE Z-SCORE BY NETSHORT\n");   
  /*
  printf("NODES\n");
  for(k=0;k<size;k++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((nodescr_dummy[k] - mean)/rms +protein[k].copy.score));
   }
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (edgescr_dummy[k] - mean_edge)/rms_edge +interaction[k].association));
   }}
   */
  printf("\nAVERAGE OF SCORES AND FLUCTUATION\n");
  printf("NODES\n");
  for(k=0;k<size;k++){
  if (protein[k].copy.score>0){
    printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,((average_nodescr_dummy[k] - mean_avg)/ rms_avg +protein[k].copy.score),rms_nodescr_dummy[k]);
  }}
  /*
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,( (average_edgescr_dummy[k] - mean_avg_edge)/ rms_avg_edge +interaction[k].association),rms_edgescr_dummy[k]);
  }}
  */
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("End calculation by NetShort. Time %e ms\n",diff);

  } // End Calculation with NetShort


// ******************************************************************************** 
// Final average on scores

  for (i=0;i<size;i++){
    nodescr[i] = nodescr[i]/3.0;
    average_nodescr[i]  = average_nodescr[i]/3.0;
    rms_nodescr[i] = sqrt(rms_nodescr[i]/3.0);
   }
  if (maxmin[3]>0 ){
    for (i=0;i<edges;i++){
      edgescr[i] = edgescr[i]/3.0;
      average_edgescr[i] = average_edgescr[i]/3.0;
      rms_edgescr[i] = sqrt(rms_edgescr[i]/3.0);
      }
   }else{
    for (i=0;i<edges;i++){
      edgescr[i] = 0;
      average_edgescr[i] = 0;
      rms_edgescr[i] = 0;
      }
   }


  printf("\nNETCOMBO SCORE \n");   
  /*
  printf("NODES\n");
  for(k=0;k<size;k++){
     printf("%10d\t%10s [%10.5f] => [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(nodescr[k]+protein[k].copy.score));
   }
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(edgescr[k]+interaction[k].association));
   }}
   */
  printf("\nAVERAGE OF SCORES AND FLUCTUATION\n");
  printf("NODES\n");
  for(k=0;k<size;k++){
  if (protein[k].copy.score>0){
    printf("%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(average_nodescr[k]+protein[k].copy.score),rms_nodescr[k]);
  }}
  /*
  if (maxmin[3]>0 ){
  printf("EDGES\n");
  for(k=0;k<edges;k++){
     printf("%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(average_edgescr[k]+interaction[k].association),rms_edgescr[k]);
  }}
  */
 printf("Done Calculation with NetCombo\n");

// Free memory

 free_fvector(error_dummy,0,MAXERROR);
 free_fvector(nodescr_dummy,0,size);
 free_fvector(edgescr_dummy,0,edges);
 free_fvector(average_nodescr_dummy,0,size);
 free_fvector(average_edgescr_dummy,0,edges);
 free_fvector(rms_nodescr_dummy,0,size);
 free_fvector(rms_edgescr_dummy,0,edges);
 free_ivector(n_nodescr_dummy,0,size);
 free_ivector(n_edgescr_dummy,0,edges);
 

 return error;

}
