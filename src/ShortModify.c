#include "ppin.h"
float *ShortModify(iter,diter,interaction,protein,dnode,dedge,edges,size,dedges,dsize,linker,edgescr,nodescr,dnodescr,maxmin)
int       iter,diter,edges,size,dsize,dedges,linker[];
float    *edgescr,*nodescr,*dnodescr;
float     maxmin[];
node     *protein,*dnode;
edge     *interaction,*dedge;
{
 float *error;
 float *fvector();
 float  ShortModifyNode(),ShortModifyDual();
 float  frac,nfrac,max,min, maxe, mine, mind,minde,max_dummy, sigma, max2[5];
 int    i,j,n;
 int    linker_node,linker_edge,indirect;
 clock_t   launch, done;
 float    diff;

 error=fvector(0,5);

 
  maxe=maxmin[3];
 

 linker_node=linker[0];
 linker_edge=linker[1];

 max_dummy=max+1.0;  // This is to avoid duplicating the transfer of NOde-Score to dual-edge
                     // Once is transferred by DualEdge it can not be transferred again in ModifyNode 
                     //
 printf("\nScore-Flow on Nodes (protein) NetShort Iteration=%d\n",iter);
 launch = clock();
 error[0]=ShortModifyNode(iter,protein,interaction,size,edges,linker_node,nodescr,edgescr,maxmin);
 done = clock();
 diff = 1000*(done - launch) / CLOCKS_PER_SEC;
 printf("Calculation of ShortModifyNode. Time %e ms\n",diff);

 if (maxe>0 ){

  frac=modff((float)(iter)/(float)(diter),&nfrac);
  if (iter==diter*(int)(nfrac)){n=0;}else{n=iter;}

 // n=iter;

  printf("\nScore-Flow on Dual Graph (edges)\n");
  launch = clock();
  error[1]=ShortModifyNode(n,dnode,dedge,edges,dedges,linker_edge,edgescr,dnodescr,maxmin);
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("Calculation of ShortModifyNode. Time %e ms\n",diff);

/*
  error[2]=ShortModifyDual(dedge,dnode,interaction,protein,dedges,edges,size,linker,nodescr,edgescr,max2);
*/
  error[2]=0.0;


 }else{
   error[1]=error[2]=0.0;
 }
 
 return error;

}

