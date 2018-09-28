#include "ppin.h"
float *ModifyZScoreD(iter,diter,interaction,protein,dnode,dedge,edges,size,dedges,dsize,linker,edgescr,nodescr,dnodescr,d_edgescr,n_edgescr,d_nodescr,n_nodescr,distribution,ndstr,edistribution,nedstr,rnd,maxmin)
int       iter,diter,edges,size,dsize,dedges,linker[],rnd[];
int     **distribution,*ndstr,**edistribution,*nedstr;
float    *edgescr,*nodescr,*dnodescr;
int      *n_edgescr,*n_nodescr;
gradient **d_edgescr,**d_nodescr;
float     maxmin[];
node     *protein,*dnode;
edge     *interaction,*dedge;
{
 float *error;
 float *fvector();
 float  ModifyNodeZScoreD(),ModifyDual();
 float  frac,nfrac,max,min, maxe, mine, mind,minde,max_dummy, sigma, max2[MAXPARSCR];
 int    i,j,n;
 int    linker_node,linker_edge,linker_degree,indirect;

 error=fvector(0,5);

  max=maxmin[0];
  min=maxmin[1];
  sigma=maxmin[2];
  maxe=maxmin[3];
  mine=maxmin[4];
  mind=maxmin[5];
  minde=maxmin[6];

 linker_node=linker[0];
 linker_edge=linker[1];
 indirect=linker[4];
 linker_degree=linker[5];

 max_dummy=max+1.0;  // This is to avoid duplicating the transfer of NOde-Score to dual-edge
                     // Once is transferred by DualEdge it can not be transferred again in ModifyNode 
                     //
 printf("\nScore-Flow on Nodes (protein) NetZscore Iteration= %d\n",iter);

 max2[0]=max;
 max2[1]=min;
 max2[2]=sigma;
 max2[3]=minde;
 max2[4]=MAXID;
 error[0]=ModifyNodeZScoreD(iter,protein,interaction,size,edges,indirect,linker_node,linker_degree,nodescr,d_nodescr,n_nodescr,edgescr,distribution,ndstr,rnd,max2);


   for (i=0;i<edges;i++) edgescr[i]=0.0;
   error[1]=error[2]=0.0;
 
 return error;

}

