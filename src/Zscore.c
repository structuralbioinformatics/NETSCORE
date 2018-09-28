#include "ppin.h"

float  Zscore(iter,d_scr,n_scr,dnode,inode,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,gz)
int      iter,*n_scr;
gradient **d_scr;
float     maxmin[];
int       size,edges,linker,linker_degree,inode,*rnd;
int     **distribution,*ndstr;
node     *dnode;
edge     *dedge;
gauss    *gz;
{
 int   i,j,k,jj,m,n,skip,round,iterations,random_size,rand_type,degree_i,degree_k;
 int   k_rand,k_node,k_bind[MAXI];
 float max,min,sigma,x,y,z,*xr;
 float ran0(),ran3();
 int   idum,kdum,max_degree,mind;
 float *fvector();
 float Transform(),CumulativeMatrix(),RandomCumulativeMatrix();
 void  free_fvector();
 gauss g;
 gauss normal();

 max=maxmin[0];
 min=maxmin[1];
 sigma=maxmin[2];
 mind=maxmin[5];
 max_degree=(int)(maxmin[4]);
 iterations=rnd[0];
 rand_type=rnd[1];
 degree_i=dnode[inode].degree;
 x=0.0; 
 y=x;
 n=0;
 z=0.0;

 y=Transform(inode,maxmin,n_scr,d_scr);
 if (degree_i>0  && linker_degree>0)y=y/degree_i;
 if (y<=0){return z;}

 
 xr=fvector(0,iterations);

 if (rand_type == 1){
  m=1;
  while(n<iterations-1){
    xr[n]=x;
    for (j=0;j<degree_i;j++){k_bind[j]=0;}
    j=0;
    skip=0;
    round=0;
    while (skip==0 && j < degree_i ) {
       if (dedge[dnode[inode].interact[j]].i == inode){k=dedge[dnode[inode].interact[j]].j;}
       if (dedge[dnode[inode].interact[j]].j == inode){k=dedge[dnode[inode].interact[j]].i;}
       if (dnode[k].degree<MAXD-1){degree_k=dnode[k].degree;}else{degree_k=MAXD-1;}
       if (ndstr[degree_k]<max_degree){random_size=ndstr[degree_k];}else{random_size=max_degree-1;}
       idum   = (degree_i + 1)  * m ;
       kdum   = (int) ( 1000.0 * ran3(&idum) );
       k_node = (int) ( random_size * ran3(&kdum) );
       if (k_node>random_size){printf("Error K_NODE %d > %d\n",k_node,random_size);}
       k_rand = distribution[degree_k][k_node];
       if (k_rand==inode){
         if(k_node>0){k_rand=distribution[degree_k][k_node-1];}
         else        {k_rand=distribution[degree_k][k_node+1];}
       }
       for (jj=0;jj< j;jj++){ 
          if (k_rand == k_bind[jj]){skip=1; round++ ; break;}
       }
       if (round>0) { k_bind[j]=k_rand; j++; round=0;}
       if (skip==0 ) { k_bind[j]=k_rand; j++; round=0;}else{skip=0;}
        m++;
        if (m>1000){m=1;}
    }
    xr[n]=RandomCumulativeMatrix(iter,inode,maxmin,dnode,dedge,linker,size,d_scr,n_scr,k_bind);
    if (degree_i>0 && linker_degree >0) xr[n]=xr[n]/degree_i;
    n++;
  }
 }else{
  m=1;
  while(n<iterations-1){
    xr[n]=x;
    for (j=0;j<degree_i;j++){k_bind[j]=0;}
    for (j=0;j<degree_i;j++){
      idum   = (degree_i + 1) * m ;
      kdum   = (int) ( 1000.0 * ran3(&idum) );
      k_rand = (int) ( (size-1) * ran3(&kdum) );
      if (k_rand==inode){ k_rand=k_rand+1; }
      if (k_rand>=size){printf("Error K_NODE %d > %d\n",k_rand,size);k_rand=size-1;}
      k_bind[j]=k_rand;
      m++;
      if (m>1000){m=1;}
    }
    xr[n]=RandomCumulativeMatrix(iter,inode,maxmin,dnode,dedge,linker,size,d_scr,n_scr,k_bind);
    if (degree_i>0 && linker_degree >0) xr[n]=xr[n]/degree_i;
    n++;
  }
 }

 if (iterations>0){
  xr[iterations-1]=y;
  g = normal(xr,iterations);
  if (g.rmsd > 0) { z = ( y - g.average ) / g.rmsd ; }else{ if (y==g.average){z=0.0;g.rmsd=1.0;}else{z=y;g.average=0;g.rmsd=1;}}
  *gz = g;
 }else{
  g.average=0;
  g.rmsd=1;
  *gz=g;
  z=y;
 }


 free_fvector(xr,0,iterations);

 return z;
}


