#include "ppin.h"
float  ZscoreS(iter,d_scr,n_scr,dnode,inode,dedge,size,edges,linker,linker_degree,swap,rnd,maxmin,gz)
int      iter,*n_scr,*rnd;
gradient **d_scr;
float     maxmin[];
int       size,edges,linker,linker_degree,inode;
int     **swap;
node     *dnode;
edge     *dedge;
gauss    *gz;
{
 int   i,j,k,jj,m,n,skip,round,iterations,random_size,rand_type,degree_i,degree_j;
 float max,min,sigma,y,z,*xr;
 int   jnode,max_degree,mind,*k_bind;
 int   *ivector();
 float *fvector();
 float Transform(),CumulativeMatrix(),RandomCumulativeMatrix();
 void  free_fvector(),free_ivector();
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
 y=0.0;
 n=0;
 z=0.0;
 g.average=0;
 g.rmsd=1;
 *gz=g;

 y=Transform(inode,maxmin,n_scr,d_scr);
 if (degree_i>0  && linker_degree >0)y=y/degree_i;
 if (y<=0){return z;}
 
 xr=fvector(0,iterations);

 while(n<iterations-1){
    jnode=swap[inode][n];
    degree_j=dnode[jnode].degree;
    k_bind=ivector(0,degree_j);
    for (j=0;j<degree_j;j++){
       if (dedge[dnode[jnode].interact[j]].i == jnode){k=dedge[dnode[jnode].interact[j]].j;}
       if (dedge[dnode[jnode].interact[j]].j == jnode){k=dedge[dnode[jnode].interact[j]].i;}
       k_bind[j]=k;
      }
    xr[n]=RandomCumulativeMatrix(iter,jnode,maxmin,dnode,dedge,linker,size,d_scr,n_scr,k_bind);
    if (degree_j>0 && linker_degree >0) xr[n]=xr[n]/degree_j;
    free_ivector(k_bind,0,degree_j);
    n++;
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


