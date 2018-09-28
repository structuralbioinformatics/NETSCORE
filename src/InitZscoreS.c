#include "ppin.h"
float  InitZscoreS(dnode,inode,dedge,size,edges,linker,linker_degree,swap,rnd,gz)
int       size,edges,linker,linker_degree,inode,*rnd;
int     **swap;
node     *dnode;
edge     *dedge;
gauss    *gz;
{
 int   i,j,k,jj,n,skip,round,iterations,random_size,rand_type,degree_i,degree_j;
 float y,z,r,*xr;
 int   jnode;
 float LocalProduct();
 float *fvector();
 void  free_fvector();
 gauss g;
 gauss normal();

 iterations=rnd[0];
 rand_type=rnd[1];
 y=0.0;
 n=0;
 
 degree_i=dnode[inode].degree;

 if (degree_i <= 0){
  g.average=0;
  g.rmsd=1;
  *gz=g;
  z=y;
  return z;
 }


 xr=fvector(0,iterations);

 for (j=0;j<degree_i;j++){
   if (dedge[dnode[inode].interact[j]].i == inode){k=dedge[dnode[inode].interact[j]].j;}
   if (dedge[dnode[inode].interact[j]].j == inode){k=dedge[dnode[inode].interact[j]].i;}
   y+=dnode[k].copy.score*LocalProduct(inode,k,dnode)*dedge[dnode[inode].interact[j]].association/linker ;
 }
 if (degree_i>0 && linker_degree>0) y=y/degree_i;

 //printf("InitZscore[%d]= %f\n",inode,y);

 while(n<iterations-1){
    xr[n]=0.0;
    jnode=swap[inode][n];
    degree_j=dnode[jnode].degree;
    for (j=0;j<degree_j;j++){
      if (dedge[dnode[jnode].interact[j]].i == jnode){k=dedge[dnode[jnode].interact[j]].j;}
      if (dedge[dnode[jnode].interact[j]].j == jnode){k=dedge[dnode[jnode].interact[j]].i;}
       //printf("Random InitZscoreS[%d] XR[%d]= %f add[%d] Score[%d] %f x Local %f x Edge %f / linker %d\n",jnode,n,xr[n],j,k,dnode[k].copy.score,LocalProduct(jnode,k,dnode),dedge[dnode[jnode].interact[j]].association,linker);
       xr[n]+=dnode[k].copy.score*LocalProduct(jnode,k,dnode)*dedge[dnode[jnode].interact[j]].association/linker ;
     }
    if (degree_j>0 && linker_degree>0) xr[n]=xr[n]/degree_j;
    n++;
  }

 if (iterations>0){
  xr[iterations-1]=y;
  g = normal(xr,iterations);
  if (g.rmsd > 0) { z = ( y - g.average ) / g.rmsd ; }else{ if (y==g.average){z=0.0;g.rmsd=1.0;}else{z=y;g.average=0;g.rmsd=1.0;}}
  *gz = g;
  //printf("Modified InitZscore[%d]= %f\n",inode,z);
 }else{
  g.average=0;
  g.rmsd=1;
  *gz=g;
  z=y;
  //printf("Modified InitZscore[%d]= %f\n",inode,z);
 } 


 free_fvector(xr,0,iterations);

 return z;
}


