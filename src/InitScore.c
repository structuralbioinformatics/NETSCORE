#include "ppin.h"
float  InitScore(dnode,inode,dedge,edges,linker,linker_degree)
int       edges,linker,linker_degree,inode;
node     *dnode;
edge     *dedge;
{
 int   j,k,degree_i;
 float y,z;
 float LocalProduct();
 gauss g;

 degree_i=dnode[inode].degree;
 y=z=0.0;

 for (j=0;j<degree_i;j++){
   if (dedge[dnode[inode].interact[j]].i == inode){k=dedge[dnode[inode].interact[j]].j;}
   if (dedge[dnode[inode].interact[j]].j == inode){k=dedge[dnode[inode].interact[j]].i;}
   y+=dnode[k].copy.score*LocalProduct(inode,k,dnode)*dedge[dnode[inode].interact[j]].association/linker ;
 }

 if (degree_i>0 && linker_degree >0) {z=y/degree_i;}else{z=y;}


 return z;
}


