#include "ppin.h"
#define eps 1.0e-6

void WGraph(g,dnode,interaction,scr,escr,edges,size)
float **g;
node  *dnode;
edge  *interaction;
float *scr,*escr;
int   edges,size;
{
 int i,j,ii,jj;
 float a,b,e;

 for (i=0;i<size;i++){
 for (j=0;j<size;j++){
     g[i][j]=0.0;
 }}
 for (i=0;i<edges;i++){
   ii=interaction[i].i;
   jj=interaction[i].j;
   a=dnode[ii].copy.score+scr[ii];
   b=dnode[jj].copy.score+scr[jj];
   e=interaction[i].association+escr[i];
   if (a<MINSHORT) a=MINSHORT;
   if (b<MINSHORT) b=MINSHORT;
   if (e > eps){
     g[ii][jj]=2/(e*(a+b));
     g[jj][ii]=2/(e*(a+b));
   }
 }
 for (i=0;i<size;i++){g[i][i]=0.0;}

}

#undef eps

