#include "ppin.h"

int  **Graph(interaction,edges,size)
edge  *interaction;
int   edges,size;
{
 int i,j,**g;
 int **imatrix();

 g=imatrix(0,size,0,size);
 for (i=0;i<size;i++){
 for (j=0;j<size;j++){
     g[i][j]=0;
 }}
 for (i=0;i<edges;i++){
   g[interaction[i].i][interaction[i].j]=1;
   g[interaction[i].j][interaction[i].i]=1;
 }
 for (i=0;i<size;i++){g[i][i]=0;}
 return g;
}
 
