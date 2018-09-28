#include "ppin.h"

int **XFilterMatrix(dnode,dedge,size,threshold)
node    *dnode;
edge    *dedge;
int      size;
float    threshold;
{
 int   **c;
 int     n,k,m;
 int   **imatrix();

 c=imatrix(0,size,0,size);
 for (m=0;m<size;m++){for (n=0;n<size;n++){ c[m][n]=c[n][m]=0;}}
 for (m=0;m<size;m++){
   for (n=0;n<dnode[m].degree;n++){
       k=-1;
       if (dedge[dnode[m].interact[n]].i == m){k=dedge[dnode[m].interact[n]].j;}
       if (dedge[dnode[m].interact[n]].j == m){k=dedge[dnode[m].interact[n]].i;}
       if (k>=0 && dedge[dnode[m].interact[n]].association>threshold) c[m][k]=c[k][m]=1;
       }
   }

 return c;
}
