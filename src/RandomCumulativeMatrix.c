#include "ppin.h"
float RandomCumulativeMatrix(iter,i,maxmin,dnode,dedge,linker,size,d_scr,n_scr,k_random)
int       iter,i,size,linker;
float     maxmin[];
int      *n_scr;
gradient **d_scr;
int       k_random[];
node     *dnode;
edge     *dedge;
{
 float e,add,max,min;
 int j,k,kk;


  max=maxmin[0];
  min=maxmin[1];
  e=0.0;

  for (j=0;j<dnode[i].degree;j++){
    k=k_random[j];
    add=0;
    for (kk=0;kk<n_scr[k];kk++){
      if (d_scr[k][kk].j!=i){
        add += d_scr[k][kk].size;
      }
    }
    if (add>max){add=max;}
    if (add>min){e+=add;}
  }
  //e=e/pow(linker,(iter));
  e=e/(exp(iter*log((float)linker)));
  if (e>max){e=max;}
  if (e<min){e=0;}

 return e;
}
