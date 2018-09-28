#include "ppin.h"
void DualNodeSCR(dnodescr,nodescr,dedges,dedge,interaction)
float  *dnodescr,*nodescr;
int    dedges;
edge   *interaction,*dedge;
{
 
 int k,j,i;
 float association,add_i,add_j;

 for (k=0;k<dedges;k++){
        association=add_i=add_j=0;
        i=dedge[k].i;
        j=dedge[k].j;
        add_i=nodescr[interaction[i].i];
        add_j=nodescr[interaction[i].j];
        if (interaction[i].i==interaction[j].i){association+=add_i;}
        if (interaction[i].j==interaction[j].i){association+=add_j;}
        if (interaction[i].i==interaction[j].j){association+=add_i;}
        if (interaction[i].j==interaction[j].j){association+=add_j;}
        dnodescr[k]=association;
  } 

}
