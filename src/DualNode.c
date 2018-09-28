#include "ppin.h"
int   DualNode(dualnode,interaction,protein,edges,minim,nodescr)
int       edges;
node     *dualnode,*protein;
edge     *interaction;
float     minim,*nodescr;
{

  int       i,j,k;
  float     min();
  void      nrerror();
  float     association, add_i, add_j;

  for (i=0;i<edges;i++){
    strcpy(dualnode[i].name1,interaction[i].a.name1);  
    strcpy(dualnode[i].name2,interaction[i].b.name1);  
    dualnode[i].copy.cell=min(interaction[i].a.copy.cell,interaction[i].b.copy.cell);
    dualnode[i].copy.rms=min(interaction[i].a.copy.rms,interaction[i].b.copy.rms);
    dualnode[i].copy.score=interaction[i].association;
    //printf("Dual Node Score [%d] = %f\n",i,dualnode[i].copy.score);
    for (j=0;j<LOCAL;j++){
      dualnode[i].copy.local[j]=interaction[i].a.copy.local[j]*interaction[i].b.copy.local[j];
    }
  }

  k=0;
  for (i=0;i<edges;i++){
  for (j=i+1;j<edges;j++){
    if (   interaction[i].i==interaction[j].i || interaction[i].i==interaction[j].j
        || interaction[i].j==interaction[j].i || interaction[i].j==interaction[j].j){
        association=add_i=add_j=0;
        if (nodescr[interaction[i].i]>minim){add_i=nodescr[interaction[i].i];}else{add_i=0;}
        if (nodescr[interaction[i].j]>minim){add_j=nodescr[interaction[i].j];}else{add_j=0;}
        //if (add_i>0){printf("I %d J %d Add_i %f scr %f MIN %f I Node %d - Node %d & J Node %d - Node %d\n",i,j,add_i,nodescr[interaction[i].i],minim,interaction[i].i,interaction[i].j,interaction[j].j,interaction[j].i);}
        //if (add_j>0){printf("I %d J %d Add_j %f scr %f MIN %f I Node %d - Node %d & J Node %d - Node %d\n",i,j,add_j,nodescr[interaction[i].j],minim,interaction[i].j,interaction[i].i,interaction[j].i,interaction[j].j);}
        if (interaction[i].i==interaction[j].i){association+=protein[interaction[i].i].copy.score+add_i;}
        if (interaction[i].j==interaction[j].i){association+=protein[interaction[i].j].copy.score+add_j;}
        if (interaction[i].i==interaction[j].j){association+=protein[interaction[i].i].copy.score+add_i;}
        if (interaction[i].j==interaction[j].j){association+=protein[interaction[i].j].copy.score+add_j;}
        if (association > 0 ){ k++; }
    }
  }}    

  return k;

}

float min(a,b)
float a;
float b;
{ if (a<=b){return a;}else{return b;} }


