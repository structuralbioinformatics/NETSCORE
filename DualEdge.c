#include "ppin.h"
int     DualEdge(dualedge,dualnode,interaction,protein,size,edges,min,nodescr)
edge   *interaction,*dualedge;
node   *protein,*dualnode;
int     size,edges;
float   min,*nodescr;
{

  int         i,j,k,ii,skip;
  void        nrerror();
  float       association, add_i, add_j;

  ii=0;
  for (i=0;i<edges;i++){
    for (j=i+1;j<edges;j++){
     if (   interaction[i].i==interaction[j].i || interaction[i].i==interaction[j].j
         || interaction[i].j==interaction[j].i || interaction[i].j==interaction[j].j){
       skip=0;
/*
       for (k=0;k<ii;k++){
         if ( (dualedge[k].i==i && dualedge[k].j==j) || (dualedge[k].j==i && dualedge[k].i==j) ){skip++;printf("Found redundancy Dedge[%d]I=%d Dedge[%d]J=%d and  Edges [%d][%d]\n",k,dualedge[k].i,k,dualedge[k].j,i,j);break;}
       }
*/
       if (skip==0){
        dualedge[ii].i=i;
        dualedge[ii].j=j;
        dualedge[ii].a=dualnode[dualedge[ii].i];
        dualedge[ii].b=dualnode[dualedge[ii].j];
        association=add_i=add_j=0;
        if (nodescr[interaction[i].i]>min){add_i=nodescr[interaction[i].i];}else{add_i=0;}
        if (nodescr[interaction[i].j]>min){add_j=nodescr[interaction[i].j];}else{add_j=0;}
        if (interaction[i].i==interaction[j].i){association+=protein[interaction[i].i].copy.score+add_i;}
        if (interaction[i].j==interaction[j].i){association+=protein[interaction[i].j].copy.score+add_j;}
        if (interaction[i].i==interaction[j].j){association+=protein[interaction[i].i].copy.score+add_i;}
        if (interaction[i].j==interaction[j].j){association+=protein[interaction[i].j].copy.score+add_j;}
        if (association > 0 ){ 
          dualedge[ii].association=association;
        //printf("DualEdge.association[%d]= %f from I[%d] - J[%d] in nodes N[%d]= %s and N[%d]= %s with nodes N[%d]= %s and N[%d]= %s\n",ii,dualedge[ii].association,i,j,interaction[i].i,protein[interaction[i].i].name1, interaction[i].j, protein[interaction[i].j].name1,interaction[j].i,protein[interaction[j].i].name1, interaction[j].j, protein[interaction[j].j].name1);
         ii++;
        }
        if (ii>size){nrerror("Unexpected number of edges in DualEdge");}
       }
     }
    }
  }

  for (i=0;i<edges;i++){
  dualnode[i].degree=0;
  //printf("DUALNODE %s-%s %d \n",dualnode[i].name1,dualnode[i].name2,i);
  for (j=0;j<ii;j++){
    if (  dualedge[j].i == i || dualedge[j].j == i){
       dualnode[i].interact[dualnode[i].degree]=j;
       //printf("DUAL INTERACTION %s-%s %d and %s-%s %d %s-%s %d EDGE %d Degree[%d]=%d\n",dualnode[i].name1,dualnode[i].name2,i,dualnode[dualedge[j].i].name1,dualnode[dualedge[j].i].name2,dualedge[j].i,dualnode[dualedge[j].j].name1,dualnode[dualedge[j].j].name2,dualedge[j].j,j,i,dualnode[i].degree);
       dualnode[i].degree++;
       if (dualnode[i].degree>MAXI){nrerror("Degree larger than MAXI");}
    }
  }}    


  return ii;
 
}     
  

