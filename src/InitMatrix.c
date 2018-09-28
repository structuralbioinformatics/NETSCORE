#include "ppin.h"
float InitMatrix(i,maxmin,z,dnode,dedge,size,linker,d_scr,n_scr)
int       i,size,linker,*n_scr;
float     z,maxmin[];
gradient **d_scr;
node     *dnode;
edge     *dedge;
{
 float e,mass,max,min;
 int   j,k,kk;
 float Transform();
 float LocalProduct();
 void  nrerror();

    max=maxmin[0]; 
    min=maxmin[1];
    mass=0;
    e=0;

    if (z>0){

    for (j=0;j<dnode[i].degree;j++){
      if (dedge[dnode[i].interact[j]].i == i){k=dedge[dnode[i].interact[j]].j;}
      if (dedge[dnode[i].interact[j]].j == i){k=dedge[dnode[i].interact[j]].i;}
      mass +=  dnode[k].copy.score*dedge[dnode[i].interact[j]].association ;
    }
    if (mass>0){
       for (j=0;j<dnode[i].degree;j++){
          if (dedge[dnode[i].interact[j]].i == i){k=dedge[dnode[i].interact[j]].j;}
          if (dedge[dnode[i].interact[j]].j == i){k=dedge[dnode[i].interact[j]].i;}
          if (n_scr[i] == MAXWN){printf("Warning Dimension Node/edge SCR[%d] (%d) > %d\n",i,n_scr[i],MAXWN);}
          if (n_scr[i] == MAXWE){printf("Error Dimension Node/edge SCR[%d] (%d) exceeded %d, change MAXWE\n",i,n_scr[i],MAXWE);nrerror("STOP\n");}
          if ((dnode[k].copy.score*dedge[dnode[i].interact[j]].association*LocalProduct(i,k,dnode) )>min){
            kk=n_scr[i];
            d_scr[i][kk].j=k;
            d_scr[i][kk].size = ( z/mass ) *  (dnode[k].copy.score*LocalProduct(i,k,dnode)*dedge[dnode[i].interact[j]].association ) ;
            if (d_scr[i][kk].size>max){d_scr[i][kk].size=max;}
            //printf("d_scr[%d]I[%d]J[%d]= %e N_SCR=%d\n",kk,i,d_scr[i][kk].j,d_scr[i][kk].size,n_scr[i]);
            n_scr[i]=n_scr[i]+1;
          }
       }
    } 
    e=Transform(i,maxmin,n_scr,d_scr);

    }

    return e;

}
