#include "ppin.h"
float  ModifyNodeZScoreD(iter,dnode,dedge,size,edges,indirect,linker,linker_degree,scr,d_scr,n_scr,escr,distribution,ndstr,rnd,maxmin)
int       iter,size,edges,linker,linker_degree,indirect,*rnd;
int     **distribution,*ndstr,*n_scr;
float    *scr,*escr,*maxmin;
gradient **d_scr;
node     *dnode;
edge     *dedge;
{
  int   i,j,k,kk,*n_scr_old,n_max;
  float dummy,diff,error,max,min,sigma,mind;
  gradient **d_old;
  gradient **gmatrix();
  float InitZscorePopulation(), InitMatrix(),CmpMatrix(),CumulativeMatrix();
  int   *ivector();
  float *fvector();
  void  free_gmatrix(),free_ivector(),free_fvector(), ZscorePopulation();

  max=maxmin[0];
  min=maxmin[1];
  sigma=maxmin[2];
  mind=maxmin[3];

  if (linker==0){linker=MINLINK;} 
  diff=0.0;
  if (iter==0){
     error=InitZscorePopulation(d_scr,n_scr,dnode,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,scr);
  }else{
    n_scr_old=ivector(0,size);
    n_max=0;
    for (i=0;i<size;i++){
      n_scr_old[i]=n_scr[i];
      if (n_scr[i]>n_max){n_max=n_scr[i];}
    }
    d_old=gmatrix(0,size,0,n_max);
    for (i=0;i<size;i++){
     kk=0;
     for (j=0;j<n_scr_old[i];j++){
         if (d_scr[i][j].j>size||d_scr[i][j].j<0){
            printf("Error d_scr[%d][%d].j = %d OUT OF RANGE\n");
            n_scr_old[i]=n_scr_old[i]-1;
         }else{
            d_old[i][kk].j=d_scr[i][kk].j;
            d_old[i][kk].size=d_scr[i][kk].size;
            kk++;
         }
     }
    }
    for (i=0;i<size;i++){
       dummy=CumulativeMatrix(iter,maxmin,indirect,linker,i,dnode,dedge,size,d_old,d_scr,n_scr,n_scr_old,mind,escr);
    }
    diff=CmpMatrix(d_old,d_scr,n_scr,size);
    error=diff;
    ZscorePopulation(iter,d_scr,n_scr,dnode,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,scr);
    free_gmatrix(d_old,0,size,0,n_max);
    free_ivector(n_scr_old,0,size);
  }

 

  return error;
 
}


