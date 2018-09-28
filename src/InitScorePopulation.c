#include "ppin.h"
float  InitScorePopulation(d_scr,n_scr,dnode,dedge,size,edges,linker,linker_degree,maxmin,scr)
gradient **d_scr;
int       size,edges,linker,linker_degree;
int      *n_scr;
float    *scr, *maxmin;
node     *dnode;
edge     *dedge;
{

 float    zi,max,error,diff,min;
 int      i,degree_i;
 float    InitMatrix(),InitScore();

 max=maxmin[0];
 min=maxmin[1];

    for (i=0;i<size;i++){
       degree_i=dnode[i].degree;
       zi=InitScore(dnode,i,dedge,edges,linker,linker_degree);
       if (zi>max){zi=max;}
       if (zi<min){zi=min;}
       scr[i]=InitMatrix(i,maxmin,zi,dnode,dedge,size,linker,d_scr,n_scr);
       if (degree_i>0 && linker_degree>0){scr[i]=scr[i]/degree_i;}
       if (scr[i]>max){scr[i]=max;}
       diff+=scr[i];
      }
    error=diff*diff;


 return error;

}
