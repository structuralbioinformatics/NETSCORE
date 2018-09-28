#include "ppin.h"
void  ScorePopulation(d_scr,n_scr,dnode,dedge,size,maxmin,scr,linker,linker_degree)
gradient **d_scr;
int       size;
int      *n_scr, linker_degree,linker;
float    *scr, *maxmin;
node     *dnode;
edge     *dedge;
{

 float    zi,max,min;
 int      i,j,degree_i;
 float    Score();

  max=maxmin[0];
  min=maxmin[1];
  for (i=0;i<size;i++){
       degree_i=dnode[i].degree;
       zi=Score(d_scr,n_scr,i,maxmin);
       if (degree_i>0 && linker_degree>0){zi=zi/degree_i;}
       if (zi>max){zi=max;}
       if (zi<min){zi=0.0;}
       scr[i]=zi;
      }


}
