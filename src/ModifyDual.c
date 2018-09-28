#include "ppin.h"

float  ModifyDual(dualedge,dnode,interaction,protein,dedges,edges,size,linker,nodescr,edgescr,maxmin)
int       dedges,edges,size,linker[];
float    *nodescr,*edgescr,maxmin[];
edge     *dualedge,*interaction;
node     *dnode,*protein;
{
  int i,j;
  float init,diff,error,errord,max,min,add;
  float funcE(),funcN();

  max=maxmin[0];
  min=maxmin[1];

  error=0;
  errord=0;

  if (linker[2]>0){
   for (i=0;i<size;i++){
    init= nodescr[i];
    add=0;   
    for (j=0;j<dedges;j++){
     if (!strcmp(dualedge[j].a.name1,protein[i].name1) || !strcmp(dualedge[j].a.name1,protein[i].name2)){add+=funcN(dualedge[j].association,nodescr[i],maxmin,linker)/2;}
     if (!strcmp(dualedge[j].a.name2,protein[i].name1) || !strcmp(dualedge[j].a.name2,protein[i].name2)){add+=funcN(dualedge[j].association,nodescr[i],maxmin,linker)/2;}
     if (!strcmp(dualedge[j].b.name1,protein[i].name1) || !strcmp(dualedge[j].b.name1,protein[i].name2)){add+=funcN(dualedge[j].association,nodescr[i],maxmin,linker)/2;}
     if (!strcmp(dualedge[j].b.name2,protein[i].name1) || !strcmp(dualedge[j].b.name2,protein[i].name2)){add+=funcN(dualedge[j].association,nodescr[i],maxmin,linker)/2;}
    }
 
    nodescr[i]+=add;
    if (nodescr[i]>max){nodescr[i]=max;} 

    diff=nodescr[i]-init;
    error+=diff*diff;
   }
  }

  if (linker[3]>0){
    for (i=0;i<edges;i++){
     init= edgescr[i];
     edgescr[i]+=funcE(dnode[i].copy.score,edgescr[i],maxmin,linker); 
     if (edgescr[i]>max){edgescr[i]=max;}
     diff=edgescr[i]-init;
     errord+=diff*diff;
    }
  }

/*
  for (i=0;i<edges;i++){
     interaction[i].a=protein[interaction[i].i];
     interaction[i].b=protein[interaction[i].j];
  }
*/

  if (size>0 && error>0){error=sqrt(error/size);}
  if (edges>0 && errord>0){errord=sqrt(errord/size);}
  if (error<errord){error=errord;}

  return error;
  

}

float funcN(x,z,maxmin,linker)
float x,z, maxmin[];
int   linker[];
{
 float y,max,min;

  max=maxmin[0];
  min=maxmin[1];

  y=0;
  if (linker[2]>0){if ((x+z)>max){y=max/linker[2];}else{y=(x+z)/linker[2];}}

 return y;
}

float funcE(x,z,maxmin,linker)
float x, maxmin[];
int   linker[];
{
 float y,max,min;

  max=maxmin[0];
  min=maxmin[1];

  y=0;
  if (linker[3]>0){if ((x+z)>max){y=max/linker[3];}else{y=(x+z)/linker[3];}}

 return y;
}


