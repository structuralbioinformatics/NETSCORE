#include "ppin.h"
void  ZscorePopulation(iter,d_scr,n_scr,dnode,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,scr)
gradient **d_scr;
int       iter,size,edges,linker_degree,linker,*rnd;
int     **distribution,*ndstr,*n_scr;
float    *scr, *maxmin;
node     *dnode;
edge     *dedge;
{

 float    *z,zi,rmsdi,s,average,rmsd,max,degree_i,sigma,min;
 gauss    gz;
 int      i,j;
 float    *fvector(), Zscore();
 void     free_fvector();


  z=fvector(0,size);
  max=maxmin[0];
  sigma=maxmin[2];
  min=maxmin[1];

  switch (rnd[2]){
  case 0:
    average=rmsdi=rmsd=0;
    for (i=0;i<size;i++){
       zi=Zscore(iter,d_scr,n_scr,dnode,i,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,&gz);
       z[i]=zi*gz.rmsd + gz.average ;
       if (gz.average==0){
          average+=z[i];
          rmsdi+=z[i]*z[i];
         }else{
          average+= gz.average ;
          rmsdi+=gz.rmsd*gz.rmsd + gz.average*gz.average;
         }
      }
    if (size>0)average=average/size;
    if (size>0)rmsdi=rmsdi/size;
    s=rmsdi-average*average;
    if (s>0) rmsd=sqrt(s);
    for (i=0;i<size;i++){
       if (rmsd>0){zi=(z[i]-average)/rmsd;}else{if(z[i]==average){zi=0.0;}else{zi=z[i];}}
       if (zi>sigma*rmsd){zi=sigma*rmsd;}
       if (zi>max){zi=max;}
       if (zi<min){zi=0.0;}
       scr[i]=zi;
       if (scr[i]>max){scr[i]=max;}
      }
  break;
  case 1:
    for (i=0;i<size;i++){
       zi=Zscore(iter,d_scr,n_scr,dnode,i,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,&gz);
       average=gz.average;
       rmsd=gz.rmsd;
       degree_i=dnode[i].degree;
       if (rmsd>0){zi=(z[i]-average)/rmsd;}else{if(z[i]==average){zi=0.0;}else{zi=z[i];}}
       if (zi>sigma*rmsd){zi=sigma*rmsd;}
       if (zi>max){zi=max;}
       if (zi<min){zi=0.0;}
       scr[i]=zi;
       if (scr[i]>max){scr[i]=max;}
      }
  break;
  default:
    average=rmsdi=rmsd=0;
    for (i=0;i<size;i++){
       zi=Zscore(iter,d_scr,n_scr,dnode,i,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,&gz);
       z[i]=zi*gz.rmsd + gz.average ;
       if (gz.average==0){
          average+=z[i];
          rmsdi+=z[i]*z[i];
         }else{
          average+= gz.average ;
          rmsdi+=gz.rmsd*gz.rmsd + gz.average*gz.average;
         }
      }
    if (size>0)average=average/size;
    if (size>0)rmsdi=rmsdi/size;
    s=rmsdi-average*average;
    if (s>0) rmsd=sqrt(s);
    for (i=0;i<size;i++){
       if (rmsd>0){zi=(z[i]-average)/rmsd;}else{if(z[i]==average){zi=0.0;}else{zi=z[i];}}
       if (zi>sigma*rmsd){zi=sigma*rmsd;}
       if (zi>max){zi=max;}
       if (zi<min){zi=0.0;}
       scr[i]=zi;
       if (scr[i]>max){scr[i]=max;}
      }
  }


  free_fvector(z,0,size);

}
