#include "ppin.h"
float CumulativeMatrix(iter,maxmin,indirect,linker,i,dnode,dedge,size,d_old,d_scr,n_scr,dim,mind,escr)
int       i,size,iter,linker,indirect;
float     maxmin[];
int      *n_scr,*dim;
gradient **d_old,**d_scr;
float     mind,*escr;
node     *dnode;
edge     *dedge;
{
 float e,max,min,weight;
 int j,k,kk,jj,skip;
 float Transform();
 void  nrerror();

  max=maxmin[0];
  min=maxmin[1];

   for (j=0;j<dnode[i].degree;j++){
   if (dedge[dnode[i].interact[j]].i == i){k=dedge[dnode[i].interact[j]].j;}
   if (dedge[dnode[i].interact[j]].j == i){k=dedge[dnode[i].interact[j]].i;}
     for (kk=0;kk<dim[i];kk++){
       if (d_scr[i][kk].j== k && d_old[i][kk].size>min){
         if (escr[dnode[i].interact[j]]>mind){
           weight=dedge[dnode[i].interact[j]].association+escr[dnode[i].interact[j]];
           if (dedge[dnode[i].interact[j]].association>0){
            weight=weight/dedge[dnode[i].interact[j]].association;
           }
         }else{
           weight=1.0;
         }
         d_scr[i][kk].size=weight*d_old[i][kk].size;
         if (d_scr[i][kk].size>min && weight!=1.0){
         printf("FIX Original Flow[%d] I[%d]J[%d] %f X %f = %f\n",kk,i,k,weight,d_old[i][kk].size,d_scr[i][kk].size);
         }
         if (d_scr[i][kk].size>max){d_scr[i][kk].size=max;}
       }
     }
     for (jj=0;jj<dim[k];jj++){
        if (d_old[k][jj].j!=i && d_old[k][jj].size>min && d_old[k][jj].j<size && d_old[k][jj].j>=0){
          skip=0;
          if (indirect>0){
          for (kk=0;kk<n_scr[i];kk++){
             if (skip==1){break;}
             if (d_scr[i][kk].j==d_old[k][jj].j){
                if ( d_old[k][jj].size/pow(linker,(iter)) > d_scr[i][kk].size){d_scr[i][kk].size=d_old[k][jj].size/pow(linker,(iter));}
                if (d_scr[i][kk].size>max){d_scr[i][kk].size=max;}
                skip=1;
             }
           }}
           if (skip==0){
            if (n_scr[i] == MAXWN){printf("Warning Dimension Node/edge SCR[%d] (%d) > %d\n",i,n_scr[i],MAXWN);}
            if (n_scr[i] == MAXWE){printf("Error Dimension Node/edge SCR[%d] (%d) exceeded %d, change MAXWE\n",i,n_scr[i],MAXWE);nrerror("STOP\n");}
            if (n_scr[i] == size){printf("Error Dimension Node/edge SCR[%d] (%d) exceeded SIZE %d, use -dxi larger than 0\n",i,n_scr[i],size);
                                  for (kk=0;kk<n_scr[i];kk++){printf("d_scr[%d] I[%d]J[%d] \n",kk,i,d_scr[i][kk].j);};
                                  nrerror("STOP\n");}
            kk=n_scr[i];
            d_scr[i][kk].j=d_old[k][jj].j;
            d_scr[i][kk].size=d_old[k][jj].size/pow(linker,(iter));
            if (d_scr[i][kk].size>max){d_scr[i][kk].size=max;}
            n_scr[i]=n_scr[i]+1;
          } 
        }
     }
  }

  e=Transform(i,maxmin,dim,d_old);

 return e;
}
