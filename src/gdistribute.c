#include "ppin.h"
void gdistribute(d,nd,maxk,maxx,min)
int **d;
int  *nd;
int   maxk,maxx,min;
{
 int  i,k,m,n,nn,j,jj,*ini,dim;
 int *ivector();
 void free_ivector();

   for (k=0;k<maxk;k++){
   dim=nn=0;
   while(dim<min && (k+nn)<maxk){dim+=nd[k+nn];nn++;}
   if (dim>maxx){dim=maxx;}
   m=0;
   ini=ivector(0,nn);
   ini[0]=nd[k];
   for (j=1;j<nn;j++){
     ini[j]=ini[j-1]+nd[k+j];
     for (jj=0;jj<nd[k+j];jj++){
         if ( (nd[k]+m)<dim){
           d[k][nd[k]+m]=d[k+j][jj];
           m++;
         }else{
           break;
         }
     }
   }
   nd[k]=dim;
   for (j=1;j<nn;j++){
     m=0;
     for (jj=0;jj<dim;jj++){
        if (jj<ini[j-1] || jj>=ini[j]){
         d[k+j][nd[k+j]+m]=d[k][jj];
         m++;
        }
     }
     nd[k+j]=dim;
   }
   k=k+nn-1; 
   free_ivector(ini,0,nn);
  }

}
