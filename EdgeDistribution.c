#include "ppin.h"
void    EdgeDistribution(dnode,dsize,edistribution,nedstr,maxmin)
node  *dnode;
int    dsize;
int  **edistribution,*nedstr;
float  maxmin[];
{
 int  i,j,k,degree,dim, idum, kdum, m, max_degree,min_dstr,maxx;
 float r;
 float ran0(),ran3();
 void gdistribute();

 max_degree=0;
 for (k=0;k<MAXD;k++){nedstr[k]=1;}
 for (i=0;i<dsize;i++){
   degree=(int)(dnode[i].degree/5);
   m=1;
   if (degree<MAXD-1){
    if ( (degree+1) > max_degree){max_degree=degree+1;}
    if (nedstr[degree]>MAXIDE-1){
         idum   = (degree)  * m ;
         kdum   = (int) ( 1000.0 * ran3(&idum) );
         r      =  ran3(&kdum);
         if (r>1) {r = r - (int)(r);}
         dim    =  (int) ((MAXIDE-1) *  r );
	 //printf("1-EdgeDistribution random selection: idum %d Kdum %d R %f DIM %d\n",idum,kdum,r,dim);
    }else{
         dim=nedstr[degree];
    }
    edistribution[degree][dim]=i;
    nedstr[degree]++;
    if (dim>nedstr[degree]){printf("Error on DIM vs NEDSTR\n");}
    //printf("1-EdgeDsitribution size (%d) %d is larger than MAXIDE (Edistr[%d][%d]=%d\n",degree,nedstr[degree],degree,dim,i);
    //printf("1-Edistr[%d][%d]=%d NEDSTR=%d Real degree %d\n",degree,dim,i,nedstr[degree],dnode[i].degree);
   }else{
    max_degree=MAXD;
    if (nedstr[MAXD-1]>MAXIDE-1){
         idum   = (degree)  * m ;
         kdum   = (int) ( 1000.0 * ran3(&idum) );
         r      =  ran3(&kdum);
         if (r>1) {r = r - (int)(r);}
         dim    =  (int) ((MAXIDE-1) * r );
	 //printf("2-EdgeDistribution random selection: idum %d Kdum %d R %f DIM %d\n",idum,kdum,r,dim);
    }else{
         dim=nedstr[MAXD-1];
    }
    edistribution[MAXD-1][dim]=i;
    nedstr[MAXD-1]++;
    //printf("2-EdgeDsitribution size (%d) %d is larger than MAXIDE (Edistr[%d][%d]=%d\n",MAXD-1,nedstr[MAXD-1],MAXD-1,dim,i);
    //printf("2-Edistr[%d][%d]=%d NEDSTR=%d Real degree %d \n",degree,dim,i,nedstr[MAXD-1],dnode[i].degree);
   }
 }

  min_dstr=(int)maxmin[8];
  maxx=MAXIDE;
  gdistribute(edistribution,nedstr,max_degree,maxx,min_dstr);


}
