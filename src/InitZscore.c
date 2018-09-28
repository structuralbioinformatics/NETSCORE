#include "ppin.h"
float  InitZscore(dnode,inode,dedge,size,edges,linker,linker_degree,distribution,ndstr,rnd,maxmin,gz)
int       size,edges,linker,linker_degree,inode,rnd[];
int     **distribution,*ndstr;
float     maxmin[];
node     *dnode;
edge     *dedge;
gauss    *gz;
{
 int   i,j,k,jj,m,n,skip,round,iterations,random_size,rand_type,degree_i,degree_k;
 int   k_rand,k_node,k_bind[MAXI];
 float max,min,sigma,x,y,z,r,*xr;
 float ran0(),ran3();
 int   idum,kdum,max_degree,mind;
 float LocalProduct();
 float *fvector();
 void  free_fvector();
 gauss g;
 gauss normal();

 max=maxmin[0];
 min=maxmin[1];
 sigma=maxmin[2];
 mind=maxmin[3];
 max_degree=maxmin[4];
 iterations=rnd[0];
 rand_type=rnd[1];
 degree_i=dnode[inode].degree;
 /* x=dnode[inode].copy.score; */
 x=0.0; 
 y=x;
 n=0;
 xr=fvector(0,iterations);

 for (j=0;j<degree_i;j++){
   if (dedge[dnode[inode].interact[j]].i == inode){k=dedge[dnode[inode].interact[j]].j;}
   if (dedge[dnode[inode].interact[j]].j == inode){k=dedge[dnode[inode].interact[j]].i;}
   y+=dnode[k].copy.score*LocalProduct(inode,k,dnode)*dedge[dnode[inode].interact[j]].association/linker ;
   //printf("RealScore_add[%d]_from[%d][%d]= %f * %f / %d gives Y= %f\n",j,k,dnode[inode].interact[j],dnode[k].copy.score,dedge[dnode[inode].interact[j]].association,linker,y); 
 }
 if (degree_i>0 && linker_degree>0) y=y/degree_i;

 if (rand_type == 1){
  m=1;
  while(n<iterations-1){
    xr[n]=x;
    for (j=0;j<degree_i;j++){k_bind[j]=0;}
    j=0;
    skip=0;
    round=0;
    while (skip==0 && j < degree_i ) {
       if (dedge[dnode[inode].interact[j]].i == inode){k=dedge[dnode[inode].interact[j]].j;}
       if (dedge[dnode[inode].interact[j]].j == inode){k=dedge[dnode[inode].interact[j]].i;}
       if (dnode[k].degree<MAXD-1){degree_k=dnode[k].degree;}else{degree_k=MAXD-1;}
       if (ndstr[degree_k]<max_degree){random_size=ndstr[degree_k];}else{random_size=max_degree-1;}
       idum   = (degree_i + 1)  * m ;
       kdum   = (int) ( 1000.0 * ran3(&idum) );
       r=ran3(&kdum); 
       k_node = (int) ( random_size * r );
       if (k_node>=max_degree){printf("Knode %d Random_size %d Idum %d kdum %d randm %f\n",k_node,random_size,idum,kdum,r);nrerror("EMERGENCY STOP\n");}
       k_rand = distribution[degree_k][k_node];
       for (jj=0;jj< j;jj++){ 
          if (k_rand == k_bind[jj]){skip=1; round++ ; break;}
       }
       if (round>10) { k_bind[j]=k_rand; j++; round=0;}
       if (skip==0 ) { k_bind[j]=k_rand; j++; round=0;}else{skip=0;}
       /*printf("N %d J %d R_SIZE %d Degree_K %d  Degree_i %d k_node %d K_rand %d Skip %d Round %d Name_k %s\n",n,j,random_size,degree_k,degree_i,k_node,k_rand,skip,round,dnode[k_rand].name1);*/
       m++;
       if (m>1000){m=1;}
    }
    for  (j=0;j<degree_i;j++){
       k_rand = k_bind[j];
       xr[n] += dnode[k_rand].copy.score*LocalProduct(inode,k_rand,dnode)*dedge[dnode[inode].interact[j]].association/linker ;
/*       printf("Random Score_J[%d] Add_from[%d][%d] = Node %f * Edge %f / Linker %d  is %f\n",j,k_rand,dnode[inode].interact[j],dnode[k_rand].copy.score,dedge[dnode[inode].interact[j]].association,linker,xr[n]); */
    }
    if (degree_i>0 && linker_degree>0) xr[n]=xr[n]/degree_i;
/*    printf("IterationR %d Score[%d] %f Random %f Name1 %s Name2 %s \n",n,inode,y,xr[n],dnode[inode].name1,dnode[inode].name2);
    printf("K_RAND:\n");
    for  (j=0;j<degree_i;j++){printf("%d ",k_bind[j]);}
    printf("\n");
 */
    n++;
  }
 }else{
  m=1;
  while(n<iterations-1){
    xr[n]=x;
    for (j=0;j<degree_i;j++){
      idum   = (degree_i + 1)  * m ;
      kdum   = (int) ( 1000.0 * ran3(&idum) );
      r=ran3(&kdum);
      k_rand = (int) ( (size-1) * r);
      if (k_rand>=size){printf("Knode %d size %d Idum %d kdum %d randm %f\n",k_rand,size,idum,kdum,r);nrerror("EMERGENCY STOP\n");}
      xr[n]+= dnode[k_rand].copy.score*LocalProduct(inode,k_rand,dnode)*dedge[dnode[inode].interact[j]].association/linker ;
      m++;
      if (m>1000){m=1;}
    }
    if (degree_i>0 && linker_degree>0) xr[n]=xr[n]/degree_i;
    n++;
  }
  /*  printf("IterationR %d Score %f Random %f\n",n,y,xr[n]);*/
 }
 //printf("\n");

 if (iterations>0){
  xr[iterations-1]=y;
  g = normal(xr,iterations);
  if (g.rmsd > 0) { z = ( y - g.average ) / g.rmsd ; }else{ if (y==g.average){z=0.0;g.rmsd=1.0;}else{z=y;g.average=0;g.rmsd=1;}}
  *gz = g;
 }else{
  g.average=0;
  g.rmsd=1;
  *gz=g;
  z=y;
 } 

/*
 if (z>sigma*g.rmsd){z=sigma*g.rmsd;}
 if (z>max){z=max;}
 if (z<min){z=0;}
 if (degree_i>MAXD){z=MAXD*z/degree_i;}
*/

 free_fvector(xr,0,iterations);

 return z;
}


