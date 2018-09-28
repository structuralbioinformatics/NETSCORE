#include "ppin.h"

void   RandomSwapNodes(size2,interaction,protein,distribution,ndstr,max_degree,rnd,swap)
int       **distribution,*ndstr;
int       *size2;
node      *protein;
edge      *interaction;
int       max_degree;
int       rnd[],**swap;
{
 int   i,j,k,jj,m,n,skip,round,iterations,random_size,rand_type,degree_i,degree_k,size,inode;
 int   k_rand,k_node,k_bind[MAXI];
 float ran0(),ran3(),free_imatrix();
 int   iterswap,idum,kdum,mind,**imatrix();

 iterations=rnd[0];
 rand_type=rnd[1];
 size=size2[0];
 check=imatrix(0,size,0,iterations);

 for (i=0;i<size;i++){
 for (j=0;j<iterations;j++){
    check[i][j]=0;
 }}

 printf("RandomSwap  %d x %d\n",size,iterations);
 
 if (rand_type == 1){
  for (inode=0;inode<size;inode++){
      degree_i=protein[inode].degree;
      if (degree_i>=MAXD-1) degree_i=MAXD-1;
      //printf("Inode=%d Degree = %d NDSTR= %d \n",inode,degree_i,ndstr[degree_i]);
      n=0;
      while(n<iterations){
         if (check[inode][n] == 0 ){
           if (ndstr[degree_i]<max_degree){random_size=ndstr[degree_i];}else{random_size=max_degree-1;}
           if (random_size>=1){
             idum   = (degree_i + 1)  * (n + 1) ;
             kdum   = (int) ( 1000.0 * ran3(&idum) );
             k_node = (int) ( random_size * ran3(&kdum) );
             if (k_node>random_size){printf("Error K_NODE %d > %d\n",k_node,random_size);}
             if (k_node<0)k_node=0;
             k_rand = distribution[degree_i][k_node];
             if (k_rand==inode){
               if(k_node>0)                   {k_rand=distribution[degree_i][k_node-1];}
               else if((k_node+1)<random_size){k_rand=distribution[degree_i][k_node+1];}
               }
             iterswap=0;
             while (check[k_rand][n] == 1 && iterswap < 100 && iterswap < random_size){
                idum   = (degree_i + 1)  * (n + 1) ;
                kdum   = (int) ( 1000.0 * ran3(&idum) );
                k_node = (int) ( random_size * ran3(&kdum) );
                if (k_node>random_size){printf("Error K_NODE %d > %d\n",k_node,random_size);}
                if (k_node<0)k_node=0;
                k_rand = distribution[degree_i][k_node];
                if (k_rand==inode){
                  if(k_node>0)                   {k_rand=distribution[degree_i][k_node-1];}
                  else if((k_node+1)<random_size){k_rand=distribution[degree_i][k_node+1];}
                  }
                if (check[k_rand][n] == 0){
                  swap[inode][n]=k_rand;
                  swap[k_rand][n]=inode;
                  check[inode][n]=1;
                  check[k_rand][n]=1;
                }
                iterswap++;
             }
             if (iterswap==100 || iterswap == random_size){swap[inode][n]=inode;}
             n++;
            }else{
             //printf("swap[Inode=%d][Iteration=%d]= %d Random_size=%d \n",inode, n, swap[inode][n]);
             swap[inode][n]=inode;
             n++;
            }
          }
        }
      }
 }else{
  for (inode=0;inode<size;inode++){
      n=0;
      degree_i=protein[inode].degree;
      if (degree_i>=MAXD-1) degree_i=MAXD-1;
      while(n<iterations){
          if (check[inode][n] == 0 ){
           idum   = (degree_i + 1)  * (n + 1) ;
           kdum   = (int) ( 1000.0 * ran3(&idum) );
           k_node = (int) ( (size-1) * ran3(&kdum) );
           if (k_node>size){printf("Error K_NODE %d > %d\n",k_node,size);}
           if (k_node<0)k_node=0;
           k_rand = k_node;
           if (k_rand==inode){
               if(k_node>0)            {k_rand=k_node-1;}
               else if((k_node+1)<size){k_rand=k_node+1;}
               }
           iterswap=0;
           while (check[k_rand][n] == 1 && iterswap < 100 && iterswap < random_size){
              idum   = (degree_i + 1)  * (n + 1) ;
              kdum   = (int) ( 1000.0 * ran3(&idum) );
              k_node = (int) ( (size-1) * ran3(&kdum) );
              if (k_node>size){printf("Error K_NODE %d > %d\n",k_node,size);}
              if (k_node<0)k_node=0;
              k_rand = k_node;
              if (k_rand==inode){
                  if(k_node>0)            {k_rand=k_node-1;}
                  else if((k_node+1)<size){k_rand=k_node+1;}
                  }
              if (check[k_rand][n] == 0){
               swap[inode][n]=k_rand;
               swap[k_rand][n]=inode;
               check[inode][n]=1;
               check[k_rand][n]=1;
               }
              iterswap++;
            }
           if (iterswap==100 || iterswap == random_size){swap[inode][n]=inode;}
           n++;
          }
        }
      }
 }

 free_imatrix(check,0,size,0,iterations);

}
