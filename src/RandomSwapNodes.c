#include "ppin.h"

void   RandomSwapNodes(size2,interaction,protein,distribution,ndstr,max_degree,rnd,swap)
int       **distribution,*ndstr;
int       *size2;
node      *protein;
edge      *interaction;
int       max_degree;
int       rnd[],**swap;
{
 int   i,j,k,n,iteration,random_size,rand_type,degree_i,size,inode,dim;
 int   *list;
 int   *ivector();
 float ran3();
 int   idum,init;
 void  free_ivector();
 void shuffle();

 iteration=rnd[0];
 rand_type=rnd[1];
 size=size2[0];

 //printf("RandomSwap  %d x %d Max-degree %d Max-iteration %d Type %d\n",size,iteration,max_degree,iteration,rand_type);
 
 for (inode=0;inode<size;inode++){
 for (n=0;n<iteration;n++){
     swap[inode][n]=inode;
     }}

 //printf("Got space for SWAP\n");

 if (rand_type == 1){
  list=ivector(0,max_degree);
  for (n=0;n<iteration;n++){
    for (j=1;j<MAXD;j++){
	//printf("Degree: %d\n",j);
	//printf("NDSTR: %d\n",ndstr[j]);
        if (ndstr[j]<max_degree){random_size=ndstr[j];}else{random_size=max_degree;}
        for (k=0;k<random_size;k++){list[k] = distribution[j][k];}
        idum=j*(n+1)*(random_size+1);
        init=(int) (idum*ran3(&idum));
        shuffle(list,random_size,init);
        for (k=0;k<random_size;k++){
		//printf("Iteration %d Random %d degree %d\n",n,k,j);
		//printf("Selected in distribution %d\n",distribution[j][k]);
		//printf("List %d\n",list[k]);
		swap[distribution[j][k]][n] = list[k];
	       }
        }
     }
  //printf("Free LIST\n");
  free_ivector(list,0,max_degree);
  //printf("Done\n");
 }else{
  list=ivector(0,size);
  for (n=0;n<iteration;n++){
      for (k=0;k<size;k++){list[k]=k;}
      idum=(n+1)*(size+1); 
      init=(int) (idum*ran3(&idum));
      shuffle(list,size,init);
      for (k=0;k<size;k++){swap[k][n]=list[k];}
     }
  free_ivector(list,0,size);
 }


}


void shuffle(theArr,size,init)
int *theArr;
int size;
int init;
{
	   int temporary, randomNum, last, kdum,idum;
           float ran3(),test;
	   for (last = size; last > 1; last--)
	   {
              idum   = last*init;
              test   = ran3(&idum);
              kdum   = (int) (test  * last);
              test   = ran3(&kdum);
	      randomNum =  (int) (test * last) % last;
	      temporary = theArr[randomNum];
	      theArr[randomNum] = theArr[last - 1];
	      theArr[last - 1] = temporary;
	   }
}
// end shuffleElements( )

