#include "ppin.h"
// DensitySearch: total numebr of groups
// **s: scoring matrix of similarity
// **f: filter to check distance similarities (0 exclude, 1 include)
// n: size of **s is n X n
// **g: groups of indexes that were clustered
// d[i]: size of g[i]
// p[0]: minimum score to link a putative new member to a group	
// p[1]: maximum deviation of the mean produced by the inclusion of a new member
// p[2]: minimum ratio of common elements of the smallest group to be merged (0 forces all groups to have null intersection)
// p[3]: Maximum of ratio of deviation of the mean to accept merged groups
// p[4]: Minimum cluster mean cut-off to stop cluster growth (0 forces infinite growth)
// 
int DensitySearch(s,f,n,g,d,p)
float **s;
int   **f;
int     n;
int   **g;
int    *d;
float  *p;
{
 int    i,j,ng,merge; 
 int    MergeGroups();
 int   *skip;
 int   **cluster,*dcluster;
 int   *ivector(),**imatrix();
 void   free_ivector(),free_imatrix();
 void   nrerror(),InitGroup(),ExpandGroup();

 
 skip=ivector(0,n);
 cluster=imatrix(0,MAXXG,0,MAXXSG);
 dcluster=ivector(0,MAXXG);
 for (i=0;i<n;i++){skip[i]=1;}
 ng=0;
 while (sum(skip,n)>0){
   for (i=0;i<n;i++){
   if (skip[i]==1){
      InitGroup(i,g,ng,d,skip);
      ExpandGroup(i,g,ng,skip,s,f,n,d,p);
      ng++;
      if (ng>MAXXG)nrerror("Too many groups in Density Search Clustering, increase MAXXG");
   }}
 }
 if (p[2]!=0)merge=1;
 else        merge=0;
 while (merge>0){merge=MergeGroups(g,ng,d,cluster,dcluster,s,f,n,p);ng=ng-merge;}
 free_imatrix(cluster,0,MAXXG,0,MAXXSG);
 free_ivector(dcluster,0,MAXXG);
 free_ivector(skip,0,n);
 return ng;
}

int MergeGroups(g,ng,d,cluster,dcluster,s,f,n,p)
int     ng;
int   **g;
int   **cluster;
float **s;
int   **f;
int     n;
int    *d;
int    *dcluster;
float  *p;
{
 int   i,j,k,m,nn,ii,jj,*skip,kk,redundant,**joint;
 float c,Common();
 float x,y,z,Average(),JointAverage();
 int   **imatrix(),*ivector();
 void  free_ivector(),free_imatrix();

 skip=ivector(0,ng);
 joint=imatrix(0,ng,0,2);

 m=0;
 for (i=0;i<ng;i++){skip[i]=1;}
 for (i=0;i<ng;i++){
 for (j=i+1;j<ng;j++){
   c=Common(g,i,j,d);
   if ((c/d[i]) > (c/d[j])){x=(c/d[i]);}else{x=(c/d[j]);}
   y=JointAverage(g,i,j,s,d)-Average(g,i,s,d);
   y=sqrt(y*y)/Average(g,i,s,d);
   z=JointAverage(g,i,j,s,d)-Average(g,j,s,d);
   z=sqrt(z*z)/Average(g,j,s,d);
   //printf("Merging groups %d & %d => %f commons (Ratio=%f Deviation= (%f , %f) )\n",i,j,c,x,y,z);
   if (z>y)y=z;
   if ( x>=p[2] && y<=p[3] && (d[i]+d[j]-(int)c)<MAXXSG && skip[i]>0 && skip[j]>0){
    joint[m][0]=i;
    joint[m][1]=j;
    skip[i]=0;
    skip[j]=0;
    //printf("Accepted merge Joint[%d] %d and %d\n",m, joint[m][0],joint[m][1]);
    m++;
   }
 }}

     /*
     printf("\nInitial Clusters\n");
     for (i=0;i<ng;i++){
        printf("Group[%d]:\t",i);
        for (k=0;k<d[i];k++){
             printf("%5d ",g[i][k]);
             }
         printf("\n");
         }
      */

 nn=0;
 dcluster[nn]=0;
 for (k=0;k<m;k++){
   ii=joint[k][0];
   jj=joint[k][1];
   for (i=0;i<d[ii];i++){
    cluster[nn][dcluster[nn]]=g[ii][i];
    dcluster[nn]+=1;
    }
   for (j=0;j<d[jj];j++){
     redundant=0;
     for (kk=0;kk<dcluster[nn];kk++){
         if (cluster[nn][kk]==g[jj][j]){redundant++;break;}}
     if (redundant==0){
       cluster[nn][dcluster[nn]]=g[jj][j];
       dcluster[nn]+=1;
     }
    }
   nn++;
   dcluster[nn]=0;
 }
 for (i=0;i<ng;i++){
 if (skip[i]==1){
    dcluster[nn]=d[i];
    for (j=0;j<d[i];j++){
      cluster[nn][j]=g[i][j];
    }
    nn++;
 }}

     /*
     printf("\nNew Clusters\n");
     for (i=0;i<nn;i++){
        printf("Group[%d]:\t",i);
        for (k=0;k<dcluster[i];k++){
             printf("%5d ",cluster[i][k]);
             }
         printf("\n");
         }
     */


 for (i=0;i<nn;i++){
 d[i]=dcluster[i];
 for (j=0;j<d[i];j++){
   g[i][j]=cluster[i][j];
 }} 
 for (i=nn;i<ng;i++){d[i]=0;}
 

 free_ivector(skip,0,ng);
 free_imatrix(joint,0,ng,0,2);
 return  m;
}

float  Common(g,i,j,d)
int   **g;
int     i,j;
int    *d;
{
 int   k,m;
 float x;
 x=0.0;
 for (k=0;k<d[i];k++){
 for (m=0;m<d[j];m++){
   if (g[i][k]==g[j][m])x++; 
 }}
 return  x;
}

void ExpandGroup(init,g,ng,skip,s,f,n,d,p)
int     init,ng;
int   **g;
int   *skip;
float **s;
int   **f;
int     n;
int    *d;
float  *p;
{
 int   **fg,check,i,j;
 int     TestGroup();
 int   **imatrix();
 void    free_imatrix();

  fg=imatrix(0,n,0,n);

  for (i=0;i<n;i++){
   for (j=i;j<n;j++){
     fg[i][j]=f[i][j];
     fg[j][i]=f[i][j];
    }}

/*
  if (p[2]==0.0){
    for (i=0;i<n;i++){
      if (skip[i]==0 && i!=init){
        for (j=i;j<n;j++){
           fg[i][j]=fg[j][i]=0;
        }}}
  }
*/

  check=0;
  while (check==0){ check=TestGroup(g,ng,skip,s,fg,n,d,p); } 

  free_imatrix(fg,0,n,0,n);

}

int TestGroup(g,ng,skip,s,f,n,d,p)
int     ng;
int   **g;
int    *skip;
float **s;
int   **f;
int     n;
int    *d;
float  *p;
{
 int    *nn,check;
 int    *nearest2group();
 int    *ivector();
 float   max,mean,new,diff,d1,d2;
 float   Average(),NewAverage(),NewDiff();
 void    free_ivector(),nrerror();

  check=0;
  new=0.0;
  mean=0.0;
  nn=ivector(0,2);
  nn=nearest2group(g,ng,s,f,n,d,skip,p);
  check=1;
  if (nn[1]>=0){
   max=s[nn[0]][nn[1]];
   if (d[ng]>1){
    mean=Average(g,ng,s,d);
    new=NewAverage(g,ng,s,d,nn);
    d1=sqrt((new-mean)*(new-mean));
    d2=NewDiff(g,ng,s,d,nn);
    if (d1<d2) diff=d2; else diff=d1;
   }else{
    diff=0.0;
    mean=max;
    new=max;
   }
   if (max>=p[0] && diff<=p[1] && (p[4]==0 || new>=p[4]) ){
     skip[nn[1]]=0;
     g[ng][d[ng]]=nn[1];
     printf("New member of G[%d] is G[%d][%d]=%d (Old Mean: %f & New Mean: %f differ by %f (D1 %f D2 %f)\n",ng,ng,d[ng],g[ng][d[ng]],mean,new,diff,d1,d2);
     d[ng]++;
     if (d[ng]>MAXXSG)nrerror("Too many elements of a group in Density Search Clustering, increase MAXXSG");
     check=0;
   }
  }
  free_ivector(nn,0,2);

  return check;

} 
void InitGroup(i,g,ng,d,skip)
int     i;
int     ng;
int   **g;
int    *d;
int   *skip;
{
   g[ng][0]=i;
   //printf("First member of G[%d] is G[%d][%d]=%d\n",ng,ng,d[ng],g[ng][d[ng]]);
   d[ng]=1; 
   skip[i]=0;
}

float Average(g,n,s,d)
int     n;
int   **g;
float **s;
int    *d;
{
 float  x;
 int    i,j,m;
  m=0;
  x=0.0;
  for (i=0;i<d[n];i++){
    for (j=i+1;j<d[n];j++){
      x+=s[g[n][i]][g[n][j]];
      m++;
    }}
  x=x/m;
  return x;
}

float JointAverage(g,n,m,s,d)
int     n,m;
int   **g;
float **s;
int    *d;
{
 float  x;
 int    i,j,k;
  k=0;
  x=0.0;
  for (i=0;i<d[n];i++){
    for (j=i+1;j<d[n];j++){
      x+=s[g[n][i]][g[n][j]];
      k++;
    }}
  for (i=0;i<d[m];i++){
    for (j=i+1;j<d[m];j++){
      x+=s[g[m][i]][g[m][j]];
      k++;
    }}
  for (i=0;i<d[n];i++){
  for (j=0;j<d[m];j++){
      x+=s[g[n][i]][g[m][j]];
      k++;
  }}
  x=x/k;
  return x;
}


float NewAverage(g,ng,s,d,nn)
int     ng;
int   **g;
float **s;
int    *d;
int    *nn;
{
 float  x;
 int    i,j,n,new;
  new=nn[1];
  n=0;
  x=0.0;
  for (i=0;i<d[ng];i++){
  for (j=i+1;j<d[ng];j++){
      x+=s[g[ng][i]][g[ng][j]];
      n++;
  }}
  for (i=0;i<d[ng];i++){
      x+=s[g[ng][i]][new];
      n++;
  }
  x=x/n;
  return x;
}

float NewDiff(g,ng,s,d,nn)
int     ng;
int   **g;
float **s;
int    *d;
int    *nn;
{
 float  x,y,z,dx,mx,dy,my,sx,sy;
 int    i,j,n,m,new;
  new=nn[1];
  x=dx=mx=0.0;
  y=dy=my=0.0;
  z=0.0;
  n=0;
  m=0;
  for (i=0;i<d[ng];i++){
  for (j=i+1;j<d[ng];j++){
      if (s[g[ng][i]][g[ng][j]]>x) x=s[g[ng][i]][g[ng][j]];
      dx+=s[g[ng][i]][g[ng][j]]*s[g[ng][i]][g[ng][j]];
      mx+=s[g[ng][i]][g[ng][j]];
      dy=dx;
      my=mx;
      n++;
  }}
  y=x;
  m=n;
  for (i=0;i<d[ng];i++){
      if (s[g[ng][i]][new]>y)y=s[g[ng][i]][new];
      dy+=s[g[ng][i]][new]*s[g[ng][i]][new];
      my+=s[g[ng][i]][new];
      m++;
  }
  sx=dx/n - mx*mx;
  sy=dy/m - my*my;
  if (sx>0) sx=sqrt(sx);
  if (sy>0) sy=sqrt(sy);
  if ((y-x)>(sy-sx)) z=y-x;else z=sy-sx;
  
  return z;
}


int sum(a,n)
int  *a;
int   n;
{
 int i,s;
 s=0;
 for(i=0;i<n;i++){s+=a[i];}
 return s; 
}

int nearest_neighbor(i,s,f,n,skip,p)
int     i;
float **s;
int   **f;
int     n;
int    *skip;
float  *p;

{
 int    j,k;
 float  max;
  max=0.0;
  j=-1;
  for (k=0;k<n;k++){
     if (s[i][k]>max && f[i][k]==1 && i!=k && (p[2]>0 || skip[k]>0)){max=s[i][k];j=k;}
  }
 return j;
}

int *nearest2group(g,ng,s,f,n,d,skip,p)
int   **g,*d;
int     n,ng;
float **s;
int   **f;
int    *skip;
float  *p;
{
 int     i,j,k,m,h,ne,*nn;
 int   **fg;
 int     nearest_neighbor();
 int   **imatrix(),*ivector();
 float   max;
 void    free_imatrix();

  fg=imatrix(0,n,0,n);
  nn=ivector(0,2);
  ne=d[ng];
  max=0.0;
  nn[0]=nn[1]=-1;
  for (i=0;i<n;i++){
  for (j=i;j<n;j++){
   fg[i][j]=f[i][j];
   fg[j][i]=f[i][j];
  }}
  for (i=0;i<ne;i++){
  for (j=i;j<ne;j++){
   m=g[ng][i];
   h=g[ng][j];
   fg[m][h]=fg[h][m]=0;
  }}
  for (k=0;k<ne;k++){
    i=g[ng][k];
    j=nearest_neighbor(i,s,fg,n,skip,p);
    //printf("Check Nearest to G[%d][%d]=%d is %d \n",ng,k,i,j);
    if (j>=0){
     if (s[i][j]>max && f[i][j]==1 ){
       max=s[i][j];
       nn[0]=i;
       nn[1]=j;
     }}
  }
/*
  if (nn[0]>=0 && nn[1]>=0){
   printf("Final nearest to G[%d] = %d (closest from G is %d, distance= %e)\n",ng,nn[1],nn[0],s[nn[0]][nn[1]]);
  }else{
   printf("Final nearest to G[%d]  DOES NOT EXIST\n",ng);
  }
*/
  free_imatrix(fg,0,n,0,n);

 return nn;
}


