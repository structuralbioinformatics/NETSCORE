#include "ppin.h"

void WDijkstra(A,N,dist,s)
float  **A,*dist;
int  N,s;
{
  int i, j, u, v, count ;
  float  k;
  int *ivector(),*mark,*predecessor, Wminimum();
  void free_ivector(),Wprintpath();
 
 mark=ivector(0,N);
 predecessor=ivector(0,N);
 k=(float)N/MINSHORT;


  for(v = 0; v < N; v++){
    mark[v]=0;
    predecessor[v]=-1;
    if(v == s){
      dist[v] = 0.0;
    }else{
      dist[v] = k;
    }
  }
  count=0;
  while(count<N){
    u=Wminimum(dist,mark,N);
    mark[u]=1;
    for(i=0;i<N;i++)
    {
     if(A[u][i]>0) {
      if(mark[i]!=1) {
      if(dist[i]>=dist[u]+A[u][i]) {
       dist[i]=dist[u]+A[u][i];
       predecessor[i]=u;
      }}
     }
    }
    count++;
   }

/*
   for(i=0;i<N;i++)
   {
   Wprintpath(s,i,predecessor);
   if(dist[i]!=k) printf("->(%f)\n",dist[i]);
   }
*/

   free_ivector(mark,0,N);
   free_ivector(predecessor,0,N);

}

int Wminimum(a,m,k)
float *a;
int *m;
int k;
{
 float mi;
 int i,t;
 mi=(float)k/MINSHORT;
 for(i=0;i<k;i++)
 {
  if(m[i]!=1)
  {
  if(mi>=a[i])
  {
   mi=a[i];
   t=i;
  }
 }}
 return t;
}

void Wprintpath(x,i,p)
int  x,i,*p;
{
 printf("\n");
 if(i==x)
 {
 printf("%d",x);
 }else if(p[i]<0){
 printf("no path from %d to %d",x,i);
 }else{
 printpath(x,p[i],p);
 printf("..%d",i);
 }
}
