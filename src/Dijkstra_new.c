#include "ppin.h"

void Dijkstra(A,N,dist,s)
int  **A,*dist,N,s;
{
  int i, j, u, v, count ;
  int *ivector(),*mark,*predecessor;
  void free_ivector(),swap();
 
 mark=ivector(0,N);
 predecessor=ivector(0,N);


  for(v = 0; v < N; v++){
    mark[j]=0;
    predecessor[j]=0;
    if(v == s){
      dist[v] = 0;
    }else{
      dist[v] = N;
    }
  }
  count=0;
  while(count<N){
    u=minimum(dist,mark,N);
    mark[u]=1;
    for(i=0;i<N;i++)
    {
     if(A[u][i]>0)
     {
      if(mark[i]!=1)
      {
      if(dist[i]>dist[u]+A[u][i])
      {
       dist[i]=dist[u]+A[u][i];
       predecessor[i]=u;
      }
      }
     }
    }
    count++;
   }
   for(i=0;i<N;i++)
   {
   printpath(s,i,predecessor);
   if(dist[i]!=N) printf("->(%d)\n",dist[i]);
   }
   free_ivector(mark,0,N);
   free_ivector(predecessor,0,N);

}

int minimum(a,m,k)
int *a;
int *m;
int k;
{
 int mi;
 int i,t;
 mi=k;
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
void printpath(x,i,p)
int  x,i,*p;
{
 printf("\n");
 if(i==x)
 {
 printf("%d",x);
 }
 else if(p[i]==0)
 printf("no path from %d to %d",x,i);
 else
 {
 printpath(x,p[i],p);
 printf("..%d",i);
 }
}
