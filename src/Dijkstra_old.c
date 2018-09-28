#include "ppin.h"

void Dijkstra(A,N,dist,s)
int  **A,*dist,N,s;
{
  int i, j, u, v, c ;
  int *PQ;
  int *ivector();
  void free_ivector(),swap();
  
  for (i=0;i<N;i++){ printf("Source=%d Shortestpath[%d]=%d\n",s,i,dist[i]);}

  PQ=ivector(0,N);

  for(v = 0; v < N; v++){
    if(v == s){
      dist[v] = 0;
    }else{
      dist[v] = N;
    }
    PQ[v] = dist[v];
    printf("test V=%d\n",v);
  }
  for (i=0;i<N;i++){ printf("New Shortestpath[%d]=%d\n",i,dist[i]);}
  i = N-1;
  while(i > 0){
    for(j = 0; j < i; j++){if(PQ[j] < PQ[i]) swap(PQ,i,j);}
    u = PQ[i];
    for(v = 0; v < N, A[u][v] != 0; j++){
      c = dist[u] + A[u][v];
      if(c < dist[v]) dist[v] = c;
    }
    i = i-1;
  }
 
  free_ivector(PQ,0,N);

}

void swap(a,i,j)
int *a,i,j;
{
 int x;
 x=a[i];
 a[i]=a[j];
 a[j]=x;
}

