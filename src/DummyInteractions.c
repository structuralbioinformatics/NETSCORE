#include "ppin.h"
int     DummyInteractions(fcopy,flocal,protein,size)
char    fcopy[MAXS],flocal[MAXS];
node   *protein;
int     size;
{
  FILE       *INP;
  char        proteome[MAXP][MAXNAME];
  int         i,edges,proteosize;
  abundance   copy[MAXP];
  int         ReadNode();
  abundance   AssignNode();
  void        nrerror();


  proteosize=ReadNode(fcopy,flocal,copy,proteome);

  for (i=0;i<proteosize;i++){
    strcpy(protein[i].name1,proteome[i]);
    strcpy(protein[i].name2,proteome[i]);
    protein[i].copy=AssignNode(copy,proteome,proteosize,proteome[i]); 
  }
  edges=0;
  return edges;
 
}
  

