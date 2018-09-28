#include "ppin.h"
static char* GetLine(FILE* inputfile)
/* The function GetLine reads one line from the inputfile, and returns it as a
 * null-terminated string. If inputfile is at EOF, a null pointer is returned.
 * The calling routine should free the char* returned by GetLine.
 */
{ int c;
  int n = 0;
  int size = 1023;
  char* line = malloc((size+1)*sizeof(char));
  while ((c = getc(inputfile))!=EOF && c!='\r' && c!='\n')
  { if (n == size)
    { size *= 2;
      line = realloc(line,(size+1)*sizeof(char));
    }
    line[n] = (char)c;
    n++;
  }
  if (c=='\r')
  { c = getc(inputfile);
    if (c!='\n' && c!=EOF) ungetc(c,inputfile);
  }
  if (n==0 && c==EOF)
  { free(line);
    return 0;
  }
  line[n] = '\0';
  line = realloc(line,(n+1)*sizeof(char));
  return line;
}

static char* tokenize(char* s)
{
  char* p = s;
  while (1)
  {
    if (*p=='\0') return NULL;
    if (*p=='\t')
    {
      *p = '\0';
      return p+1;
    }
    p++;
  }
  /* Never get here */
  return NULL;
}


int ReadClusterFile(x,l,f,g,d)
node   *x;
int     l;
char*   f;
int   **g;
int    *d;
{
 FILE *INP;
 int   ng,i;
 const int ln = strlen(f) + 1;
 const char text[] = "Error: Attempt to read file ";
 const int lm = strlen(text) + 1;
 char* error = malloc((lm + ln)*sizeof(char));
 char* line;
 char* s;
 
 sprintf(error,"%s %s ",text,f);
 ng=0;
 
 INP=fopen(f,"r");
 if (!INP)nrerror(error);
 else printf("Reading Cluster file %s\n", f);

 while((line = GetLine(INP))) 
 {
  if (strlen(line) > 1)  /* Ignore completely empty lines */
  {
    if(!line)
    { const char text[] = "Error: Attempt to read empty file";
      const int m = strlen(text) + 1;
      char* error = malloc(m*sizeof(char));
      strcpy(error,text);
      nrerror(error);
    }
    s=line;
    d[ng]=0;
    while (s)
    { char* token = s;
      s = tokenize(s);
      for (i=0;i<l;i++){
       if(!strcmp(x[i].name1,token)){g[ng][d[ng]]=i;break;}
      } 
      d[ng]++;
      if (d[ng]>MAXXSG)nrerror("Too many elements of a group in Cluster file, increase MAXXSG");
    }
    ng++;
    if (ng>MAXXG)nrerror("Too many groups in Cluster file, increase MAXXG");
  }
 }

 fclose(INP);

 return ng;
}

