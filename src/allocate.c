#include  "ppin.h"


 void nrerror(error_text)
 char error_text[];
 {
 void exit();
 fprintf(stderr,"Numerical Recipes run-time error...\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"...now exiting to system...\n");
 exit(1);
 }



void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((int *) (m[i]+ncl));
   free((int *) (m+nrl));
}


int  **imatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 int  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(int **)  malloc((unsigned) (nrh-nrl+1)*sizeof(int *));
 if (!m) nrerror("allocation failure 1 in imatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
  if (!m[i]) nrerror("allocation failure 2 in imatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

float  **matrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 float  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(float **)  malloc((unsigned) (nrh-nrl+1)*sizeof(float *));
 if (!m) nrerror("allocation failure 1 in matrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
  if (!m[i]) nrerror("allocation failure 2 in matrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((float *) (m[i]+ncl));
   free((float *) (m+nrl));
}

char   **cmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 char   **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(char  **)  malloc((unsigned) (nrh-nrl+1)*sizeof(char  *));
 if (!m) nrerror("allocation failure 1 in cmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(char  *) malloc((unsigned) (nch-ncl+1)*sizeof(char ));
  if (!m[i]) nrerror("allocation failure 2 in cmatrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_cmatrix(m,nrl,nrh,ncl,nch)
char  **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((char  *) (m[i]+ncl));
   free((char  *) (m+nrl));
}



node   *avector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   node  *v;
   void nrerror();
   v=(node *)malloc((unsigned) (nh-nl+1)*sizeof(node));
   if (!v) nrerror("allocation failure in avector");
   return v-nl;
 }

void free_avector(v,nl,nh)
node  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((node*) (v+nl));
}




express   *xvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   express  *v;
   void nrerror();
   v=(express *)malloc((unsigned) (nh-nl+1)*sizeof(express));
   if (!v) nrerror("allocation failure in xvector");
   return v-nl;
 }

void free_xvector(v,nl,nh)
express  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((express*) (v+nl));
}





edge   *evector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   edge  *v;
   void nrerror();
   v=(edge *)malloc((unsigned) (nh-nl+1)*sizeof(edge));
   if (!v) nrerror("allocation failure in evector");
   return v-nl;
 }

void free_evector(v,nl,nh)
edge  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((edge*) (v+nl));
}


int   *ivector(nl,nh)
 int nl, nh;
 /* Allocates a sequence vector with range [nl...nh].  */
 {
   int  *v;
   void nrerror();
   v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
   if (!v) nrerror("allocation failure in ivector");
   return v-nl;
 }

void free_ivector(v,nl,nh)
int  *v;
int nl, nh;
/* Frees an sequence vector allocated by vector() */
{
 free((int*) (v+nl));
}


char  *cvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   char  *v;
   void nrerror();
   v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
   if (!v) nrerror("allocation failure in cvector");
   return v-nl;
 }

void free_cvector(v,nl,nh)
char  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((char*) (v+nl));
}

float  *fvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   float  *v;
   void nrerror();
   v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
   if (!v) nrerror("allocation failure in fvector");
   return v-nl;
 }

void free_fvector(v,nl,nh)
float  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((float*) (v+nl));
}


gradient  *gvector(nl,nh)
 int nl, nh;
 /* Allocates an alignment vector with range [nl...nh].  */
 {
   gradient  *v;
   void nrerror();
   v=(gradient *)malloc((unsigned) (nh-nl+1)*sizeof(gradient));
   if (!v) nrerror("allocation failure in gvector");
   return v-nl;
 }

void free_gvector(v,nl,nh)
gradient  *v;
int nl, nh;
/* Frees an alignment vector allocated by vector() */
{
 free((gradient*) (v+nl));
}


float  ***kform(nrl,nrh,ncl,nch,ndl,ndh)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch][ndl..ndh]  */
int nrl, nrh, ncl, nch, ndl, ndh;
{
 int i,j;
 float  ***m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(float ***)  malloc((unsigned) (nrh-nrl+1)*sizeof(float **));
 if (!m) nrerror("allocation failure 1 in kform");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(float **) malloc((unsigned) (nch-ncl+1)*sizeof(float *));
  if (!m[i]) nrerror("allocation failure 2 in kform");
  m[i] -= ncl;
  for (j=ndl;j<=ndh;j++){
    m[i][j]=(float *) malloc((unsigned) (ndh-ndl+1)*sizeof(float ));
    if (!m[i][j]) nrerror("allocation failure 3 in kform");
    m[i][j] -= ndl;
  }
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_kform(m,nrl,nrh,ncl,nch,ndl,ndh)
float ***m;
int nrl, nrh, ncl, nch, ndl, ndh;
/* Frees a matrix allocated with matrix */
{
   int i,j;
   for (i=nrh;i>=nrl;i--) for (j=ndh;j>=ndl;j--) free((float **) (m[i][j]+ndl));
   for (i=nrh;i>=nrl;i--) free((float *) (m[i]+ncl));
   free((float *) (m+nrl));
}


gradient  **gmatrix(nrl,nrh,ncl,nch)
 /* Allocate a Structure=SIMIL  matrix with range [nrl..nrh][ncl..nch]  */
int nrl, nrh, ncl, nch;
{
 int i;
 gradient  **m;
 void nrerror();

 /* Allocate pointers to rows. */
 m=(gradient **)  malloc((unsigned) (nrh-nrl+1)*sizeof(gradient *));
 if (!m) nrerror("allocation failure 1 in gmatrix");
 m -= nrl;
 /* Allocate rows ans set pointers to them . */
 for(i=nrl;i<=nrh;i++)
 {
  m[i]=(gradient *) malloc((unsigned) (nch-ncl+1)*sizeof(gradient));
  if (!m[i]) nrerror("allocation failure 2 in matrix");
  m[i] -= ncl;
 }
 /* Return pointer to array of pointers to rows. */
 return m;
}

void free_gmatrix(m,nrl,nrh,ncl,nch)
gradient **m;
int nrl, nrh, ncl, nch;
/* Frees a matrix allocated with matrix */
{
   int i;
   for(i=nrh;i>=nrl;i--) free((gradient *) (m[i]+ncl));
   free((gradient *) (m+nrl));
}



