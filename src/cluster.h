/******************************************************************************/
/* The C Clustering Library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon 'AT' gsc.riken.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

#ifdef WINDOWS
#  include <windows.h>
#endif

#define CLUSTERVERSION "1.50"

/* Chapter 2 */
float clusterdistance (int nrows, int ncolumns, float** data, int** mask,
  float weight[], int n1, int n2, int index1[], int index2[], char dist,
  char method, int transpose);
float** distancematrix (int ngenes, int ndata, float** data,
  int** mask, float* weight, char dist, int transpose);

/* Chapter 3 */
int getclustercentroids(int nclusters, int nrows, int ncolumns,
  float** data, int** mask, int clusterid[], float** cdata, int** cmask,
  int transpose, char method);
void getclustermedoids(int nclusters, int nelements, float** distance,
  int clusterid[], int centroids[], float errors[]);
void kcluster (int nclusters, int ngenes, int ndata, float** data,
  int** mask, float weight[], int transpose, int npass, char method, char dist,
  int clusterid[], float* error, int* ifound);
void kmedoids (int nclusters, int nelements, float** distance,
  int npass, int clusterid[], float* error, int* ifound);

/* Chapter 4 */
typedef struct {int left; int right; float distance;} tNode;
/*
 * A Node struct describes a single node in a tree created by hierarchical
 * clustering. The tree can be represented by an array of n Node structs,
 * where n is the number of elements minus one. The integers left and right
 * in each Node struct refer to the two elements or subnodes that are joined
 * in this node. The original elements are numbered 0..nelements-1, and the
 * nodes -1..-(nelements-1). For each node, distance contains the distance
 * between the two subnodes that were joined.
 */

tNode* treecluster (int nrows, int ncolumns, float** data, int** mask,
  float weight[], int transpose, char dist, char method, float** distmatrix);
void cuttree (int nelements, tNode* tree, int nclusters, int clusterid[]);

/* Chapter 5 */
void somcluster (int nrows, int ncolumns, float** data, int** mask,
  const float weight[], int transpose, int nxnodes, int nynodes,
  float inittau, int niter, char dist, float*** celldata,
  int clusterid[][2]);

/* Chapter 6 */
int pca(int m, int n, float** u, float** v, float* w);

/* Utility routines, currently undocumented */
void sort(int n, const float data[], int index[]);
float mean(int n, float x[]);
float median (int n, float x[]);

float* calculate_weights(int nrows, int ncolumns, float** data, int** mask,
  float weights[], int transpose, char dist, float cutoff, float exponent);
