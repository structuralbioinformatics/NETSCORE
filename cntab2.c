#include <math.h>

#define TINY 1.0e-30

void cntab2(nn,ni,nj,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
int ni,nj,**nn;
float *h,*hx,*hy,*hygx,*hxgy,*uygx,*uxgy,*uxy;
{
int i,j,k;
float sum,p,*sumi,*sumj;
float *fvector();
void free_fvector();
     
        sum=0.0;
	sumi=fvector(0,ni+1);
	sumj=fvector(0,nj+1);
	for (i=1;i<=ni;i++) {
		sumi[i]=0.0;
		for (j=1;j<=nj;j++) {
			sumi[i] += nn[i][j];
			sum += nn[i][j];
                        //printf("SUMI[%d]=%f SUM=%f\n",i,sumi[i],sum);
		}
	}
	for (j=1;j<=nj;j++) {
		sumj[j]=0.0;
		for (i=1;i<=ni;i++)
			sumj[j] += nn[i][j];
	}
        //for (j=1;j<=nj;j++)printf("SUMJ[%d]=%f\n",j,sumj[j]);
        //for (j=1;j<=ni;j++)printf("SUMI[%d]=%f\n",j,sumi[j]);
	*hx=0.0;
	for (k=1; k<=ni; k++){
                //printf("test %d NI %d S %f\n",k,ni,sumi[k]);
		if (sumi[k]!=0) {
			p=sumi[k]/sum;
                        //printf("Pi=%f I=%d\n",p,k);
			*hx -= p*log(p);
                        //printf("HX=%f\n",*hx);
		}
        }
	*hy=0.0;
	for (j=1;j<=nj;j++)
		if (sumj[j]) {
			p=sumj[j]/sum;
			*hy -= p*log(p);
		}
	*h=0.0;
	for (i=1;i<=ni;i++)
		for (j=1;j<=nj;j++)
			if (nn[i][j]) {
				p=nn[i][j]/sum;
				*h -= p*log(p);
			}
	*hygx=(*h)-(*hx);
	*hxgy=(*h)-(*hy);
	*uygx=(*hy-*hygx)/(*hy+TINY);
	*uxgy=(*hx-*hxgy)/(*hx+TINY);
	*uxy=2.0*(*hx+*hy-*h)/(*hx+*hy+TINY);
	free_fvector(sumj,0,nj+1);
	free_fvector(sumi,0,ni+1);
}

#undef TINY
