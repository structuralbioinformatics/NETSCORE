#include <math.h>

#define TINY 1.0e-30

void cntab1(nn,ni,nj,chisq,df,prob,cramrv,ccc)
int **nn,ni,nj;
float *chisq,*df,*prob,*cramrv,*ccc;
{
	int nnj,nni,j,i,minij;
	float sum=0.0,expctd,*sumi,*sumj,temp,gammq(),*fvector();
	void free_fvector();

	sumi=fvector(1,ni+1);
	sumj=fvector(1,nj+1);
	nni=ni;
	nnj=nj;
	for (i=1;i<=ni;i++) {
		sumi[i]=0.0;
		for (j=1;j<=nj;j++) {
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
		if (sumi[i] == 0.0) --nni;
	}
	for (j=1;j<=nj;j++) {
		sumj[j]=0.0;
		for (i=1;i<=ni;i++) sumj[j] += nn[i][j];
		if (sumj[j] == 0.0) --nnj;
	}
	*df=nni*nnj-nni-nnj+1;
	*chisq=0.0;
	for (i=1;i<=ni;i++) {
		for (j=1;j<=nj;j++) {
			expctd=sumj[j]*sumi[i]/sum;
			temp=nn[i][j]-expctd;
			*chisq += temp*temp/(expctd+TINY);
		}
	}
	*prob=gammq(0.5*(*df),0.5*(*chisq));
	minij = nni < nnj ? nni-1 : nnj-1;
	*cramrv=sqrt(*chisq/(sum*minij));
	*ccc=sqrt(*chisq/(*chisq+sum));
	free_fvector(sumj,1,nj+1);
	free_fvector(sumi,1,ni+1);
}

#undef TINY