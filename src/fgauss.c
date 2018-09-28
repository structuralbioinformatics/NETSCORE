#include <math.h>

void fgauss(x,a,y,dyda,na)
float x,a[],*y,dyda[];
int na;
{
	int i;
	float fac,ex,aa,arg;

	*y=0.0;
	for (i=1;i<=na-1;i+=3)
	{
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
                aa=fabs(a[i]);
		fac=aa*ex*2.0*arg;
		*y += aa*ex;
		dyda[i]=a[i]*ex/aa;
		dyda[i+1]=fac/a[i+2];
		dyda[i+2]=fac*arg/a[i+2];
	}
}
