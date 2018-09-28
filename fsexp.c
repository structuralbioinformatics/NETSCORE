#include <math.h>

void fsexp(x,a,y,dyda,na)
float x,a[],*y,dyda[];
int na;
{
	int i;
	float ex,arg;

	*y=0.0;
	for (i=1;i<=na-1;i+=2)
	{
		arg= a[i+1]*x;
		ex=exp(-arg);
		*y += a[i]*ex;
		dyda[i]=ex;
		dyda[i+1]=-a[i+1]*a[i]*ex;
	}
}
