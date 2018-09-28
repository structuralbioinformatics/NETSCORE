#include <math.h>

float ffgauss(x,a,na)
float x,a[];
int na;
{
	int i;
	float y,ex,arg;

	y=0.0;
	for (i=1;i<=na-1;i+=3)
	{
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		y += fabs(a[i])*ex;
	}

        return y;

}
