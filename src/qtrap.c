#include <math.h>
#include <stdio.h>

#define EPS 1.0e-6
#define JMAX 20

float qtrap(func,a,b,p,na)
float a,b,*p;
int   na;
float (*func)();
{
	int j;
	float s,olds,trapzd();
	void nrerror();

        if (a==b) return 0.0;
	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j,p,na);
                //printf("QTRAP[%d] Pre-Integral= %e Integral= %e\n",j,olds,s);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	printf("Too many steps in routine QTRAP\n");
        return s;
}

#undef EPS
#undef JMAX
