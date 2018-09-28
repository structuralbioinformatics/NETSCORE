#include <math.h>
#include <stdio.h>

#define EPS 1.0e-6
#define JMAX 14
#define JMAXP JMAX+1
#define K 5

float qromo(func,a,b,choose,p,na)
float a,b,*p;
int   na;
float (*func)();
float (*choose)();	/* ANSI: float choose(float(*)(float),float,float,int); */
{
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];
	void polint(),nrerror();

        if (a==b) return 0.0;
	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j,p,na);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        //printf("QROMO[%d] Err[%d]= %f Integral= %e\n",j,K,dss,ss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	printf ("Too many steps in routing QROMO, skip integral\n");
        return ss;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K
