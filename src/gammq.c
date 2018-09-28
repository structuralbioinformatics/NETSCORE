#include        <stdio.h>

float gammq(a,x)
float a,x;
{
	float gamser,gammcf,gln;
	void gcf(),gser(),nrerror();

	if (x < 0.0 || a <= 0.0) {printf("Arguments in GAMMQ: X= %f and A= %f \n",a,x);nrerror("Invalid arguments in routine GAMMQ");}
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
