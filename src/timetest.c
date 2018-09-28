#include <time.h>
#include        <stdio.h>
#include        <ctype.h>
#include        <string.h>
#include        <sys/types.h>
#include        <sys/wait.h>
#include        <malloc.h>
#include        <openssl/des.h>
#include        <stdlib.h>
#include        <unistd.h>
#include        <errno.h>

main()
{
	int i;
clock_t launch = clock();
//do work
printf("Test \n");
for (i=0;i<10000000;i++)
	if (i<1000000) printf("Test \n");

clock_t done = clock();
float  diff = (done - launch) / CLOCKS_PER_SEC;
printf("Test Time %f %d %d %d \n",diff,launch,done,CLOCKS_PER_SEC);

}
