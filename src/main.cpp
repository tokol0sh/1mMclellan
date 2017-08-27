#include <cstdio>

#include <stdio.h>      /* puts, printf */
#include <time.h>  


struct tPoint {
	double MA;
	double ME;
};

struct telescope {
	double HA_observed;
	double DEC_observed;
	double HA_actual;
	double DEC_actual;
	double LST;
	tPoint Model;
};

int main(int argc, char** argv) {

	time_t rawtime;
	struct tm * ptm;

	time(&rawtime);

	ptm = gmtime(&rawtime);

	puts("Current time around the World:");
	while (1) {
		time(&rawtime);
		ptm = gmtime(&rawtime);
		printf("current time :  %2d:%02d:%02d\n", (ptm->tm_hour + 12) % 24, ptm->tm_min, ptm->tm_sec);
		//printf("Reykjavik (Iceland) : %2d:%02d\n", (ptm->tm_hour + 13) % 24, ptm->tm_min);
		//printf("Beijing (China) :     %2d:%02d\n", (ptm->tm_hour + 2) % 24, ptm->tm_min);
	}
	return 0;
}