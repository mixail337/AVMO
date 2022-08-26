#include "drobs.h"

int drobi::NOD(int a, int b) {
	int t;
	if (a < b) {
		t = a;
		a = b;
		b = t;
	}
	while (b != 0) {
		t = b;
		b = a % b;
		a = t;
	}
	if (!a) return 1;
	return a;
}

int drobi::NOK(int a, int b) {
	return (a * b) / NOD(a, b);
}

void drobi::sokr() {
	int nod = NOD(chisl, znam);
	chisl /= nod;
	znam /= nod;
	if (znam < 0) {
		znam *= -1;
		chisl *= -1;
	}
}
