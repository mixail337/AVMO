#pragma once

class drobi {
public:
	int chisl = 0;
	int znam = 1;

	drobi(){}
	drobi(int ch): chisl(ch){}
	drobi(int ch, int zn): chisl(ch), znam(zn) {
		sokr();
	}

	void sokr();

	static int NOD(int a, int b);
	static int NOK(int a, int b);

};

