#include <iostream>
#include "DBL.h"
using namespace jys;
using namespace std;

dbl<4> c[4], r;

int main (int argc, char ** argv)
{
	c[0]=c[1]=c[2]= -1.0;
	c[3]=1.0; // x^3 - x^2 - x - 1 = 0

	if (argc>1) cout << hex; // invoke as "./sample h" to see the elements in hex

	r = polyroot<4> (c, 4, 1.0, 2.0);
	r.printElm (cout) << endl;
	cout << r << endl;
	cout << sqr (r) << endl;
	cout << r*r*r - r*r - r - 1.0 << endl;
	cout << (r=poly (c, 4, r)) << endl;
	r.printElm (cout) << endl;

	cout << "new value\n";
	r="3.99999999999999999999999999999999999999999999999e-124";
	r.printElm (cout) << endl;
	cout << "3.99999999999999999999999999999999999999999999999e-124" << endl;
	cout << r << endl;
	cout << sqrt(r) << endl;
}
