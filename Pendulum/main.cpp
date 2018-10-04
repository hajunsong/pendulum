#include <iostream>
#include "pendulum.h"

int main() {
	for (uint i = 0; i <= 1e10*0; i++) {
		Matrix3 a(1), b;
		cout << "a : " << a << endl;

		b(0, 0) = 1;
		b(0, 1) = 2;
		b(0, 2) = 3;
		cout << "b : " << b << endl;

		a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		cout << "a : " << a << endl;

		Matrix3 c = a*a*a*a*a*a*a;
		cout << "c = a*a*a*a*a*a*a : " << c << endl;

		Matrix3 d = a + c + c + c + c + c + c + c;
		cout << "a + c : " << d << endl;

		Matrix3 A;
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		Vector3 B;
		B << 1, 3, 4;
		cout << A * B << endl;

		Matrix6 A2;
		Vector6 B2;
		cout << A2 * B2 << endl;
	}

	system("pause");
	return 0;
}