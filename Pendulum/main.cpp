#include <iostream>
#include "pendulum.h"

int main() {
	for (uint i = 0; i <= 1e10*0; i++) {
		Matrix a(3,3), b(3,3);
		cout << "a : " << a << endl;

		b(0, 0) = 1;
		b(0, 1) = 2;
		b(0, 2) = 3;
		cout << "b : " << b << endl;

		a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		cout << "a : " << a << endl;

		Matrix c = a*a*a*a*a*a*a;
		cout << "c = a*a*a*a*a*a*a : " << c << endl;

		Matrix d = a + c + c + c + c + c + c + c;
		cout << "a + c : " << d << endl;

		Matrix A(3,3);
		A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		Vector B(3);
		B << 1, 3, 4;
		cout << "A*B : " << A * B << endl;

		Matrix A2(6,6);
		Vector B2(6);
		cout << "A2*B2 : " << A2 * B2 << endl;

		Matrix C(3, 6);
		C = (A, A) << (A, A);
		cout << "C : " << C << endl;

		Matrix aa(3,3);
		aa << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		Vector bb(3);
		bb << 11, 22, 33;
	}

	Pendulum pen;

	system("pause");
	return 0;
}