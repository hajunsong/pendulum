#include <iostream>
#include "pendulum.h"

int main() {
	uint num_body = 6;
	Pendulum pendulum(num_body);
	pendulum.run();

	system("pause");
	return 0;
}