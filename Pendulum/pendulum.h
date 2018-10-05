#pragma once

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

typedef unsigned int uint;

class Matrix {
public:
	uint rows = 0, cols = 0, len = 0;
	vector<double> data;

	Matrix(uint row, uint col=1) {
		rows = row;
		cols = col;
		len = rows * cols;
		data.assign(len, 0);
	}
	Matrix() {}
	~Matrix() {
		data.clear();
	}
	Matrix t() {
		Matrix mat(cols, rows);
		for (uint i = 0; i < rows; i++) {
			for (uint j = 0; j < cols; j++) {
				mat(j, i) = data[i*cols + j];
			}
		}
		return mat;
	}
	static Matrix zeros(uint row, uint col=1) {
		Matrix mat(row, col);
		for (uint i = 0; i < mat.rows; i++) {
			for (uint j = 0; j < mat.cols; j++) {
				mat(j, i) = 0;
			}
		}
		return mat;
	}
	static Matrix eye(uint m_size) {
		Matrix mat(m_size, m_size);
		for (uint i = 0; i < m_size; i++) {
			for (uint j = 0; j < m_size; j++) {
				mat(j, i) = i == j ? 1 : 0;
			}
		}
		return mat;
	}
	Matrix col(uint col_indx) {
		Matrix mat(rows, 1);
		for (uint i = 0; i < rows; i++) {
			mat(i, 0) = data[i*cols + col_indx];
		}
		return mat;
	}
	Matrix row(uint row_indx) {
		Matrix mat(1, cols);
		for (uint i = 0; i < cols; i++) {
			mat(row_indx, i) = data[row_indx*cols + i];
		}
	}
	void insertCol(Matrix mat_in, uint start_indx, uint end_indx) {
		for (uint i = 0; i < rows; i++) {
			for (uint j = start_indx, j2 = 0; j <= end_indx; j++, j2++) {
				operator()(i, j) = mat_in(i, j2);
			}
		}
	}
	Matrix insertRow(Matrix mat_in, uint start_indx, uint end_indx) {
		for (uint i = start_indx; i <= end_indx; i++) {
			for (uint j = 0; j <= cols; j++) {
				this->operator()(i, j) = mat_in(i, j);
			}
		}
		return *this;
	}
	Matrix tilde() {
		Matrix mat(3, 3);
		mat(0, 0) = 0;
		mat(0, 1) = -operator()(2);
		mat(0, 2) = operator()(1);;
		mat(1, 0) = operator()(2);
		mat(1, 1) = 0;
		mat(1, 2) = -operator()(0);
		mat(2, 0) = -operator()(1);
		mat(2, 1) = operator()(0);
		mat(2, 2) = 0;
		return mat;
	}
	double& operator()(uint row_indx, uint col_indx=0) {
		return data[row_indx*cols + col_indx];
	}
	friend ostream& operator<<(ostream& os, Matrix rhs) {
		cout << endl;
		for (uint i = 0; i < rhs.rows; i++) {
			for (uint j = 0; j < rhs.cols; j++) {
				std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
				std::cout.precision(5);
				cout << rhs(i, j) << "    ";
			}
			cout << endl;
		}
		return os;
	}
	friend Matrix operator-(Matrix rhs) {
		return rhs * (-1);
	}
	Matrix& operator<<(double rhs) {
		if (data.size() == len) {
			data.clear();
		}
		data.push_back(rhs);
		return *this;
	}
	Matrix& operator,(double rhs) {
		data.push_back(rhs);
		return *this;
	}
	Matrix operator*(Matrix rhs) {
		Matrix mat(rows, rhs.cols);
		for (uint i = 0; i < rows; i++) {
			for (uint j = 0; j < rhs.cols; j++) {
				double temp = 0;
				for (uint k = 0; k < cols; k++) {
					temp += data[i*cols + k] * rhs.data[k*rhs.cols + j];
				}
				mat(i, j) = temp;
			}
		}
		return mat;
	}
	template <typename T> Matrix operator*(T value) {
		Matrix mat(rows, cols);
		for (uint i = 0; i < len; i++) {
			mat.data[i] = data[i] * value;
		}
		return mat;
	}
	template <typename T> friend Matrix operator*(T value, Matrix mat) {
		return mat * value;
	}
	template <typename T> Matrix operator/(T value) {
		Matrix mat(rows, cols);
		for (uint i = 0; i < len; i++) {
			mat.data[i] = data[i] / value;
		}
		return mat;
	}
	template <typename T> friend Matrix operator/(T value, Matrix mat) {
		return mat / value;
	}
	Matrix operator+(Matrix rhs) {
		Matrix mat(rows, cols);
		for (uint i = 0; i < mat.len; i++) {
			mat.data[i] = data[i] + rhs.data[i];
		}
		return mat;
	}
	Matrix operator+=(Matrix rhs) {
		for (uint i = 0; i < len; i++) {
			data[i] += rhs.data[i];
		}
		return *this;
	}
	Matrix operator-(Matrix rhs) {
		Matrix mat(rows, cols);
		for (uint i = 0; i < mat.len; i++) {
			mat.data[i] = data[i] - rhs.data[i];
		}
		return mat;
	}
	template <typename T> Matrix operator=(T rhs) {
		this->rows = rhs.rows;
		this->cols = rhs.cols;
		this->len = rhs.len;
		this->data = rhs.data;
		return *this;
	}
	Matrix operator<<(Matrix rhs) {
		Matrix mat(rows+rhs.rows, rhs.cols);
		for (uint i = 0; i < mat.rows; i++) {
			for (uint j = 0; j < mat.cols; j++) {
				mat(i, j) = i < rows ? data[i*cols + j] : rhs(i - rows, j);
			}
		}
		return mat;
	}
	Matrix operator,(Matrix rhs) {
		Matrix mat(rhs.rows, cols+rhs.cols);
		for (uint i = 0; i < mat.rows; i++) {
			for (uint j = 0; j < mat.cols; j++) {
				mat(i, j) = j < cols ? data[i*cols + j] : rhs(i, j - cols);
			}
		}
		return mat;
	}
};
typedef Matrix Vector;

class Integrator {
public:
	Integrator(double h, uint n) {
		step_size = h;

		Y_next = Vector(n);
		AW = Matrix(n, 2);
		AW1 = Matrix(n, 2);
	};

	void absh3(Vector Y, Vector Yp, double t_current) {
		/* ABSH3 : constant step Adams Bashforth 3rd order formulation.
	written by Sung-Soo Kim
	Date: Oct. 19, 1998
	copyright reserved by Sung-Soo Kim

	input variables
	t_current: current time
	Y : current state
	Yp : current derivative of state
	step_size: integration step_size

	output variables
	Y_next : state at next time step
	t_next : Next time

	STARTER: upto 2h, i.e., derivatives are stored for the initial  time steps at 0, h, 2h, to form
	3rd order Adams Bashforth formula */

		switch (intcount) {
			case 1:
				// Forward Euler method with 0.25 step_size for initial step
				// use derivative information at 0 step
				// y=y+step_size*Yp/4.0;
				Y_next = Y + step_size * Yp / 4.0;
				// w(:,2) = Yp;
				AW.insertCol(Yp, 1, 1);
				// w1(:,2) = Yp;
				AW1.insertCol(Yp, 1, 1);
				t_next = t_current + step_size / 4.0;
				break;
			case 2:
				// Adams Bashforth 2nd order method with 0.25 step_size for 2nd step
				// use derivative information at 0, h/4
				// y = y + step_size_h * ( 3.0*Yp - w1(:,2))/8.0;
				Y_next = Y + step_size * (3.0 * Yp - AW1.col(1)) / 8.0;
				// w1(:,1) = Yp;
				AW1.insertCol(Yp, 0, 0);
				t_next = t_current + step_size / 4.0;
				break;
			case 3:
				// Adams Bashforth 3rd order method with 0.25 step_size for 3rd step
				// use derivative information at 0, h/4, h/2
				// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/48.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW1.col(0) + 5.0*AW1.col(1)) / 48.0;
				// w1(:,2) = w1(:,1);
				AW1.insertCol(AW1.col(0), 1, 1);
				// w1(:,1) = Yp;
				AW1.insertCol(Yp, 0, 0);
				t_next = t_current + step_size / 4.0;
				break;
			case 4:
				// Adams Bashforth 3rd order method with 0.25 step_size for 4th step
				// use derivative information at h/4, h/2, 3h/4
				// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/48.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW1.col(0) + 5.0*AW1.col(1)) / 48.0;
				// w1(:,2) = w(:,2);
				AW1.insertCol(AW.col(1), 1, 1);
				t_next = t_current + step_size / 4.0;
				break;
			case 5:
				// Adams Bashforth 3rd order method with 0.5 step_size for 5th step
				// use derivative information at 0, h/2, h
				// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/24.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW1.col(0) + 5.0*AW1.col(1)) / 24.0;
				// w(:,1) = Yp;
				AW.insertCol(Yp, 0, 0);
				// w1(:,2) = w1(:,1);
				AW1.insertCol(AW1.col(0), 1, 1);
				// w1(:,1) = Yp;
				AW1.insertCol(Yp, 0, 0);
				t_next = t_current + step_size / 2.0;
				break;
			case 6:
				// Adams Bashforth 3rd order method with 0.5 step_size for 6th step
				// use derivative information at h/2, h,  3h/2
				// y = y + step_size * ( 23.0*Yp - 16.0*w1(:,1) + 5.0*w1(:,2))/24.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW1.col(0) + 5.0*AW1.col(1)) / 24.0;
				// w1(:,2) = w1(:,1);
				AW1.insertCol(AW1.col(0), 1, 1);
				// w1(:,1) = Yp;
				AW1.insertCol(Yp, 0, 0);
				t_next = t_current + step_size / 2.0;
				break;
			case 7:
				// Adams Bashforth 3rd order method with step_size for 7th step
				// use derivative information at 0,  h,  2h
				// y = y + step_size * ( 23.0*Yp - 16.0*w(:,1) + 5.0*w(:,2))/12.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW.col(0) + 5.0*AW.col(1)) / 12.0;
				// w(:,2) = w(:,1);
				AW.insertCol(AW.col(0), 1, 1);
				// w(:,1) = Yp;
				AW.insertCol(Yp, 0, 0);
				t_next = t_current + step_size;
				break;
			default:
				// Adams Bashforth 3rd order method with step_size for more than 8th step
				// use derivative information t_current-2h, t_current-h, t_current
				// y = y + step_size * ( 23.0*Yp - 16.0*w(:,1) + 5.0*w(:,2))/12.0;
				Y_next = Y + step_size * (23.0*Yp - 16.0*AW.col(0) + 5.0*AW.col(1)) / 12.0;
				// w(:,2) = w(:,1);
				AW.insertCol(AW.col(0), 1, 1);
				// w(:,1) = Yp;
				AW.insertCol(Yp, 0, 0);
				t_next = t_current + step_size;
				break;
		}
		intcount++;
	}
	// void rectangle();
	// void trapezoidal();

	Vector Y_next;
	double t_next = 0, step_size;
	int intcount = 1;
private:
	Matrix AW, AW1;
};

class Body
{
public:
	Body() {
		u_vec << 0, 0, 1;
	};
	// base body information
	Matrix A0 = Matrix(3,3), C01 = Matrix(3,3);
	Vector s01p = Vector(3);
	// body initial data
	double qi, qi_dot, m;
	Vector ri = Vector(3), ri_dot = Vector(3), wi = Vector(3), rhoip = Vector(3), sijp = Vector(3);
	Matrix Jip = Matrix(3, 3), Cii = Matrix(3, 3), Cij = Matrix(3, 3);
	// Orientation
	Matrix Aijpp = Matrix(3, 3), Ai = Matrix(3, 3);
	Vector Hi = Vector(3), u_vec = Vector(3);
	// Position
	Vector sij = Vector(3), rhoi = Vector(3), ric = Vector(3);
	Matrix rit = Matrix(3, 3);
	// Velocity State
	Vector Bi = Vector(6), Yih = Vector(6);
	// Cartesian velocity state
	Matrix Ti = Matrix(6, 6), wit = Matrix(3, 3);
	Vector Yib = Vector(6), ric_dot = Vector(3);
	// Mass & Force
	Matrix Jic = Matrix(3, 3), rict = Matrix(3, 3), rict_dot = Matrix(3, 3), Mih = Matrix(6, 6);
	Vector Fic = Vector(3), Tic = Vector(3), Qih = Vector(6);
	// Velocity Coupling
	Matrix rit_dot = Matrix(3, 3);
	Vector dHi = Vector(3), Di = Vector(6);
	// System EQM
	Matrix Ki = Matrix(6, 6);
	Vector Li = Vector(6);
	// Acceleration
	double qi_ddot;
};

class Pendulum
{
public:
	Pendulum(uint numbody);
	~Pendulum();
	void run();

private:
	uint num_body;
	Body *body;
	Integrator *integr;
	// system variable
	double start_time, end_time, h, g, t_current;
	// state vector
	Vector Y, Yp;
	// file
	char file_name[256];
	FILE *fp;

	void analysis();
		Matrix ludcmp(Matrix a, uint* indx);
		Vector lubksb(Matrix fac, uint* indx, Vector b);
	void save_data();
};