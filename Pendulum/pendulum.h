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

template <uint row, uint col>
class Matrix {
public:
	uint rows = 0, cols = 0, len = 0;
	vector<double> data;

	Matrix(double init_value = 0) {
		rows = row;
		cols = col;
		len = rows * cols;
		data.assign(len, init_value);
	}
	~Matrix() {
		data.clear();
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
	template <typename T> Matrix& operator<<(T rhs) {
		if (data.size() == len) {
			data.clear();
		}
		data.push_back(rhs);
		return *this;
	}
	template <typename T> Matrix& operator,(T rhs) {
		data.push_back(rhs);
		return *this;
	}
	Matrix operator*(Matrix rhs) {
		Matrix<row, col> mat;
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
		Matrix<row, col> mat;
		for (uint i = 0; i < len; i++) {
			mat.data[i] += value;
		}
		return mat;
	}
	template <typename T> friend Matrix operator*(T value, Matrix mat) {
		return mat * value;
	}
	Matrix operator+(Matrix rhs) {
		Matrix<row, col> mat;
		for (uint i = 0; i < mat.len; i++) {
			mat.data[i] = data[i] + rhs.data[i];
		}
		return mat;
	}
	Matrix operator-(Matrix rhs) {
		Matrix<row, col> mat;
		for (uint i = 0; i < mat.len; i++) {
			mat.data[i] = data[i] - rhs.data[i];
		}
		return mat;
	}
	Matrix& operator=(Matrix& rhs) {
		this->rows = rhs.rows;
		this->cols = rhs.cols;
		this->len = rhs.len;
		this->data = rhs.data;
	}
	//Matrix operator<<(Matrix rhs) {
	//	Matrix<row, col> mat;

	//	return mat;
	//}
	//Matrix operator,(Matrix rhs) {

	//}
};

typedef Matrix<3, 3> Matrix3;
typedef Matrix<6, 6> Matrix6;
typedef Matrix<3, 1> Vector3;
typedef Matrix<6, 1> Vector6;

inline Vector3 operator*(Matrix3 mat, Vector3 vec) {
	Vector3 vec_result;
	for (uint i = 0; i < mat.rows; i++) {
		double temp = 0;
		for (uint j = 0; j < mat.cols; j++) {
			temp += mat(i, j)*vec(j);
		}
		vec_result(i) = temp;
	}
	return vec_result;
}
inline Vector6 operator*(Matrix6 mat, Vector6 vec) {
	Vector6 vec_result;
	for (uint i = 0; i < mat.rows; i++) {
		double temp = 0;
		for (uint j = 0; j < mat.cols; j++) {
			temp += mat(i, j)*vec(j);
		}
		vec_result(i) = temp;
	}
	return vec_result;
}

//template <typename T, typename U> T operator*(U mat, T vec) {
//	T vec_result;
//	for (uint i = 0; i < mat.rows; i++) {
//		double temp = 0;
//		for (uint j = 0; j < mat.cols; j++) {
//			temp += mat(i, j)*vec(j);
//		}
//		vec_result(i) = temp;
//	}
//	return vec_result;
//}

class Pendulum
{
public:
	Pendulum();
	~Pendulum();
};