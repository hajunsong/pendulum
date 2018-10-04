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
		Matrix mat(rows, cols);
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
			mat(i, col_indx) = data[i*cols + col_indx];
		}
		return mat;
	}
	Matrix row(uint row_indx) {
		Matrix mat(1, cols);
		for (uint i = 0; i < cols; i++) {
			mat(row_indx, i) = data[row_indx*cols + i];
		}
	}
	Matrix insertCol(Matrix mat_in, uint start_indx, uint end_indx) {
		for (uint i = 0; i < rows; i++) {
			for (uint j = start_indx; j <= end_indx; j++) {
				this->operator()(i, j) = mat_in(i, j);
			}
		}
		return *this;
	}
	Matrix insertRow(Matrix mat_in, uint start_indx, uint end_indx) {
		for (uint i = start_indx; i <= end_indx; i++) {
			for (uint j = 0; j <= cols; j++) {
				this->operator()(i, j) = mat_in(i, j);
			}
		}
		return *this;
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

class Pendulum
{
public:
	Pendulum();
	~Pendulum();
};