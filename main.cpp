/*
 * main.cpp
 *
 *  Created on: 20 janv. 2014
 *      Author: ludotravail
 */


#include "lazy_matrix.h"
#include <iostream>
#include <memory>

using namespace lazy;

int main(void){
	lazy::matrix<double> m1;
	lazy::matrix<double> m2;

	int n = 3, m = 3;

	m1.init(n,m);
	m2.init(m,n);

	auto m3 = m1 * m2;
	auto m4 = m1 + m2;
	auto m5 = m3 + m4;

	for(int i = 0; i < m1.row(); ++i){
		for (int j = 0; j < m1.col(); ++j){
			m1.setValue(i,j,double(i+j));
		}
	}

	for(int i = 0; i < m2.row(); ++i){
		for (int j = 0; j < m2.col(); ++j){
			m2.setValue(i,j,2.0*double(i+j));
		}
	}

	std::cout << "m1 is " << m1.row() << " x " << m1.col() << " :" << std::endl;
	disp(m1);
	std::cout << "m2 is " << m2.row() << " x " << m2.col() << " :" << std::endl;
	disp(m2);
	std::cout << "m3 is " << m3.row() << " x " << m3.col() << " :" << std::endl;
	disp(m3);
	std::cout << "m4 is " << m4.row() << " x " << m4.col() << " :" << std::endl;
	disp(m4);
	std::cout << "m5 is " << m5.row() << " x " << m5.col() << " :" << std::endl;
	disp(m5);

	lazy::matrix<double> m6; m6 = m5;
	std::cout << "m6 is " << m6.row() << " x " << m6.col() << " :" << std::endl;
	disp(m6);

}



