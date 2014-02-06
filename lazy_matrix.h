/*
 * lazy_matrix.h
 *
 *  Created on: 20 janv. 2014
 *      Author: Ludovic Henno
 */

#ifndef LAZY_MATRIX_H_
#define LAZY_MATRIX_H_

#include <memory>
#include <boost/shared_array.hpp>
#include <iostream>
#include <unordered_map>


namespace lazy {

template<typename T>
class expr_mat {
public:
	virtual ~expr_mat() {
	}
	;
	virtual T operator()(const int i, const int j) const = 0;
	virtual int row() const = 0;
	virtual int col() const = 0;
};

template<typename T>
class matrix: public expr_mat<T> {

private:
	int _row, _col;
	boost::shared_array<T> _body;
public:
	matrix() {
		init(0,0);
	}
	;

	explicit matrix(matrix<T> &o) :
		_col(o.col()), _row(o.row()){
		_body = o._body;
	}
	;

	void init(const int &new_row, const int &new_col) {
		_col = new_col;
		_row = new_row;
		_body = boost::shared_array<T>(new T[_col * _row]);
	}

	virtual ~matrix() {
	}
	;

	T getValue(const int &i, const int &j) const {
		return _body[i * _row + j];
	}

	void setValue(const int &i, const int &j, const T &value) {
		if (i < _row && j < _col)
			_body.get()[i * _col + j] = value;
	}

	virtual T operator()(const int i, const int j) const {
		return getValue(i, j);
	}

	bool operator==(const matrix<T> &o) const {
		if (_row != o._row || _col != o._col)
			return false;
		int i;
		for (i = 0; i < _row * _col; ++i) {
			if (_body[i] != o._body[i])
				return false;
		}
		return true;
	}

	template<typename Mat>
	const matrix<T> & operator=(const Mat &o) {
		_row = o.row();
		_col = o.col();

		_body.reset(new T[_col * _row]);

		for (int i = 0; i < _row; ++i) {
			for (int j = 0; j < _col; ++j) {
				_body[i * _col + j] = o(i,j);
			}
		}

		return *this;
	}



	virtual int col() const {
		return _col;
	}

	virtual int row() const {
		return _row;
	}

};

template<typename T,typename M>
class step_matrix{
private:
	int _row, _col;
	std::unordered_map<int,T> _body;

public:
	step_matrix(){};
	template <typename M>
	step_matrix(const M &m){
		_row = m.row();
		_col = m.col();
	}
	virtual ~step_matrix(){};

	/*void init(const int &row, const int &col){
		_row = row;
		_col = col;
	}*/
};

template<typename Left, typename Right>
class added_matrix {
private:
	const Left &_lhs;
	const Right &_rhs;

public:
	added_matrix() {
	}
	;
	added_matrix(const Left &lhs, const Right &rhs) :
			_lhs(lhs), _rhs(rhs) {
	}

	added_matrix(const added_matrix<Left, Right> &o) :
			_lhs(o._lhs), _rhs(o._rhs) {
	}

	virtual ~added_matrix() {

	}

	auto getValue(const int i,
			const int j) const -> decltype(_lhs(i,j) + _rhs(i,j)) {

		return _lhs(i, j) + _rhs(i, j);
	}

	auto operator()(const int i,
			const int j) const -> decltype(_lhs(i,j) + _rhs(i,j)) {
		return getValue(i, j);
	}

	int col() const {
		return _lhs.col();
	}

	int row() const {
		return _lhs.row();
	}

};

template<typename Left, typename Right>
class subed_matrix {
private:
	const Left &_lhs;
	const Right &_rhs;

public:
	subed_matrix() {
	}
	;
	subed_matrix(const Left &lhs, const Right &rhs) :
			_lhs(lhs), _rhs(rhs) {
	}

	subed_matrix(const subed_matrix<Left, Right> &o) :
			_lhs(o._lhs), _rhs(o._rhs) {
	}

	virtual ~subed_matrix() {

	}

	auto getValue(const int i,
			const int j) const -> decltype(_lhs(i,j) - _rhs(i,j)) {
		return _lhs(i, j) - _rhs(i, j);
	}

	auto operator()(const int i,
			const int j) const -> decltype(_lhs(i,j) - _rhs(i,j)) {
		return getValue(i, j);
	}

	int col() const {
		return _lhs.col();
	}

	int row() const {
		return _lhs.row();
	}

};

template<typename Left, typename Right>
class mult_matrix {
private:
	const Left &_lhs;
	const Right &_rhs;

public:
	mult_matrix() {
	}
	;
	mult_matrix(const Left &lhs, const Right &rhs) :
			_lhs(lhs), _rhs(rhs) {
	}

	mult_matrix(const mult_matrix<Left, Right> &o) :
			_lhs(o._lhs), _rhs(o._rhs) {
	}

	virtual ~mult_matrix() {

	}

	auto getValue(const int i,
			const int j) const -> decltype(_lhs(i,j) * _rhs(i,j)) {
		decltype(_lhs(i,j) * _rhs(i,j)) res = 0;

		for (int k = 0; k < _lhs.col(); ++k) {
			res += _lhs(i, k) * _rhs(k, j);
		}

		return res;
	}

	auto operator()(const int i,
			const int j) const -> decltype(_lhs(i,j) - _rhs(i,j)) {
		return getValue(i, j);
	}

	int col() const {
		return _rhs.col();
	}

	int row() const {
		return _lhs.row();
	}

};

template<typename Mat>
class trans_matrix {
private:
	const Mat &_m;

public:
	trans_matrix() {
	}
	;
	trans_matrix(const Mat &m) :
			_m(m) {
	}
	trans_matrix(const trans_matrix<Mat> &o) :
			_m(o._m) {
	}

	virtual ~trans_matrix() {

	}

	auto getValue(const int i, const int j) const -> decltype(_m(j,i)) {

		return _m(j, i);
	}

	auto operator()(const int i, const int j) const -> decltype(_m(j,i)) {
		return getValue(i, j);
	}

	int col() const {
		return _m.row();
	}

	int row() const {
		return _m.col();
	}

};

}
;

template<typename Left, typename Right>
const lazy::added_matrix<Left, Right> operator+(const Left &lhs,
		const Right &rhs) {
	return lazy::added_matrix<Left, Right>(lhs, rhs);

}

template<typename Left, typename Right>
const lazy::subed_matrix<Left, Right> operator-(const Left &lhs,
		const Right &rhs) {
	return lazy::subed_matrix<Left, Right>(lhs, rhs);
}

template<typename Left, typename Right>
const lazy::mult_matrix<Left, Right> operator*(const Left &lhs,
		const Right &rhs) {
	return lazy::mult_matrix<Left, Right>(lhs, rhs);
}

template<typename Mat>
const lazy::trans_matrix<Mat> trans(const Mat &m) {
	return lazy::trans_matrix<Mat>(m);
}


template<typename Mat>
auto tr(Mat &m) ->decltype(m(0,0)) {
	decltype(m(0,0)) res = 0;
	for (int i = 0; i < m.row(); ++i) {
		for (int j = 0; j < m.col(); ++j) {
			res += m(i, j);
		}
	}
	return res;
}

template<typename Mat>
void disp(Mat &m, std::ostream &out = std::cout){
	for(int i=0; i< m.row(); ++i){
		for(int j = 0; j< m.col(); ++j){
			std::cout << m(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

#endif /* LAZY_MATRIX_H_ */
