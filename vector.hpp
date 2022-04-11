#pragma once
#include <vector>
#include "/scripts/code/dcheck.hpp"
#include <stdexcept>

template<class T>
class Vector : public std::vector<T>
{
	public:
		using std::vector<T>::vector;
#ifndef NDEBUG
	T& operator[](size_t n) /*override*/ { 
		DCHECK_LT(n, this->size());
		return std::vector<T>::operator[](n); 
	}
	const T& operator[](size_t n) const /*override*/ { 
		DCHECK_LT(n, this->size());
		return std::vector<T>::operator[](n);
	}
#endif
};
