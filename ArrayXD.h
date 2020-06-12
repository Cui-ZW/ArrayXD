// ArrayXD.h
#ifndef ARRAYXD_H_
#define ARRAYXD_H_
#include <iostream>
#include <complex>
#include <math.h>
using namespace std;
typedef unsigned int uint;

template <typename T>
class ArrayXD
{
public:
	uint N0,N1,N2,N3;
	uint Dim;
	uint Size;
	T *m;
public:
	ArrayXD(uint n0,uint n1=1,uint n2=1,uint n3=1);  
	ArrayXD();									    
	ArrayXD(const ArrayXD&);
	~ArrayXD();
	void Show()const;
//  Get
	inline uint Get_Dim(){return Dim;}
	inline T    Get_EndValue(){return m[Size-1];}
	inline T    Get_StartValue(){return m[0];}



// operator
	T & operator()(uint ind0, uint ind1=0, uint ind2=0,uint int3=0);
	ArrayXD operator()(const ArrayXD<int> & ind0, const ArrayXD<int> &ind1=1, const ArrayXD<int> &ind2=1,const ArrayXD<int> &ind3=1);
	
	ArrayXD& operator=(const ArrayXD&);
	ArrayXD& operator=(const T &val);

	ArrayXD& operator+=(const ArrayXD&);
    	ArrayXD& operator-=(const ArrayXD&);
    	ArrayXD& operator*=(const ArrayXD&);
	ArrayXD& operator/=(const ArrayXD&);
	ArrayXD& operator+=(const T &val);
	ArrayXD& operator-=(const T &val);
   	ArrayXD& operator*=(const T &val);
	ArrayXD& operator/=(const T &val);

	ArrayXD operator+(const ArrayXD&);
	ArrayXD operator-(const ArrayXD&);
	ArrayXD operator*(const ArrayXD&);
	ArrayXD operator/(const ArrayXD&);
    	ArrayXD operator+(const T &val);
	ArrayXD operator-(const T &val);
	ArrayXD operator*(const T &val);
	ArrayXD operator/(const T &val);

	template <typename T> friend ArrayXD<T>  operator+(const T &val,const ArrayXD<T>&);
	template <typename T> friend ArrayXD<T>  operator-(const T &val,const ArrayXD<T>&);
	template <typename T> friend ArrayXD<T>  operator*(const T &val,const ArrayXD<T>&);
	template <typename T> friend ArrayXD<T>  operator/(const T &val,const ArrayXD<T>&);


};

template <typename T> ArrayXD<T> operator+(const T &val,const ArrayXD<T>&);		
template <typename T> ArrayXD<T> operator-(const T &val,const ArrayXD<T>&);	
template <typename T> ArrayXD<T> operator/(const T &val,const ArrayXD<T>&);


// double/int to complex<double>
template<typename T> inline complex<double> dcplx(T real, T imag=0. ){
	complex<double> temp((double)real,(double)imag);
	return temp;
};



ArrayXD<int> range(const int &ind0, const int &ind1,const int &dind=1);
ArrayXD<int> range(const int &ind0);
#endif
