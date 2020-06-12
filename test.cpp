//ArrayXD.cpp
#include<iostream>
#include <complex>
#include <valarray>
#include"ArrayXD.hpp"


int main()
{

	ArrayXD<double> A(3,3);
	A=3.;
	A.Show();
	ArrayXD<double> B=1;
	B.Show();
	A*=2.;
	ArrayXD<double> C;
	C=1.+1./A+4.;
	C.Show();
	ArrayXD<double> D;
	C(range(2),range(2))=0;
	
	C.Show();
	valarray<double> a(2,2);
	cout<<asin(a[1])<<endl;
	return 0;
}
