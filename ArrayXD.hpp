// ArrayXD.h
#ifndef MATRIX_H_
#define MATRIX_H_
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

template <typename T>
ArrayXD <T>::ArrayXD(uint n0,uint n1,uint n2, uint n3)
{
	N0=n0;N1=n1;N2=n2;N3=n3;
	Dim =4;
	if (N0==1){Dim--;}
	if(N1==1){Dim--;}
	if(N2==1){Dim--;}
	if(N3==1){Dim--;}
	Size=N0*N1*N2*N3;
	m=new T [Size];
	for(uint i=0;i<N0;i++)
	{
		for(uint j=0;j<N1;j++)
		{
			for(uint k=0;k<N2;k++)
			{
				for(uint r=0;r<N3;r++)
				{
					m[r+N3*(k+N2*(j+N1*i))]=0;
				}
			}
		}
	}

}
template <typename T>
ArrayXD<T>::ArrayXD()
{
	N0=1;N1=1;N2=1;N3=1;
	Dim =4;
	if (N0==1){Dim--;}
	if(N1==1){Dim--;}
	if(N2==1){Dim--;}
	if(N3==1){Dim--;}
	Size=N0*N1*N2*N3;
	m=new T [Size];
	m[0]=0;

}
template <typename T>
ArrayXD<T>::ArrayXD(const ArrayXD& _array)
{

		N0=_array.N0;N1=_array.N1;
		N2=_array.N2;N3=_array.N3;
		Dim =4;
		if(N1==1){Dim--;}
		if(N2==1){Dim--;}
		if(N3==1){Dim--;}
		Size=N0*N1*N2*N3;
		m=new T [Size];

    for (uint i = 0; i < Size; i++) {
		m[i]=_array.m[i];
    }

}
template <typename T>
ArrayXD<T>::~ArrayXD()
{
	delete [] m;

}

template <typename T>
void ArrayXD<T>::Show() const
{
	for(uint i=0;i<N0;i++)
	{   if (N1!=1){std::cout<<"("<<i<<",:)"<<endl;}
		for(uint j=0;j<N1;j++)
		{
			if (N2!=1){std::cout<<"("<<i<<","<<j<<",:)"<<endl;}
			for(uint k=0;k<N2;k++)
			{
				if (N3!=1){cout<<"("<<i<<","<<j<<","<<k<<":)"<<endl;}
				for(uint r=0;r<N3;r++)
				{
					std::cout<<m[r+N3*(k+N2*(j+N1*i))]<<'\t';
				}
				if (N3!=1){std::cout<<endl;}
			}
			if (N2!=1){std::cout<<endl;}
		}
		if (N1!=1){std::cout<<endl;}
	}

}

template <typename T>
T & ArrayXD<T> ::operator()(uint ind0, uint ind1, uint ind2,uint ind3){
	return m[ind3+N3*(ind2+N2*(ind1+N1*ind0))];
}

template <typename T>
ArrayXD<T> ArrayXD<T> ::operator()(const ArrayXD<int> &ind0, const ArrayXD<int>  &ind1, const ArrayXD<int>  &ind2,const ArrayXD<int> &ind3){
	
	
	ArrayXD<T> temp(ind0.N0,ind1.N0,ind2.N0,ind3.N0);
	for (uint i=0;i<temp.N0;i++){
		for(uint j=0;j<temp.N1;j++){
			for (uint k = 0; k < temp.N2; k++){
				for (uint r = 0; r < temp.N3; r++){
					temp.m[r+temp.N3*(k+temp.N2*(j+temp.N1*i))]=m[ind3.m[r]+N3*(ind2.m[k]+N2*(ind1.m[j]+N1*ind0.m[i]))];
				}
			}
		}
	}
	return temp;
}


template <typename T>
ArrayXD<T> & ArrayXD<T>::operator=(const ArrayXD& _array)
{
    if (this == &_array) {
        return *this;
    }

    if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		delete [] m;
		N0=_array.N0;N1=_array.N1;
		N2=_array.N2;N3=_array.N3;
		Dim =4;
		if(N1==1){Dim--;}
		if(N2==1){Dim--;}
		if(N3==1){Dim--;}
		Size=N0*N1*N2*N3;
		m=new T [Size];
    }

    for (uint i = 0; i < Size; i++) {
		m[i]=_array.m[i];
    }
    return *this;
}

template <typename T>
ArrayXD<T> & ArrayXD<T>::operator=(const T & val)
{   
    for (uint i = 0; i < Size; i++) {
		m[i]=val;
    }
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator+=(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	else{
		for (uint i=0;i<Size;i++){
			m[i]+=_array.m[i];
		}

	}
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator-=(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	else{
		for (uint i=0;i<Size;i++){
			m[i]-=_array.m[i];
		}

	}
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator*=(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	else{
		for (uint i=0;i<Size;i++){
			m[i]*=_array.m[i];
		}

	}
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator/=(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	else{
		for (uint i=0;i<Size;i++){
			m[i]/=_array.m[i];
		}

	}
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator+=(const T & val)
{
	
		for (uint i=0;i<Size;i++){
			m[i]+=val;
		}
    return *this;
}

template <typename T>
ArrayXD<T>& ArrayXD<T>::operator-=(const T & val)
{
	
		for (uint i=0;i<Size;i++){
			m[i]-=val;
		}
    return *this;
}
template <typename T>
ArrayXD<T>& ArrayXD<T>::operator*=(const T & val)
{
	
		for (uint i=0;i<Size;i++){
			m[i]*=val;
		}
    return *this;
}
template <typename T>
ArrayXD<T>& ArrayXD<T>::operator/=(const T & val)
{
	
		for (uint i=0;i<Size;i++){
			m[i]/=val;
		}
    return *this;
}

template <typename T>
ArrayXD<T> ArrayXD<T>::operator+(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	
		ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]+_array.m[i];
		}
		 return temp;
	

}

template <typename T>
ArrayXD<T> ArrayXD<T>::operator-(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	
		ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]-_array.m[i];
		}
		 return temp;
	

}
template <typename T>
ArrayXD<T> ArrayXD<T>::operator*(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	
		ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]*_array.m[i];
		}
		 return temp;
	

}

template <typename T>
ArrayXD<T> ArrayXD<T>::operator/(const ArrayXD& _array)
{
	if (N0!= _array.N0 || N1!= _array.N1|| N2!= _array.N2|| N3!= _array.N3) {
		std::cout<<"xxxx"<<endl;
    }
	
		ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]/_array.m[i];
		}
		 return temp;	
}

template <typename T>
ArrayXD<T> ArrayXD<T>::operator+(const T & val)
{
	  ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]+val;
		}
    return temp;
}

template <typename T>
ArrayXD<T> ArrayXD<T>::operator-(const T & val)
{
	  ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]-val;
		}
    return temp;
}
template <typename T>
ArrayXD<T> ArrayXD<T>::operator*(const T & val)
{
	  ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]*val;
		}
    return temp;
}
template <typename T>
ArrayXD<T> ArrayXD<T>::operator/(const T & val)
{
	  ArrayXD<T> temp(N0,N1,N2,N3);
		for (uint i=0;i<Size;i++){
			temp.m[i]=m[i]/val;
		}
    return temp;
}


//------------------------
template <typename T>
ArrayXD<T> operator+(const T & val,const ArrayXD<T>& _array)
{
	  ArrayXD<T> temp(_array);
		for (uint i=0;i<temp.Size;i++){
			temp.m[i]+=val;
		}
    return temp;
}

template <typename T>
ArrayXD<T> operator-(const T & val,const ArrayXD<T>& _array)
{
	  ArrayXD<T> temp(_array);
		for (uint i=0;i<temp.Size;i++){
			temp.m[i]-=val;
		}
    return temp;
}
template <typename T>
ArrayXD<T> operator*(const T & val,const ArrayXD<T>& _array)
{
	  ArrayXD<T> temp(_array);
		for (uint i=0;i<temp.Size;i++){
			temp.m[i]*=val;
		}
    return temp;
}
template <typename T>
ArrayXD<T> operator/(const T & val,const ArrayXD<T>& _array)
{
	  ArrayXD<T> temp(_array);
		for (uint i=0;i<temp.Size;i++){

			temp.m[i]=val/(temp.m[i]);
			
		}
    return temp;
}

ArrayXD<int> range(const int &ind0, const int &ind1,const int &dind){
		int Num=(ind1-ind0)/dind+1;
		ArrayXD<int> temp(Num);
		for (int i=0;i<Num;i++){
		temp(i)=ind0+i*dind;
		}
		return temp;

}

ArrayXD<int> range(const int &ind0){
		int Num=ind0;
		ArrayXD<int> temp(Num);
		for (int i=0;i<Num;i++){
		temp(i)=i;
		}
		return temp;

}

#endif
