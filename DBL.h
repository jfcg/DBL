#ifndef _JYS_DBL_H_
#define _JYS_DBL_H_ 1

/*	The DBL Multi-term Extended Precision Floating Point Arithmetic Library
 *	JiYuSoft Ltd. All rights reserved.
 *	http://jiyusoft.com
 *
 *	This software can be used, modified or distributed under the terms of
 *	The GNU Public License v3, available at http://www.gnu.org/licenses/gpl-3.0.html
*/

#ifdef _GLIBCXX_IOSTREAM
#include <cstdio> // for printElm()
#endif

namespace jys {

template <int n>
class dbl;

template <int n> // unnormalized temporary version of dbl<n>, for efficiency
class udbl {
	mutable double d [n+1];
	udbl & operator= (const udbl &);
public:
	explicit udbl (const double * p)
	{
		for (int i=0; i<=n; i++) d[i] = p[i];
	}

	double operator[] (int i) const // retrieve element
	{
		return d[i];
	}

	dbl<n> norm() const; // normalize

	udbl operator- () const // unary minus
	{
		udbl s;
		for (int i=0; i<=n; i++) s.d[i] = -d[i];
		return s;
	}

	dbl<n> operator+ (double r) const; // + - with double

	dbl<n> operator- (double r) const;


	template <int k>
	dbl<n> operator+ (const udbl<k> & r) const; // + - with udbl<k>

	template <int k>
	dbl<n> operator- (const udbl<k> & r) const;


	udbl operator* (double r) const // * / with double
	{ return norm() * r; }

	udbl operator/ (double r) const
	{ return norm() / r; }


	template <int k>
	udbl<n> operator* (const udbl<k> & r) const // * / with udbl<k>
	{ return norm() * r; }

	template <int k>
	udbl<n> operator/ (const udbl<k> & r) const
	{ return norm() / r; }

	// a udbl has to be saved to a dbl before comparing it to a dbl/double
}; // udbl

template <int n> // the DBL data type
class dbl {
	double d [n];
public:
	dbl () {} // constructors

	dbl (double r)
	{ *this = r; }

	dbl (long r)
	{ *this = r; }

	template <int k>
	dbl<n> (const dbl<k> & r)
	{ *this = r; }

	template <int k>
	dbl<n> (const udbl<k> & r)
	{ *this = r; }

	explicit dbl (const char * c) // construct from a string like "-3.141592653589793238462643383279502e-123"
	{ *this = c; }
	
	explicit dbl (const double * p) // p array should be normalized
	{
		for (int i=0; i<n; i++) d[i]=p[i];
	}


	double operator[] (int i) const // retrieve element
	{
		return d[i];
	}

	dbl & operator= (double r) // assignments
	{
		d[0]=r;
		for (int i=1; i<n; i++) d[i]=0;
		return *this;
	}

	dbl & operator= (long r)
	{
		if (sizeof(long)==4) // compile-time constant
		{
			d[0]=r;
			for (int i=1; i<n; i++) d[i]=0;
		}else{
			union {
				long a;
				unsigned int b[2];
			};
			bool neg = r<0;
			if (neg) r = -r;
			a=r;

			d[0]=b[1];
			d[1]=b[0];
			for (int i=2; i<n; i++) d[i]=0;
			d[0]*=0x1p32;
			double h = d[0]+d[1];
			d[1] -= h-d[0];
			d[0] = h;
			if (neg)
			{
				d[0]=-d[0];
				d[1]=-d[1];
			}
		}
		return *this;
	}

	template <int k>
	dbl<n> & operator= (const dbl<k> & r)
	{
		char dummy[n-k]; // we expect n>=k
		int i;
		for (i=0; i<k; i++) d[i]=r[i];
		for (   ; i<n; i++) d[i]=0;
		return *this;
	}

	template <int k>
	dbl<n> & operator= (const udbl<k> & r)
	{ return *this = r.norm(); }

	dbl & operator= (const char * c);


	const double * get() const // return element array
	{ return d; }


	dbl operator- () const // unary minus
	{
		dbl s;
		for (int i=0; i<n; i++) s.d[i] = -d[i];
		return s;
	}


	udbl<n> operator+ (double r) const // additions
	{
		double x[n+1];
		x[0]=d[0]; x[1]=r;
		for (int i=1; i<n; i++) x[i+1]=d[i];
		return udbl<n>(x);
	}

	void operator+= (double r)
	{ *this = *this + r; }

	template <int k>
	dbl<n> operator+  (const dbl<k> & r) const;

	template <int k>
	void   operator+= (const dbl<k> & r)
	{ *this = *this + r; }

	template <int k>
	dbl<n> operator+  (const udbl<k> & r) const;

	template <int k>
	void   operator+= (const udbl<k> & r)
	{ *this = *this + r; }


	udbl<n> operator- (double r) const // subtractions
	{
		double x[n+1];
		x[0]=d[0]; x[1]=-r;
		for (int i=1; i<n; i++) x[i+1]=d[i];
		return udbl<n>(x);
	}

	void operator-= (double r)
	{ *this = *this - r; }

	template <int k>
	dbl<n> operator-  (const dbl<k> & r) const;

	template <int k>
	void   operator-= (const dbl<k> & r)
	{ *this = *this - r; }

	template <int k>
	dbl<n> operator-  (const udbl<k> & r) const;

	template <int k>
	void   operator-= (const udbl<k> & r)
	{ *this = *this - r; }


	udbl<n> operator* (double r) const; // multiplications

	void   operator*= (double r)
	{ *this = *this * r; }

	template <int k>
	udbl<n> operator* (const dbl<k> & r) const;

	template <int k>
	void   operator*= (const dbl<k> & r)
	{ *this = *this * r; }

	template <int k>
	udbl<n> operator* (const udbl<k> & r) const
	{ return *this * r.norm(); }

	template <int k>
	void   operator*= (const udbl<k> & r)
	{ *this = *this * r; }


	udbl<n> operator/ (double r) const; // divisions

	void   operator/= (double r)
	{ *this = *this / r; }

	template <int k>
	udbl<n> operator/ (const dbl<k> & r) const;

	template <int k>
	void   operator/= (const dbl<k> & r)
	{ *this = *this / r; }

	template <int k>
	udbl<n> operator/ (const udbl<k> & r) const
	{ return *this / r.norm(); }

	template <int k>
	void   operator/= (const udbl<k> & r)
	{ *this = *this / r; }


	void ewmul (double m) // element-wise multiply (m is likely to be [the negative of] a power of two)
	{
		for (int i=0; i<n; i++) d[i]*=m;
	}


	bool operator== (double r) const // comparisons
	{
		if (d[0] != r) return false;
		for (int i=1; i<n; i++) if (d[i] != 0) return false;
		return true;
	}
	bool operator!= (double r) const
	{ return ! (*this == r); }

	template <int k>
	bool operator== (const dbl<k> & r) const
	{
		char dummy[n-k]; // we expect n>=k
		int i;
		for (i=0; i<k; i++) if (d[i] != r[i]) return false;
		for (   ; i<n; i++) if (d[i] != 0)    return false;
		return true;
	}
	template <int k>
	bool operator!= (const dbl<k> & r) const
	{ return ! (*this == r); }


	bool operator< (double r) const // <
	{
		if (d[0] < r) return true;
		if (d[0] > r) return false;
		if (d[1] < 0) return true;
		return false;
	}
	bool operator>= (double r) const
	{ return ! (*this < r); }

	template <int k>
	bool operator< (const dbl<k> & r) const
	{
		char dummy[n-k];
		for (int i=0; i<k; i++)
		{
			if (d[i] < r[i]) return true;
			if (d[i] > r[i]) return false;
		}
		if (k<n && d[k] < 0) return true;
		return false;
	}
	template <int k>
	bool operator>= (const dbl<k> & r) const
	{ return ! (*this < r); }


	bool operator> (double r) const // >
	{
		if (d[0] > r) return true;
		if (d[0] < r) return false;
		if (d[1] > 0) return true;
		return false;
	}
	bool operator<= (double r) const
	{ return ! (*this > r); }

	template <int k>
	bool operator> (const dbl<k> & r) const
	{
		char dummy[n-k];
		for (int i=0; i<k; i++)
		{
			if (d[i] > r[i]) return true;
			if (d[i] < r[i]) return false;
		}
		if (k<n && d[k] > 0) return true;
		return false;
	}
	template <int k>
	bool operator<= (const dbl<k> & r) const
	{ return ! (*this > r); }


/*	c[] should have at least 16n+10 chars. returns number of chars written (excluding null).
	writes up to 16n+2 significant digits
	flags & 1: showpos
	flags & 2: uppercase
	flags & 4: hex	// not implemented yet
	flags & 8: print trailing zeros
*/
	int write (char * c, int flags=0) const;
 
	// #include <iostream> before #include "DBL.h"; printElm() and operator<< will be available
#ifdef _GLIBCXX_IOSTREAM
	std::ostream & printElm (std::ostream & o) const // print elements
	{
		if (o.flags() & std::ios_base::hex) // providing hex output as well
		{
			char p[5] = "%la ", buf[26];
			if (o.flags() & std::ios_base::uppercase)
				p[2]='A';
			for (int i=0; i<n-1; i++)
				o.write (buf, std::sprintf (buf, p, d[i]));
			p[3]=0;
			return o.write (buf, std::sprintf (buf, p, d[n-1]));
		}
		for (int i=0; i<n-1; i++)
			o << d[i] << ' ';
		return o << d[n-1];
	}
#endif
}; // dbl


#ifdef _GLIBCXX_IOSTREAM
template <int n>
std::ostream & operator<< (std::ostream & o, const dbl<n> & r)
{
	char buf [16*n+10];
	int f=0;
	if (o.flags() & std::ios_base::showpos)   f+=1;
	if (o.flags() & std::ios_base::uppercase) f+=2;
	return o.write (buf, r.write (buf, f));
}

template <int n>
std::ostream & operator<< (std::ostream & o, const udbl<n> & r)
{ return o << r.norm(); }
#endif


template <int n, int k>
dbl<n> operator+ (const udbl<k> & a, const dbl<n> & r) // arithmetic with udbl
{ return r+a; }

template <int n, int k>
dbl<n> operator- (const udbl<k> & a, const dbl<n> & r)
{ return -r+a; }

template <int n, int k>
udbl<n> operator* (const udbl<k> & a, const dbl<n> & r)
{ return a.norm() * r; }

template <int n, int k>
udbl<n> operator/ (const udbl<k> & a, const dbl<n> & r)
{ return a.norm() / r; }


template <int n>
dbl<n> operator+ (double a, const dbl<n> & r) // arithmetic with double
{ return r+a; }

template <int n>
dbl<n> operator- (double a, const dbl<n> & r)
{ return -r+a; }

template <int n>
udbl<n> operator* (double a, const dbl<n> & r)
{ return r*a; }

template <int n>
udbl<n> operator/ (double a, const dbl<n> & r)
{ return dbl<n>(a)/r; }


template <int n>
dbl<n> operator+ (double a, const udbl<n> & r) // arithmetic with double 2
{ return r+a; }

template <int n>
dbl<n> operator- (double a, const udbl<n> & r)
{ return -r+a; }

template <int n>
udbl<n> operator* (double a, const udbl<n> & r)
{ return r*a; }

template <int n>
udbl<n> operator/ (double a, const udbl<n> & r)
{ return dbl<n>(a)/r; }


template <int n>
bool operator== (double a, const dbl<n> & r) // comparisons with double
{ return r==a; }

template <int n>
bool operator!= (double a, const dbl<n> & r)
{ return r!=a; }

template <int n>
bool operator<  (double a, const dbl<n> & r)
{ return r>a; }

template <int n>
bool operator<= (double a, const dbl<n> & r)
{ return r>=a; }

template <int n>
bool operator>  (double a, const dbl<n> & r)
{ return r<a; }

template <int n>
bool operator>= (double a, const dbl<n> & r)
{ return r<=a; }


 // simple functions
template <int n>
udbl<n> sqr  (const dbl<n> & r); // square
template <int n>
udbl<n> sqr  (const udbl<n> & r)
{ return sqr (r.norm()); }

template <int n>
udbl<n> sqrt (const dbl<n> & r); // square root
template <int n>
udbl<n> sqrt (const udbl<n> & r)
{ return sqrt (r.norm()); }

template <int n>
dbl<n> abs  (const dbl<n> & r)  // absolute value
{
	if (r[0] < 0) return -r;
	return r;
}

template <int n>
dbl<n> abs  (const udbl<n> & r)
{
	dbl<n> x=r.norm();
	if (x[0]<0) return -x;
	return x;
}

template <int n>
bool  isp2  (const dbl<n> & r)  // is (negative of) a power of 2?
{
	union {
		double b;
		int a[2];
	};
	b = r[0];
	return a[0]==0 && (a[1] & 0x000fffff)==0 && r[1]==0;
}


 // polynomial given by coeffs, c0*x^0 1st, # of coeffs, value to evaluate at
template <int n>
dbl<n> poly (const dbl<n> * c, int N, const dbl<n> & x);
template <int n>
dbl<n> poly (const dbl<n> * c, int N, const udbl<n> & x)
{ return poly (c, N, x.norm()); }

 // same input as above: derivative of polynomial given by coeffs, # of coeffs, value to evaluate the derivative at
template <int n>
dbl<n> dpoly (const dbl<n> * c, int N, const dbl<n> & x);
template <int n>
dbl<n> dpoly (const dbl<n> * c, int N, const udbl<n> & x)
{ return dpoly (c, N, x.norm()); }

 // root of a polynomial via modified Regula Falsi, two distinct initial points
template <int n>
dbl<n> polyroot (const dbl<n> * c, int N, const dbl<n> & x0, const dbl<n> & x1, double ltol=2.5);
template <int n>
dbl<n> polyroot (const dbl<n> * c, int N, const udbl<n> & x0, const dbl<n> & x1, double ltol=2.5)
{ polyroot (c, N, x0.norm(), x1, ltol); }
template <int n>
dbl<n> polyroot (const dbl<n> * c, int N, const dbl<n> & x0, const udbl<n> & x1, double ltol=2.5)
{ polyroot (c, N, x0, x1.norm(), ltol); }
template <int n>
dbl<n> polyroot (const dbl<n> * c, int N, const udbl<n> & x0, const udbl<n> & x1, double ltol=2.5)
{ polyroot (c, N, x0.norm(), x1.norm(), ltol); }


#ifdef _JYS_DBL_INLINE_
#include "DBLbody.h"
#endif

} // namespace jys

#endif // _JYS_DBL_H_
