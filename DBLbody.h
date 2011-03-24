#ifndef _JYS_DBL_H_
#error "Include DBL.h instead of this file!"
#endif

/*	The DBL Multi-term Extended Precision Floating Point Arithmetic Library
 *	JiYuSoft Ltd. All rights reserved.
 *	http://jiyusoft.com
 *
 *	This software can be used, modified or distributed under the terms of
 *	The GNU Public License v3, available at http://www.gnu.org/licenses/gpl-3.0.html
 *
 *	Base functions for multi-term fp arithmetic are described in the paper
 *	"Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates"
 *	by J. R. Shewchuk
*/

static inline void _split (double a, double & hi, double & lo) // splits a into high and low words, assumes |a| <= 0x1.ffffffp996
{
	double t = 134217729 * a; // 1+2^27
	hi = t - (t - a);
	lo = a - hi;
}

static inline double _tp (double a, double b, double & e) // Computes fl(a*b) and err(a*b)
{
	double p, ah, al, bh, bl;
	_split (a, ah, al);
	_split (b, bh, bl);
	e = ((ah * bh - (p=a*b)) + ah * bl + al * bh) + al * bl;
	return p;
}

static inline double _tps (double a, double & e) // Computes fl(a*a) and err(a*a)
{
	double p, h, l;
	_split (a, h, l);
	e = ((h * h - (p=a*a)) + 2 * h * l) + l * l;
	return p;
}

static inline double _qts (double a, double b, double & e) // Computes fl(a+b) and err(a+b), assumes exp(a)>=exp(b) or a=0
{
	double r=a+b;
	e = b - (r-a);
	return r;
}

static inline double _ts (double a, double b, double & e) // Computes fl(a+b) and err(a+b)
{
	double r=a+b, c=r-a;
	e = (a - (r-c)) + (b-c);
	return r;
}

static inline double _td (double a, double b, double & e) // Computes fl(a-b) and err(a-b)
{
	double r=a-b, c=r-a;
	e = (a - (r-c)) - (b+c);
	return r;
}

static inline bool _over (double p) // is nan or inf?
{
	union {
		double b;
		int a[2];
	};
	b=p;
	return (a[1] & 0x7ff00000)==0x7ff00000;
}

static inline void _normalize (double * x, int n, bool ch0=true) // normalizes a mixed-range array of doubles into a unique expansion
{
	if (ch0 && _over(x[0])) return;
	if (_over(x[1]))
	{
		double t = x[0];
		x[0] = x[1];
		x[1] = t;
		return;
	}
reStart:
	if (ch0 && x[0]==0) // handle 1st
	{
		int j;
		for (j=1; j<n; j++) // find next non-zero
			if (x[j]!=0) break;
		if (j>=n) return;

		x[0] = x[j];
		x[j] = 0;
	}

	for (int i=1; i<n; i++)
	{
		if (x[i]==0)
		{
			int j;
			for (j=i+1; j<n; j++) // find next non-zero
				if (x[j]!=0) break;
			if (j>=n) return;

			x[i] = x[j];
			x[j] = 0;
		}

		double r = x[i-1] + x[i];

		if (x[i-1] != r) // need a ts?
		{
			double	c = r - x[i-1],
				e = (x[i] - c) + (x[i-1] - (r-c));
			x[i-1] = r;
			x[i]   = e;
			i-=2;
			if (i<0) goto reStart;
		}
	}
}


template <int n>
dbl<n> udbl<n>::norm() const
{
	_normalize (d, n+1);
	return dbl<n>(d);
}

template <int n>
dbl<n> udbl<n>::operator+ (double r) const
{
	double x[n+2];
	x[0]=d[0]; x[1]=r;
	for (int i=1; i<=n; i++) x[i+1]=d[i];
	_normalize (x, n+2);
	return dbl<n>(x);
}

template <int n>
dbl<n> udbl<n>::operator- (double r) const
{
	double x[n+2];
	x[0]=d[0]; x[1]=-r;
	for (int i=1; i<=n; i++) x[i+1]=d[i];
	_normalize (x, n+2);
	return dbl<n>(x);
}

template <int n>
template <int k>
dbl<n> udbl<n>::operator+ (const udbl<k> & r) const
{
	char dummy[n-k]; // we expect n>=k
	double x[n+k+2];
	int i;
	for (i=0; i<=k; i++)
	{
		x[2*i]  =d[i];
		x[2*i+1]=r[i];
	}
	for (; i<=n; i++)
		x[i+k+1]=d[i];

	_normalize(x, n+k+2); // let normalize do all the work :)
	return dbl<n>(x);
}

template <int n>
template <int k>
dbl<n> udbl<n>::operator- (const udbl<k> & r) const
{
	char dummy[n-k];
	double x[n+k+2];
	int i;
	for (i=0; i<=k; i++)
	{
		x[2*i]  = d[i];
		x[2*i+1]=-r[i];
	}
	for (; i<=n; i++)
		x[i+k+1]=d[i];

	_normalize(x, n+k+2);
	return dbl<n>(x);
}


template <int n>
template <int k>
dbl<n> dbl<n>::operator+ (const dbl<k> & r) const // additions
{
	char dummy[n-k]; // we expect n>=k
	double x[n+k];
	int i;
	for (i=0; i<k; i++)
	{
		x[2*i]  =d[i];
		x[2*i+1]=r[i];
	}
	for (; i<n; i++)
		x[i+k]=d[i];

	_normalize(x, n+k);
	return dbl<n>(x);
}

template <int n>
template <int k>
dbl<n> dbl<n>::operator+ (const udbl<k> & r) const
{
	char dummy[n-k];
	double x[n+k+1];
	int i;
	for (i=0; i<k; i++)
	{
		x[2*i]  =r[i];
		x[2*i+1]=d[i];
	}
	x[2*k] = r[k];
	for (; i<n; i++)
		x[i+k+1]=d[i];

	_normalize(x, n+k+1);
	return dbl<n>(x);
}


template <int n>
template <int k>
dbl<n> dbl<n>::operator- (const dbl<k> & r) const // subtractions
{
	char dummy[n-k]; // we expect n>=k
	double x[n+k];
	int i;
	for (i=0; i<k; i++)
	{
		x[2*i]  = d[i];
		x[2*i+1]=-r[i];
	}
	for (; i<n; i++)
		x[i+k]=d[i];

	_normalize(x, n+k);
	return dbl<n>(x);
}

template <int n>
template <int k>
dbl<n> dbl<n>::operator- (const udbl<k> & r) const
{
	char dummy[n-k];
	double x[n+k+1];
	int i;
	for (i=0; i<k; i++)
	{
		x[2*i]  =-r[i];
		x[2*i+1]= d[i];
	}
	x[2*k] =-r[k];
	for (; i<n; i++)
		x[i+k+1]=d[i];

	_normalize(x, n+k+1);
	return dbl<n>(x);
}


template <int n>
static inline void _dmul (double * x, const dbl<n> & d, double r) // x = d * r
{
	double p1, p2; // for n=2..5, there are up to 6 (>=n+1) addend from degree 0,1
	// tp
	// ts x 2
	// +=
	x[0] = _tp (d[0], r , x[1]);
	p1   = _tp (d[1], r , x[2]);
	x[1] = _ts (x[1], p1, p1  ); // p1 1st ts

	if (n==2) // compile time constant
	{
		x[2] += p1;
		x[1] += x[2];

		x[0] = _qts(x[0], x[1] , x[1]); // final
		return;
	}
	x[2] = _ts (x[2], p1, x[3]); // p1 2nd ts & done

	p1   = _tp (d[2], r , p2  );
	x[2] = _ts (x[2], p1, p1  ); // p1 1st ts

	if (n==3)
		x[3] += p1 + p2;
	else{
		x[3] = _ts (x[3], p1, x[4]); // p1 2nd ts & done
		x[3] = _ts (x[3], p2, p2  ); // p2 1st ts

		if (n==4)
			x[4] += p2;
		else
			x[4] = _ts (x[4], p2, x[5]);
	}
	if (n>5) x[6]=0;

	for (int i=3; i<n; i++) // degree i
	{
		p1   = _tp (d[i], r , p2); // two new values to insert
		int t;
		for (t=i;   t<n && t<i+2; t++)
			x[t] = _ts (x[t], p1, p1);
		x[t]+=p1;

		for (t=i+1; t<n && t<i+3; t++)
			x[t] = _ts (x[t], p2, p2);
		x[t]+=p2;
	}
}

template <int n>
udbl<n> dbl<n>::operator* (double r) const // multiplications
{
	double x[n+1];
	_dmul (x, *this, r);
	return udbl<n>(x);
}

template <int n, int k>
static inline udbl<n> _mul (const double * d, const dbl<k> & r, int irst=1) // used by operator* and write(char*)
{
	double x[n+1], p1, p2; // start with degree 0,1
	// tp
	// ts x 2
	// +=
	x[0] = _tp (d[0], r[0], x[1]);
	p1   = _tp (d[0], r[1], x[2]);
	x[1] = _ts (x[1], p1  , p1  ); // p1 1st ts

	if (n==2) // compile-time constant
		x[2] += p1;
	else
		x[2] = _ts (x[2], p1, x[3]); // p1 2nd ts & done

	p1   = _tp (d[1], r[0], p2  );
	x[1] = _ts (x[1], p1  , p1  ); // p1 1st ts

	if (n==2)
	{
		x[2] += p1 + p2 + d[1] * r[1];
	//	x[1] += x[2]; // _qts(x[1], x[2] , x[2]);

	//	x[0] = _qts(x[0], x[1] , x[1]); // final
		return udbl<n>(x);
	}
	x[0] = _qts(x[0], x[1], x[1]);

	x[2] = _ts (x[2], p1  , p1  ); // p1 2nd ts
	x[3] += p1;
	x[2] = _ts (x[2], p2  , p2  ); // p2 1st ts

	if (n==3)
		x[3] += p2;
	else
		x[3] = _ts (x[3], p2 , x[4]); // p2 2nd ts & done

	if (n>4) x[5]=0;
	if (n>5) x[6]=0;

	for (int i=2; i<n; i++) // degree i
	{
		x[i-1] = _qts (x[i-1], x[i], x[i]);

		for (int ir=0; ir<=i && ir<k; ir++)
		{
			p1 = _tp (d[i-ir], r[ir], p2); // two new values to add
			int t;
			for (t=i;   t<n && t<i+2; t++)
				x[t] = _ts (x[t], p1, p1);
			x[t]+=p1;

			for (t=i+1; t<n && t<i+3; t++)
				x[t] = _ts (x[t], p2, p2);
			x[t]+=p2;
		}
	}

	x[n-1] = _qts (x[n-1], x[n], x[n]);

	for (int ir=irst; ir<k; ir++) // degree n // write(char*) uses long constants, hence will use irst=0
		x[n] += d[n-ir] * r[ir];

	return udbl<n>(x);
}

template <int n>
template <int k>
udbl<n> dbl<n>::operator* (const dbl<k> & r) const
{
	char dummy[n-k]; // we expect n>=k
	return _mul<n,k> (d, r);
}


static inline void _d2r (double * x, const double * d, double r, bool nr=false) // dbl<2>(x) = dbl<2>(d) / r
{
	double s[2];
	dbl<2> rem;

	x[0] = d[0] / r; // q1
	s[0] = _tp (x[0], r, s[1]);
	rem = dbl<2>(d) - dbl<2>(s);

	x[1] = rem[0] / r; // q2
	s[0] = _tp (x[1], r, s[1]);
	rem -= dbl<2>(s);

	if (nr)
	{
		x[2] = rem[0] / r;
		return;
	}

	s[0] = rem[0] / r; // q3
	x[0] = _qts (x[0], x[1], x[1]);
	x[1] += s[0];
	x[0] = _qts (x[0], x[1], x[1]);
}

template <int n>
udbl<n> dbl<n>::operator/ (double r) const // divisions
{
	double x[n+1]; // will be filled 2 at each iter
	if (n==2)
	{
		_d2r (x, d, r, true);
		return udbl<n>(x);
	}
	double z, y[n+3]; // will hold remainder and q[i] * -r
	for (int i=0; i<n; i++) y[i]=d[i];

	for (int i=0; i<n-1; i+=2)
	{
		_d2r (x+i, y, r);
		y[n]   = _tp (x[i]  ,-r, y[n+1]); // q[i] * -r
		z      = _tp (x[i+1],-r, y[n+2]);
		y[n+1] = _ts (y[n+1], z, z     );
		y[n+2] += z;

		z=y[1];
		y[1]=y[n];
		if (n==2) y[2]=z; // mix a little bit
		else if (n==3)
		{
			y[3]=y[4];
			y[4]=y[2];
			y[2]=z;
		}else{
			y[n]=y[2];
			y[2]=z;
			z=y[3];
			y[3]=y[n+1];
			y[n+1]=z;
		}
		_normalize(y, n+3);
	}
	if (n%2 == 0) // last single digit of x ?
		x[n] = y[0] / r;
	else
		_d2r (x+n-1, y, r);

	return udbl<n>(x);
}

template <int n>
template <int k>
udbl<n> dbl<n>::operator/ (const dbl<k> & r) const
{
	double z, x[n+1], y[n+k+1]; // q's, holds rem and r * -q_i
	for (int i=0; i<n; i++) y[i]=d[i];

	for (int i=0; i<n; i++)
	{
		x[i] = y[0] / r[0]; // q_i
		_dmul (y+n, r, -x[i]); // r * -q_i

		z=y[1];
		y[1]=y[n];
		if (n==2) y[2]=z; // mix a little bit
		else if (n==3)
		{
			y[3]=y[4];
			y[4]=y[2];
			y[2]=z;
		}else{
			y[n]=y[2];
			y[2]=z;
			z=y[3];
			y[3]=y[n+1];
			y[n+1]=z;
		}
		_normalize (y, n+k+1);
	}
	x[n] = y[0] / r[0]; // q_n

	return udbl<n>(x);
}

/*
template <int n>
static inline dbl<n> _dsum (double a, const dbl<n> & b) // non-zero a placed 1st
{
	double x[n+1];
	x[0]=a;
	for (int i=0; i<n; i++) x[i+1]=b[i];
	_normalize (x, n+1, false);
	return dbl<n>(x);
}

template <int n>
template <int k>
dbl<n> dbl<n>::operator/ (const dbl<k> & r) const // modified Goldschmidt
{
	char dummy[n-k]; // we expect n>=k
	union {
		double b;
		int a[2];
	};
	b = r[0];
	int s = a[1] & 0x80000000; // sign
	int e = a[1] & 0x7ff00000; // exp
	if (r[1]==0) // only r[0] ?
	{
		if (a[0] || s+e!=a[1])
			return *this / r[0];

		dbl<n> y = *this; // power of 2
		a[1] = s + 0x7fe00000 - e;
		y.ewmul(b); // normalize by the (negative of) power of 2 s.t. r=1
		return y;
	}

	dbl<n> y = *this, x = r; // looking for y/x
	a[0] = 0;
	a[1] = s + 0x7fd00000 - e;
	y.ewmul(b); // normalize by the (negative of) power of 2 s.t. 0.5 < x < 1
	x.ewmul(b);
	if (_over(y[0])) return y;

	x = _dsum (1.0, -x); // 0 < x < 0.5
	const double thr[5] = {0x1p-54, 0x1p-81, 0x1p-108, 0x1p-135, 0x1p-162}; // 54 n / 2
	goto lpStart; do
	{
		x = sqr (x);
	lpStart:
		y *= _dsum (1.0, x);
	}
	while (x[0]>=thr[n-2]);

	return y;
}
*/

template <int n>
udbl<n> sqr (const dbl<n> & r) // square
{
	double x[n+1], p1, p2; // for n=2..5, degree 0,1
	// tp
	// ts x 2
	// +=
	x[0] = _tps(r[0], x[1]);
	p1   = _tp (r[0], r[1], p2  );
	p1*=2; x[2]=2*p2;
	x[1] = _ts (x[1], p1  , p1  ); // p1 1st ts

	if (n==2)
		x[2] += p1;
	else
		x[2] = _ts (x[2], p1  , x[3]); // p1 2nd ts & done

	if (n==2)
	{
		x[2] += r[1] * r[1];
	//	x[1] += x[2]; // _qts(x[1], x[2] , x[2]);

	//	x[0] = _qts(x[0], x[1] , x[1]); // final
		return udbl<n>(x);
	}
	x[0] = _qts(x[0], x[1], x[1]);

	p1   = _tp (r[0], r[2], p2  );
	p1*=2; p2*=2;
	x[2] = _ts (x[2], p1  , p1  ); // p1 1st ts
	if (n==3)
		x[3] += p1+p2;
	else{
		x[3] = _ts (x[3], p1  , x[4]); // p1 2nd ts & done

		x[3] = _ts (x[3], p2  , p2  ); // p2 1st ts
		if (n==4)
			x[4] += p2;
		else
			x[4] = _ts (x[4], p2  , x[5]); // p2 2nd ts & done
	}

	x[1] = _qts(x[1], x[2], x[2]);

	p1   = _tps(r[1], p2  );
	x[2] = _ts (x[2], p1  , p1  ); // p1 1st ts
	if (n==3)
		x[3] += p1+p2;
	else{
		x[3] = _ts (x[3], p1  , p1  ); // p1 2nd ts & done
		x[4] += p1;

		x[3] = _ts (x[3], p2  , p2  ); // p2 1st ts
		if (n==4)
			x[4] += p2;
		else{
			x[4] = _ts (x[4], p2  , p2  ); // p2 2nd ts & done
			x[5] += p2;
		}
	}
	if (n>5) x[6]=0;

	for (int i=3; i<n; i++) // degree i
	{
		x[i-1] = _qts (x[i-1], x[i], x[i]); // x[i] will see about i packed values now, so ease it

		int ir;
		for (ir=0; ir<(i+1)/2; ir++)
		{
			p1 = _tp (r[i-ir], r[ir], p2); // two new values to insert
			p1*=2; p2*=2;
			int t;
			for (t=i;   t<n && t<i+2; t++)
				x[t] = _ts (x[t], p1, p1);
			x[t]+=p1;

			for (t=i+1; t<n && t<i+3; t++)
				x[t] = _ts (x[t], p2, p2);
			x[t]+=p2;
		}
		if (i%2==0) // (2*ir==i)
		{
			p1 = _tps (r[ir], p2);
			int t;
			for (t=i;   t<n && t<i+2; t++)
				x[t] = _ts (x[t], p1, p1);
			x[t]+=p1;

			for (t=i+1; t<n && t<i+3; t++)
				x[t] = _ts (x[t], p2, p2);
			x[t]+=p2;
		}
	}

	x[n-1] = _qts (x[n-1], x[n], x[n]);

	for (int ir=1; ir<(n+1)/2; ir++) // degree n
		x[n] += 2 * r[n-ir] * r[ir];
	if (n%2==0) // (2*ir==n)
		x[n] += r[n/2] * r[n/2];

	return udbl<n>(x);
}

template <int n>
udbl<n> sqrt (const dbl<n> & r) // square root: solve for r-x^-2=0 then multiply with r
{
	dbl<n> x = 1 / __builtin_sqrt (r[0]); //start with a good approximation
	dbl<n> y = r;
	y.ewmul(-0.5);

	x *= y*sqr(x)+1.5;
	x *= y*sqr(x)+1.5;
	if (n>3) x *= y*sqr(x)+1.5; // compile-time constant

	return x*r;
}


template <int n>
dbl<n> poly (const dbl<n> * c, int N, const dbl<n> & x) // polynomial given by coeffs, # of coeffs, value to evaluate at
{
	dbl<n> r = c[N-1];

	for (int i=N-2; i>=0; i--)
		r = r * x + c[i];
	return r;
}

template <int n>
dbl<n> dpoly (const dbl<n> * c, int N, const dbl<n> & x) // derivative of polynomial given by coeffs, # of coeffs, value to evaluate the derivative
{
	dbl<n> r = c[N-1] * (N-1);

	for (int i=N-2; i>0; i--)
		r = r * x + c[i] * i;
	return r;
}

template <int n>
dbl<n> polyroot (const dbl<n> * c, int N, const dbl<n> & x0,
		const dbl<n> & x1, double ltol) // root of a polynomial via Regula Falsi, two distinct initial points
{
	dbl<n> x=x0, y=x1,
		fx = poly (c, N, x),
		fy = poly (c, N, y),
		z = x - (x-y) * (fx / (fx-fy)), fz;

	for (int i=0; i<n+7; i++)
	{
		dbl<n>	_lx = z-x, // calc diffs in dbl<n> precision
			_ly = z-y;
		double	lx = __builtin_fabs (_lx[0]),
			ly = __builtin_fabs (_ly[0]);

		if (lx==0 || ly==0) break;
		fz = poly (c, N, z);

		bool bxz =  (fx[0]>0 && fz[0]>0) || (fx[0]<0 && fz[0]<0);
		if (bxz && ((fy[0]>0 && fz[0]>0) || (fy[0]<0 && fz[0]<0)))
		{
			if (lx>ly) // fx,fy,fz all have the same sign: need a big step towards beyond z
			{
				y  = z;
				fy = fz;
				x.ewmul(-0.5);
				x += z;
				x.ewmul(2);
				fx = poly (c, N, x);
			}else{
				x  = z;
				fx = fz; 
				y.ewmul(-0.5);
				y += z;
				y.ewmul(2);
				fy = poly (c, N, y);
			}
			--i; // don't count this
		}
		else if (ltol*lx < ly) // if one of x,y is significantly closer to z than the other..
		{
			y  = z;
			fy = fz;
		}
		else if (ltol*ly < lx)
		{
			x  = z;
			fx = fz; 
		}
		else if (!bxz) // pick opposite signs
		{
			y  = z;
			fy = fz;
		}else{
			x  = z;
			fx = fz; 
		}
		z = x - (x-y) * (fx / (fx-fy));
	}
	return z;
}


static const double _t46[2] = {0x1.c06a5ec5433c6p+152, 0x1.bb542c80deb48p+95};
static const double _t69[3] = {0x1.28bc8abe49f64p+229,-0x1.83b80b9aab60cp+175, -0x1.6212b659388cbp+121};
static const double _t94[4] = {0x1.32d17ed577d0cp+312,-0x1.bf02cce9d8616p+256, -0x1.e13456af99773p+201, 0x1.43016e00b2ea9p+146};
static const double _t139[3]= {0x1.adf1aea12525bp+461,-0x1.fcba562d7ba2cp+406, -0x1.dea482e152899p+351};
static const double _t164[5]= {0x1.bc8d3f7b3340cp+544,-0x1.c9035a0712651p+485, -0x1.889b050145ae8p+428, 0x1.e02500835011bp+374,-0x1.567e643425ed8p+317};
static const double _t186[6]= {0x1.d6affe45f818fp+617, 0x1.59a90cb506d15p+562,  0x1.69fa9bd46dd02p+508, 0x1.b1133b2eca5bep+454,
	-0x1.8ef1056744c95p+400, 0x1.6e335704054c8p+342}; // some dbl<5> constants will hold an extra term for fine precision

static const double _tm45[3] = {0x1.6d601ad376ab9p-150, 0x1.a27ac0f72f8cp-206 ,-0x1.796ed97d96558p-260};
static const double _tm67[5] = {0x1.59165a6ddda5bp-223, 0x1.62f139233f1e9p-277, 0x1.179fd4082810dp-333, 0x1.6de90ee04b07ep-389, -0x1.012396b87702p-444};
static const double _tm91[6] = {0x1.a12f5a0f4e3e5p-303,-0x1.4d81dfff5beccp-358, 0x1.a9e9e93281198p-412, 0x1.a16c57fa48c18p-466,
	-0x1.7218317b223acp-521,-0x1.7bf3378d7a216p-575};
static const double _tm136[6]= {0x1.29b69070816e3p-452,-0x1.0573d5bb14ba9p-511, 0x1.cdc1ae679cb4p-567 , 0x1.1a5c3d6a82badp-621,
	0x1.cd2716d51380dp-678 , 0x1.65f5421614bacp-733};
static const double _tm180[6]= {0x1.0991a9bfa58c8p-598,-0x1.89ac264c17ff8p-654, 0x1.e07375a99f794p-709, 0x1.89c739be1cfdp-764 ,
	-0x1.db3ff509895d9p-819,-0x1.b7a6607309fe8p-874};
static const double _tm217[6]= {0x1.1a66c1e139eddp-721,-0x1.a14dc769142a2p-775,-0x1.09a8cd2f453e5p-830,-0x1.01f2618c96b06p-886,
	0x1.da0544e7b6682p-942 ,-0x1.669b85e12efc7p-999};

static const double _lwp[21] = {10,100,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21};

template <int n>
static inline void _fixex (dbl<n> & y, int e) // y *= 10^-e
{
	if (e<0)
	{
		if (n==5)
		if (e<=-186)
		{
			y = _mul<n,n> (_t186, y, 0); // make sure x[n] is good
			e+=186;
		}

		if (n==3 || n==4)
		if (e<=-164)
		{
			y = _mul<n,n> (_t164, y, 0);
			e+=164;
		}

		if (n==2)
		while (e<=-139)
		{
			y = _mul<n,n> (_t139, y, 0);
			e+=139;
		}

		if (e<=-94)
		{
			if (n>3) y = _mul<n,4> (y.get(), dbl<4>(_t94));
			else y = _mul<n,n> (_t94, y, 0);
			e+=94;
		}

		if (e<=-69)
		{
			if (n>2) y = _mul<n,3> (y.get(), dbl<3>(_t69));
			else y = _mul<n,n> (_t69, y, 0);
			e+=69;
		}

		if (e<=-46)
		{
			y = _mul<n,2> (y.get(), dbl<2>(_t46));
			e+=46;
		}

		while (e<=-22)
		{
			y *= 1e22;
			e+=22;
		}
		if (e<0) y *= _lwp[-e-1];
	}
	else if (e>0)
	{
		if (n==5)
		if (e>=217)
		{
			y = _mul<n,n> (_tm217, y, 0); // make sure x[n] is good
			e-=217;
		}

	if (n>2) // compile-time constant
	{
		if (e>=180)
		{
			y = _mul<n,n> (_tm180, y, 0);
			e-=180;
		}

		if (e>=136)
		{
			y = _mul<n,n> (_tm136, y, 0);
			e-=136;
		}
	}else{
		while (e>=136)
		{
			y = _mul<n,n> (_tm136, y, 0);
			e-=136;
		}
	}

		if (e>=91)
		{
			y = _mul<n,n> (_tm91, y, 0);
			e-=91;
		}

		if (n<5)
		if (e>=67)
		{
			y = _mul<n,n> (_tm67, y, 0);
			e-=67;
		}

		if (n==2)
		if (e>=45)
		{
			y = _mul<n,n> (_tm45, y, 0);
			e-=45;
		}

		while (e>=22)
		{
			y /= 1e22;
			e-=22;
		}
		if (e>0) y /= _lwp[e-1];
	}
}

template <int n>
dbl<n> & dbl<n>::operator= (const char * c) // from a string like "-3.141592653589793238462643383279502e-123"
{
	bool neg=false;
	if (*c=='-')
	{
		c++;
		neg=true;
	}
	else if (*c=='+')
		c++;

	int ti=0, e, bd=-1; // total # of digits encountered, # of digits before dot
	*this = 0.0;
	do {
		int i=0; e=0;
		for (; i<9; c++)
		{
			if (*c>='0' && *c<='9')
			{
				e *= 10;
				e += *c - '0';
				i++;
			}
			else if (*c=='.')
			{
				if (bd>=0) break; // 2nd dot?
				bd=ti+i;
			}
			else break;
		}
		if (i==0) break;
		*this = *this * _lwp[i-1] + double(e);
		ti += i;
	}
	while(1);

	e=0;
	if (*c=='e' || *c=='E') // extract exp
	{
		c++;
		bool ng=false;
		if (*c=='-')
		{
			c++;
			ng=true;
		}
		else if (*c=='+')
			c++;
		for (int i=0; i<3 && *c>='0' && *c<='9'; i++, c++)
		{
			e *= 10;
			e += *c - '0';
		}
		if (ng) e=-e;
	}
	if (bd>=0) e -= ti-bd;
	_fixex (*this, -e);
	if (neg) *this = - *this;

	return *this;
}

/*	c[] should have at least 16n+10 chars. returns number of chars written (excluding null).
	writes up to 16n+2 significant digits
	flags & 1: showpos
	flags & 2: uppercase
	flags & 4: hex	// not implemented yet
	flags & 8: print trailing zeros
*/
template <int n>
int dbl<n>::write (char * c, int flags) const
{
	char * p = c;
	union {
		double b;
		int a[2];
	};
	b = d[0];
	bool neg = (a[1] & 0x80000000)!=0;
	if (neg)
		*p++ = '-';
	else if (flags & 1)
		*p++ = '+';

	if ((a[1] & 0x7ff00000)==0x7ff00000) // inf or nan?
	{
		if (a[0] || (a[1] & 0x000fffff))
		{
			if (flags & 2)
				*(int*)p = 'N' + 'A'*256 + 'N'*65536;
			else
				*(int*)p = 'n' + 'a'*256 + 'n'*65536;
		}else{
			if (flags & 2)
				*(int*)p = 'I' + 'N'*256 + 'F'*65536;
			else
				*(int*)p = 'i' + 'n'*256 + 'f'*65536;
		}
		return p+3-c;
	}
	if (d[0]==0)
	{
		*p++ = '0';
		*p = 0;
		return p-c;
	}

	dbl<n> y = *this;
	if (neg) y = -y;

	int e = __builtin_floor (__builtin_log10 (y[0])), eo=e;
	_fixex (y, e);
	
	if (y<1.0)
	{
		y *= 10; // correction
		eo--;
	}
	else if (y>=10.0)
	{
		y /= 10;
		eo++;
	}

	e = y[0]; // 1st digit
	*p++ = e+'0';
	*p++ = '.';

	for (int i=0; i<=16*n; i++) // totally 16n+2 significant digits
	{
		y = 10*(y - double(e));
		e = y[0];
		*p++ = e+'0';
	}

	int i=-1;
	for (; i>-16*n-1; i--) // fix possible incorrect digits
	if (p[i]<'0')
	{
		p[i]+=10;
		p[i-1]--;
	}

	if (p[i]<'0')
	{
		p[i]+=10;
		p[i-2]--;
	}

	if ((flags & 8)==0) // ignore trailing zeros?
	{
		while (*(p-1)=='0') --p;
		if (*(p-1)=='.') p++;
	}

	if (flags & 2) *p++ = 'E';
	else *p++ = 'e';

	if (eo>0 && (flags & 1)) *p++ = '+';
	if (eo<0)
	{
		*p++ = '-';
		eo = -eo;
	}

	e = eo/100;
	if (e)
	{
		*p++ = e+'0';
		eo -= 100*e;
	}

	e = eo/10;
	if (e)
	{
		*p++ = e+'0';
		eo -= 10*e;
	}
	*p++ = eo+'0';
	*p = 0;

	return p-c;
}
