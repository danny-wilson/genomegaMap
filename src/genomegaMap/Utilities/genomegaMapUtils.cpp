/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the genomegaMap library.
 *
 *  The genomegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The genomegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the genomegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/*
 *  genomegaMapUtils.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *
 */
#include <genomegaMap/Utilities/genomegaMapUtils.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_hyperg.h>
#include <cstdio>
#include <math.h>
#include <myerror.h>

namespace genomegaMapUtils {
	const double SQRT2PI = 2.5066282746310;
	const char* oneLetterCodes = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-";

	double LambertW(const double z) {
		return gsl_sf_lambert_W0(z);
		// K M Briggs
		const int dbgW=0;
		int i; 
		const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
		double p,e,t,w;
		if (dbgW) fprintf(stderr,"genomegaMapUtils::LambertW(): z=%g\n",z);
		if (z<-em1) return log(-1.);
		if (0.0==z) return 0.0;
		if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
			double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
			return 
			-1.0
			+2.331643981597124203363536062168*r
			-1.812187885639363490240191647568*q
			+1.936631114492359755363277457668*r*q
			-2.353551201881614516821543561516*q2
			+3.066858901050631912893148922704*r*q2
			-4.175335600258177138854984177460*q3
			+5.858023729874774148815053846119*r*q3
			-8.401032217523977370984161688514*q3*q;  // error approx 1e-16
		}
		/* initial approx for iteration... */
		if (z<1.0) { /* series near 0 */
			p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
			w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
		} else 
			w=log(z); /* asymptotic */
		if (z>3.0) w-=log(w); /* useful? */
		const int MAXIT = 10;
		for (i=0; i<MAXIT; i++) { /* Halley iteration */
			e=exp(w); 
			t=w*e-z;
			p=w+1.0;
			t/=e*p-0.5*(p+1.0)*t/p; 
			w-=t;
			if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
		}
		return w;
		/* should never get here */
		myutils::error("genomegaMapUtils::LambertW(): No convergence at z");
		return 0.;
	}
	
	double LambertW1(const double z) {
		return gsl_sf_lambert_Wm1(z);
		// K M Briggs
		int i; 
		const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
		double p,e,t,w,l1,l2;
		if (z<-em1 || z>=0.0 || isinf(z) || isnan(z)) { 
			myutils::error("genomegaMapUtils::LambertW1(): bad argument %g");
		}
		/* initial approx for iteration... */
		if (z<-1e-6) { /* series about -1/e */
			p=-sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
			w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
		} else { /* asymptotic near zero */
			l1=log(-z);
			l2=log(-l1);
			w=l1-l2+l2/l1;
		}
		if (fabs(p)<1e-4) return w;
		for (i=0; i<10; i++) { /* Halley iteration */
			e=exp(w); 
			t=w*e-z;
			p=w+1.0;
			t/=e*p-0.5*(p+1.0)*t/p; 
			w-=t;
			if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
		}
		/* should never get here */
		return -z;
		myutils::error("genomegaMapUtils::LambertW1(): No convergence at z=%g");
		return 0.;
	}
	
	double hypergeometric1F1(const double a, const double b, const double c) {
		double ret;
//#pragma omp critical
		ret = gsl_sf_hyperg_1F1(a,b,c);
		return ret;
	}
}

