/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the omegaMap library.
 *
 *  The omegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The omegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the omegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/*
 *  NY98_TransProbMosaicMosaic.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _NY98_TRANSPROB_MOSAIC_H_
#define _NY98_TRANSPROB_MOSAIC_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <Variables/ContinuousMosaic.h>
#include <Variables/MatrixVector.h>
#include <DAG/Transformation.h>
#include <matrix.h>
#include <vector.h>

using myutils::Matrix;
using myutils::Vector;

using namespace gcat;

namespace gcat_omegaMap {
	
class NY98_TransProbMosaic : public MatrixVectorVariable, public Transformation {
private:
	// Length of mosaic
	int _n;
	// Flags to indicate signals received from parent variables
	mutable bool _thetaT_changed, _kappa_changed, _gamma_changed, _pi_changed;
	// Storage for has_changed() method
	Vector<bool> _has_changed;
	// Storage for calculation of eigenvectors and eigenvalues
	mutable Vector< Matrix<double> > _Eigenvec, _previous_Eigenvec;
	mutable Vector< Vector<double> > _Eigenval, _previous_Eigenval;
	// Storage for transition probability matrix
	mutable Vector< Matrix<double> > _P, _previous_P;
	// Storage for the start of blocks, for efficiency
	mutable Vector< bool > _is_block_start, _previous_is_block_start;
public:
	// Constructor
	NY98_TransProbMosaic(const int n, string name="", DAG* dag=0);
	// Copy constructor
	NY98_TransProbMosaic(const NY98_TransProbMosaic& x);
	
	// Implementation of virtual functions inherited from base classes
	int length() const;
	int nrows() const;
	int ncols() const;
	double get_double(const int i, const int j, const int k) const;
	bool has_changed(const int i) const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_thetaT(ContinuousVariable* thetaT);
	void set_kappa(ContinuousVariable* kappa);
	void set_gamma(ContinuousMosaicVariable* gamma);
	void set_pi(ContinuousVectorVariable* pi);
	const ContinuousVariable* get_thetaT() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousMosaicVariable* get_gamma() const;
	const ContinuousVectorVariable* get_pi() const;
	
	// Overload this function inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
//protected:
	// (Not really const, but use of mutables justified on grounds of efficiency gains)
	void recalculate() const;
	//void recalc() const;
};
	
} // namespace gcat_omegaMap

#endif // _NY98_TRANSPROB_MOSAIC_H_

