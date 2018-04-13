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
 *  NY98.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _NY98_TRANSFORMATION_H_
#define _NY98_TRANSFORMATION_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <Variables/Matrix.h>
#include <DAG/Transformation.h>
#include <matrix.h>
#include <vector.h>

using namespace gcat;

namespace gcat_omegaMap {
	
using myutils::Matrix;
using myutils::Vector;

class NY98_TransProb : public MatrixVariable, public Transformation {
private:
	mutable bool _thetaT_changed, _kappa_changed, _gamma_changed, _pi_changed;
	// Storage for calculation of eigenvectors and eigenvalues
	mutable Matrix<double> _Eigenvec, _previous_Eigenvec;
	mutable Vector<double> _Eigenval, _previous_Eigenval;
	// Storage for transition probability matrix
	mutable Matrix<double> _P, _previous_P;
public:
	// Constructor
	NY98_TransProb(string name="", DAG* dag=0);
	// Copy constructor
	NY98_TransProb(const NY98_TransProb& x);

	// Implementation of virtual functions inherited from base classes
	int nrows() const;
	int ncols() const;
	double get_double(const int i, const int j) const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_thetaT(ContinuousVariable* thetaT);
	void set_kappa(ContinuousVariable* kappa);
	void set_gamma(ContinuousVariable* gamma);
	void set_pi(ContinuousVectorVariable* pi);
	const ContinuousVariable* get_thetaT() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousVariable* get_gamma() const;
	const ContinuousVectorVariable* get_pi() const;
	
	// Overload this function inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
};
	
} // namespace gcat_omegaMap

#endif // _NY98_TRANSFORMATION_H_
