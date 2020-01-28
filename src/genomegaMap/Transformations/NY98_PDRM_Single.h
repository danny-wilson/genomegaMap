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
 *  NY98_PDRM_Single.h
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *
 */
#ifndef _NY98_PARENT_DEPENDENT_RATE_MATRIX_SINGLE_H_
#define _NY98_PARENT_DEPENDENT_RATE_MATRIX_SINGLE_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <Variables/Matrix.h>
#include <DAG/Transformation.h>
#include <matrix.h>
#include <vector.h>

using namespace gcat;

namespace genomegaMap {
	
using myutils::Matrix;
using myutils::Vector;

class NY98_ParentDependentRateMatrix_Single : public MatrixVariable, public Transformation {
private:
	mutable bool _theta_changed, _kappa_changed, _omega_changed, _pi_changed;
	// Storage for calculate of eigenvectors and eigenvalues
	mutable Matrix<double> _Eigenvec, _previous_Eigenvec;
	mutable Vector<double> _Eigenval, _previous_Eigenval;
	// Storage for rate matrix
	mutable Matrix<double> _P, _previous_P;
public:
	// Constructor
	NY98_ParentDependentRateMatrix_Single(string name="", DAG* dag=0);
	// Constructor
	NY98_ParentDependentRateMatrix_Single(const NY98_ParentDependentRateMatrix_Single& x);
	
	// Implementation of virtual functions inherited from base classes
	int nrows() const;
	int ncols() const;
	double get_double(const int i, const int j) const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_theta(ContinuousVariable* theta);
	void set_kappa(ContinuousVariable* kappa);
	void set_omega(ContinuousVariable* omega);
	void set_pi(ContinuousVectorVariable* pi);
	const ContinuousVariable* get_theta() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousVariable* get_omega() const;
	const ContinuousVectorVariable* get_pi() const;
	
	// Overload this function inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
	//protected:
	void recalculate() const;
};
	
} // namespace genomegaMap

#endif // _NY98_PARENT_DEPENDENT_RATE_MATRIX_SINGLE_H_
