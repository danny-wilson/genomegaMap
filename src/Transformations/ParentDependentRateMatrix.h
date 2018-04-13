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
 *  ParentDependentRateMatrix.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/14/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 *	Transformation providing a parent-dependent approximation to
 *	the HKY85 codon model for use in a parent-independent likelihood.
 *
 *		theta	Twice the population mutation rate per generation per codon
 *		kappa	The transition:transversion ratio
 *		pi		Codon usage frequencies
 *
 */
#ifndef _PARENT_DEPENDENT_RATE_MATRIX_H_
#define _PARENT_DEPENDENT_RATE_MATRIX_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <Variables/Matrix.h>
#include <DAG/Transformation.h>
#include <matrix.h>
#include <vector.h>

using myutils::Matrix;
using myutils::Vector;

using namespace gcat;

namespace gcat_omegaMap {
	
class ParentDependentRateMatrix : public MatrixVariable, public Transformation {
private:
	mutable bool _theta_changed, _kappa_changed, _pi_changed;
	// Storage for calculate of eigenvectors and eigenvalues
//	Matrix<double> _Eigenvec[2];
//	Vector<double> _Eigenval[2];
	mutable Matrix<double> _Eigenvec, _previous_Eigenvec;
	mutable Vector<double> _Eigenval, _previous_Eigenval;
	// Storage for rate matrix
//	Matrix<double> _P[2];
	mutable Matrix<double> _P, _previous_P;
	// Indicator of which _Eigenval and _Eigenvec matrix is current
//	int eix;
	// Indicator of which _P matrix is current
//	int pix;
public:
	// Constructor
	ParentDependentRateMatrix(string name="", DAG* dag=0);
	// Constructor
	ParentDependentRateMatrix(const ParentDependentRateMatrix& x);
	
	// Implementation of virtual functions inherited from base classes
	int nrows() const;
	int ncols() const;
	double get_double(const int i, const int j) const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_theta(ContinuousVariable* theta);
	void set_kappa(ContinuousVariable* kappa);
	void set_pi(ContinuousVectorVariable* pi);
	const ContinuousVariable* get_theta() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousVectorVariable* get_pi() const;
	
	// Overload this function inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);

//protected:
	void recalculate() const;
	//void recalc() const;
};
	
} // namespace gcat_omegaMap

#endif // _PARENT_DEPENDENT_RATE_MATRIX_H_
