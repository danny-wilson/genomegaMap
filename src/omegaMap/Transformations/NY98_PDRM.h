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
 *  NY98_PDRM.h
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	Transformation providing a parent-dependent approximation to
 *	the NY98 codon model for use in a parent-independent likelihood.
 *
 *		theta	Twice the *neutral* population mutation rate per generation per codon
 *		kappa	The transition:transversion ratio
 *		omega	The dN:dS ratio (i.e. selection modelled as mutational bias)
 *		pi		Codon usage frequencies
 *
 */
#ifndef _NY98_PARENT_DEPENDENT_RATE_MATRIX_H_
#define _NY98_PARENT_DEPENDENT_RATE_MATRIX_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <Variables/ContinuousMosaic.h>
#include <Variables/MatrixVector.h>
#include <DAG/Transformation.h>
#include <matrix.h>
#include <vector.h>

using namespace gcat;

namespace gcat_omegaMap {
	
using myutils::Matrix;
using myutils::Vector;

class NY98_ParentDependentRateMatrix : public MatrixVectorVariable, public Transformation {
private:
	// Length of mosaic
	int _n;
	// Flags to indicate signals received from parent variables
	mutable bool _theta_changed, _kappa_changed, _omega_changed, _pi_changed;
	// Storage for has_changed() method
	Vector<bool> _has_changed;
	// Storage for calculation of eigenvectors, eigenvalues and mean rates
	mutable Vector< Matrix<double> > _Eigenvec, _previous_Eigenvec;
	mutable Vector< Vector<double> > _Eigenval, _previous_Eigenval;
	mutable Vector< double > _meanrate, _previous_meanrate;
	// Storage for rate matrix
	mutable Vector< Matrix<double> > _P, _previous_P;
	// Storage for the start of blocks, for efficiency
	mutable Vector< bool > _is_block_start, _previous_is_block_start;
public:
	// Constructor
	NY98_ParentDependentRateMatrix(const int n, string name="", DAG* dag=0);
	// Constructor
	NY98_ParentDependentRateMatrix(const NY98_ParentDependentRateMatrix& x);
	
	// Implementation of virtual functions inherited from base classes
	int length() const;
	int nrows() const;
	int ncols() const;
	double get_double(const int i, const int j, const int k) const;
	bool has_changed(const int i) const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_theta(ContinuousVariable* theta);
	void set_kappa(ContinuousVariable* kappa);
	void set_omega(ContinuousMosaicVariable* omega);
	void set_pi(ContinuousVectorVariable* pi);
	const ContinuousVariable* get_theta() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousMosaicVariable* get_omega() const;
	const ContinuousVectorVariable* get_pi() const;
	
	// Overload this function inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
	//protected:
	// (Not really const, but use of mutables justified on grounds of efficiency gains)
	void recalculate() const;
};
	
} // namespace gcat_omegaMap

#endif // _NY98_PARENT_DEPENDENT_RATE_MATRIX_H_
