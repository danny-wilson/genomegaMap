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
 *  ParentDependentRateMatrix.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/14/09.
 *
 */
#include <genomegaMap/Transformations/ParentDependentRateMatrix.h>
#include <genomegaMap/Utilities/mutation.h>
#include <limits>

using namespace gcat;

namespace genomegaMap {
	
const string ParentDependentRateMatrixParameterNames[3] = {"theta","kappa","pi"};

ParentDependentRateMatrix::ParentDependentRateMatrix(string name, DAG* dag) : DAGcomponent(name,dag,"ParentDependentRateMatrix"), Transformation(ParentDependentRateMatrixParameterNames,3), _theta_changed(true), _kappa_changed(true), _pi_changed(true) {
}

ParentDependentRateMatrix::ParentDependentRateMatrix(const ParentDependentRateMatrix& x) : DAGcomponent(x), Transformation(x), _theta_changed(x._theta_changed), _kappa_changed(x._kappa_changed), _pi_changed(x._pi_changed), _Eigenvec(x._Eigenvec), _previous_Eigenvec(x._previous_Eigenvec), _Eigenval(x._Eigenval), _previous_Eigenval(x._previous_Eigenval), _P(x._P), _previous_P(x._previous_P) {
//	_Eigenval[0] = x._Eigenval[0];
//	_Eigenval[1] = x._Eigenval[1];
//	_Eigenvec[0] = x._Eigenvec[0];
//	_Eigenvec[1] = x._Eigenvec[1];
//	_P[0] = x._P[0];
//	_P[1] = x._P[1];
}

int ParentDependentRateMatrix::nrows() const {
	return 61;
}

int ParentDependentRateMatrix::ncols() const {
	return 61;
}

void ParentDependentRateMatrix::recalculate() const {
#pragma omp critical
	if(_kappa_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const vector<double>& pi = get_pi()->get_doubles();
		// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
		NY98_61::diagonalize_symmetric(_Eigenvec,_Eigenval,kappa,1.0,pi);
	}
#pragma omp critical
	if(_theta_changed || _kappa_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double theta = get_theta()->get_double();
		Vector<double> pi(61);
		int i,j;
		for(i=0;i<61;i++) pi[i] = get_pi()->get_double(i);
		// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
		// Therefore theta is defined as twice the expected rate under neutrality
		NY98_61::build_Rt(_P,_Eigenvec,_Eigenval,pi,theta*61.0/6.0/(2.0+kappa),1.0);
		// Sweep through to build final approximate parent-independent mutation rate matrix
		for(i=0;i<61;i++) {
			double rii = _P[i][i];
			for(j=0;j<61;j++) {
				if(j!=i) _P[i][j]/=rii;
			}
			_P[i][i] = 1.0e-6; //numeric_limits<double>::min();
		}		
	}
	_theta_changed = _kappa_changed = _pi_changed = false;
	_recalculate = false;
}

double ParentDependentRateMatrix::get_double(const int i, const int j) const {
//	Bug found 6th November 2009
//	if(_kappa_changed || _kappa_changed|| _pi_changed) {
//	if(_theta_changed || _kappa_changed || _pi_changed) {
	if(_recalculate) {
		recalculate();
	}
	return _P[i][j];
}

bool ParentDependentRateMatrix::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// theta
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("ParentDependentRateMatrix::check_parameter_type(): parameter not found");
	}
	return false;
}

void ParentDependentRateMatrix::set_theta(ContinuousVariable* theta) {
	set_parameter(0,(Variable*)theta);
}

void ParentDependentRateMatrix::set_kappa(ContinuousVariable* kappa) {
	set_parameter(1,(Variable*)kappa);
}

void ParentDependentRateMatrix::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(2,(Variable*)pi);
}

const ContinuousVariable* ParentDependentRateMatrix::get_theta() const {
	return (const ContinuousVariable*)get_parameter(0);
}

const ContinuousVariable* ParentDependentRateMatrix::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousVectorVariable* ParentDependentRateMatrix::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(2);
}

void ParentDependentRateMatrix::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_theta_changed = _kappa_changed = _pi_changed = false;
	}
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_recalculate = true;
	}
	if(v==(const Value*)get_theta()) {
		if(sgl==Variable::_SET) {
			_theta_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_theta_changed = true;
			_previous_P = _P;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_P = _previous_P;
		}
	}
	if(v==(const Value*)get_kappa()) {
		if(sgl==Variable::_SET) {
			_kappa_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_kappa_changed = true;
			_previous_Eigenvec = _Eigenvec;
			_previous_Eigenval = _Eigenval;
			_previous_P = _P;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_Eigenvec = _previous_Eigenvec;
			_Eigenval = _previous_Eigenval;
			_P = _previous_P;
		}
	}
	if(v==(const Value*)get_pi()) {
		if(sgl==Variable::_SET) {
			_pi_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_pi_changed = true;
			_previous_Eigenvec = _Eigenvec;
			_previous_Eigenval = _Eigenval;
			_previous_P = _P;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_Eigenvec = _previous_Eigenvec;
			_Eigenval = _previous_Eigenval;
			_P = _previous_P;
		}
	}
	Transformation::receive_signal_from_parent(v,sgl);
}
	
} // namespace genomegaMap
