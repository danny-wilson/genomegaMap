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
 *  NY98_PDRM_Single.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/NY98_PDRM_Single.h>
#include <omegaMap/Utilities/mutation.h>
#include <limits>

using namespace gcat;

namespace gcat_omegaMap {
	
const string NY98_ParentDependentRateMatrix_SingleParameterNames[4] = {"theta","kappa","omega","pi"};

NY98_ParentDependentRateMatrix_Single::NY98_ParentDependentRateMatrix_Single(string name, DAG* dag) : DAGcomponent(name,dag,"NY98_ParentDependentRateMatrix_Single"), Transformation(NY98_ParentDependentRateMatrix_SingleParameterNames,4), _theta_changed(true), _kappa_changed(true), _omega_changed(true), _pi_changed(true) {
}

NY98_ParentDependentRateMatrix_Single::NY98_ParentDependentRateMatrix_Single(const NY98_ParentDependentRateMatrix_Single& x) : DAGcomponent(x), Transformation(x), _theta_changed(x._theta_changed), _kappa_changed(x._kappa_changed), _omega_changed(x._omega_changed), _pi_changed(x._pi_changed), _Eigenvec(x._Eigenvec), _previous_Eigenvec(x._previous_Eigenvec), _Eigenval(x._Eigenval), _previous_Eigenval(x._previous_Eigenval), _P(x._P), _previous_P(x._previous_P) {
}

int NY98_ParentDependentRateMatrix_Single::nrows() const {
	return 61;
}

int NY98_ParentDependentRateMatrix_Single::ncols() const {
	return 61;
}

void NY98_ParentDependentRateMatrix_Single::recalculate() const {
#pragma omp critical
	if(_kappa_changed || _omega_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double omega = get_omega()->get_double();
		const vector<double>& pi = get_pi()->get_doubles();
		// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
		NY98_61::diagonalize_symmetric(_Eigenvec,_Eigenval,kappa,omega,pi);
	}
#pragma omp critical
	if(_theta_changed || _kappa_changed || _omega_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double omega = get_omega()->get_double();
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
	_theta_changed = _kappa_changed = _omega_changed = _pi_changed = false;
	_recalculate = false;
}

double NY98_ParentDependentRateMatrix_Single::get_double(const int i, const int j) const {
	if(_recalculate) {
		recalculate();
	}
	return _P[i][j];
}

bool NY98_ParentDependentRateMatrix_Single::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// theta
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// omega
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 3:	// pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("NY98_ParentDependentRateMatrix_Single::check_parameter_type(): parameter not found");
	}
	return false;
}

void NY98_ParentDependentRateMatrix_Single::set_theta(ContinuousVariable* theta) {
	set_parameter(0,(Variable*)theta);
}

void NY98_ParentDependentRateMatrix_Single::set_kappa(ContinuousVariable* kappa) {
	set_parameter(1,(Variable*)kappa);
}

void NY98_ParentDependentRateMatrix_Single::set_omega(ContinuousVariable* omega) {
	set_parameter(2,(Variable*)omega);
}

void NY98_ParentDependentRateMatrix_Single::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(3,(Variable*)pi);
}

const ContinuousVariable* NY98_ParentDependentRateMatrix_Single::get_theta() const {
	return (const ContinuousVariable*)get_parameter(0);
}

const ContinuousVariable* NY98_ParentDependentRateMatrix_Single::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousVariable* NY98_ParentDependentRateMatrix_Single::get_omega() const {
	return (const ContinuousVariable*)get_parameter(2);
}

const ContinuousVectorVariable* NY98_ParentDependentRateMatrix_Single::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(3);
}

void NY98_ParentDependentRateMatrix_Single::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_theta_changed = _kappa_changed = _omega_changed = _pi_changed = false;
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
	if(v==(const Value*)get_omega()) {
		if(sgl==Variable::_SET) {
			_omega_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_omega_changed = true;
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
	
} // namespace gcat_omegaMap
