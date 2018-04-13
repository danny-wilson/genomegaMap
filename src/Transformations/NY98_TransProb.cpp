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
 *  NY98_TransProb.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/NY98_TransProb.h>
#include <omegaMap/Utilities/mutation.h>
#include <limits>

using namespace gcat;

namespace gcat_omegaMap {
	
const string NY98_TransProbParameterNames[4] = {"thetaT","kappa","gamma","pi"};

NY98_TransProb::NY98_TransProb(string name, DAG* dag) : DAGcomponent(name,dag,"NY98_TransProb"), Transformation(NY98_TransProbParameterNames,4), _thetaT_changed(true), _kappa_changed(true), _gamma_changed(true), _pi_changed(true) {
}

NY98_TransProb::NY98_TransProb(const NY98_TransProb& x) : DAGcomponent(x), Transformation(x), _thetaT_changed(x._thetaT_changed), _kappa_changed(x._kappa_changed), _gamma_changed(x._gamma_changed), _pi_changed(x._pi_changed), _Eigenvec(x._Eigenvec), _previous_Eigenvec(x._previous_Eigenvec), _Eigenval(x._Eigenval), _previous_Eigenval(x._previous_Eigenval), _P(x._P), _previous_P(x._previous_P) {
}

int NY98_TransProb::nrows() const {
	return 61;
}

int NY98_TransProb::ncols() const {
	return 61;
}

double NY98_TransProb::get_double(const int i, const int j) const {
	if(_kappa_changed || _gamma_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double gamma = get_gamma()->get_double();
		const double omega = (fabs(gamma)<1e-3) ? 1.0 : gamma/(1-exp(-gamma));
		const vector<double>& pi = get_pi()->get_doubles();
		// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
		NY98_61::diagonalize_symmetric(_Eigenvec,_Eigenval,kappa,omega,pi);
	}
	if(_thetaT_changed || _kappa_changed || _gamma_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double thetaT = get_thetaT()->get_double();
		Vector<double> pi(61);
		int i;
		for(i=0;i<61;i++) pi[i] = get_pi()->get_double(i);
		// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
		// Therefore thetaT is defined as twice the expected number of mutations under neutrality
		NY98_61::build_Pt(_P,_Eigenvec,_Eigenval,pi,thetaT*61.0/6.0/(2.0+kappa));
	}
	_thetaT_changed = _kappa_changed = _gamma_changed = _pi_changed = false;
	return _P[i][j];
}

bool NY98_TransProb::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// thetaT
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// gamma
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 3:	// pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("NY98_TransProb::check_parameter_type(): parameter not found");
	}
	return false;
}

void NY98_TransProb::set_thetaT(ContinuousVariable* thetaT) {
	set_parameter(0,(Variable*)thetaT);
}

void NY98_TransProb::set_kappa(ContinuousVariable* kappa) {
	set_parameter(1,(Variable*)kappa);
}

void NY98_TransProb::set_gamma(ContinuousVariable* gamma) {
	set_parameter(2,(Variable*)gamma);
}

void NY98_TransProb::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(3,(Variable*)pi);
}

const ContinuousVariable* NY98_TransProb::get_thetaT() const {
	return (const ContinuousVariable*)get_parameter(0);
}

const ContinuousVariable* NY98_TransProb::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousVariable* NY98_TransProb::get_gamma() const {
	return (const ContinuousVariable*)get_parameter(2);
}

const ContinuousVectorVariable* NY98_TransProb::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(3);
}

void NY98_TransProb::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_thetaT_changed = _kappa_changed = _gamma_changed = _pi_changed = false;
	}
	if(v==(const Value*)get_thetaT()) {
		if(sgl==Variable::_SET) {
			_thetaT_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_thetaT_changed = true;
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
	if(v==(const Value*)get_gamma()) {
		if(sgl==Variable::_SET) {
			_gamma_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_gamma_changed = true;
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
