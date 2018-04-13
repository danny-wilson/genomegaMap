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
 *  NY98_TransProbMosaicMosaic.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Utilities/mutation.h>
#include <limits>

using namespace gcat;

namespace gcat_omegaMap {
	
const string NY98_TransProbMosaicParameterNames[4] = {"thetaT","kappa","gamma","pi"};

NY98_TransProbMosaic::NY98_TransProbMosaic(const int n, string name, DAG* dag) : DAGcomponent(name,dag,"NY98_TransProbMosaic"), Transformation(NY98_TransProbMosaicParameterNames,4), _n(n), _thetaT_changed(true), _kappa_changed(true), _gamma_changed(true), _pi_changed(true), _has_changed(Vector<bool>(n,true)),
_Eigenvec(n), _previous_Eigenvec(n), _Eigenval(n), _previous_Eigenval(n), _P(n), _previous_P(n), _is_block_start(n,false), _previous_is_block_start(n,false) {
}

NY98_TransProbMosaic::NY98_TransProbMosaic(const NY98_TransProbMosaic& x) : DAGcomponent(x), Transformation(x), _n(x._n), _thetaT_changed(x._thetaT_changed), _kappa_changed(x._kappa_changed), _gamma_changed(x._gamma_changed), _pi_changed(x._pi_changed), _has_changed(x._has_changed), 
_Eigenvec(x._Eigenvec), _previous_Eigenvec(x._previous_Eigenvec), _Eigenval(x._Eigenval), _previous_Eigenval(x._previous_Eigenval), _P(x._P), _previous_P(x._previous_P), _is_block_start(x._is_block_start), _previous_is_block_start(x._previous_is_block_start) {
}

int NY98_TransProbMosaic::length() const {
	return _n;
}

int NY98_TransProbMosaic::nrows() const {
	return 61;
}

int NY98_TransProbMosaic::ncols() const {
	return 61;
}

void NY98_TransProbMosaic::recalculate() const {
	const ContinuousMosaicVariable& gamma = *get_gamma();
#pragma omp critical
	if(_kappa_changed || _gamma_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const vector<double>& pi = get_pi()->get_doubles();
		int pos;
		for(pos=0;pos<_n;pos++) {
			if(gamma.is_block_start(pos)) {
				_is_block_start[pos] = true;
				if((_kappa_changed || _pi_changed || gamma.has_changed(pos))) {
					const double gamma_pos = gamma.get_double(pos);
					const double omega = (fabs(gamma_pos)<1e-3) ? 1.0 : gamma_pos/(1-exp(-gamma_pos));
					// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
					NY98_61::diagonalize_symmetric(_Eigenvec[pos],_Eigenval[pos],kappa,omega,pi);
				}
			}
			else {
				_is_block_start[pos] = false;
			}
		}
	}
#pragma omp critical
	if(_thetaT_changed || _kappa_changed || _gamma_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double thetaT = get_thetaT()->get_double();
		Vector<double> pi(61);
		int ii;
		for(ii=0;ii<61;ii++) pi[ii] = get_pi()->get_double(ii);
		int pos;
		for(pos=0;pos<_n;pos++) {
			if(gamma.is_block_start(pos) && (_thetaT_changed || _kappa_changed || _pi_changed || gamma.has_changed(pos))) {
				// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
				// Therefore thetaT is defined as twice the expected number of mutations under neutrality
				NY98_61::build_Pt(_P[pos],_Eigenvec[pos],_Eigenval[pos],pi,thetaT*61.0/6.0/(2.0+kappa));
			}
		}
	}
	_thetaT_changed = _kappa_changed = _gamma_changed = _pi_changed = false;
	_recalculate = false;
}

// Position i
double NY98_TransProbMosaic::get_double(const int i, const int j, const int k) const {
//	Bug found 6th November 2009
//	if(_kappa_changed || _gamma_changed || _pi_changed || _pi_changed) {
//	if(_thetaT_changed || _kappa_changed || _gamma_changed || _pi_changed) {
	if(_recalculate) {
		recalculate();
	}
	int pos = get_gamma()->block_start(i);
	return _P[pos][j][k];
}

bool NY98_TransProbMosaic::has_changed(const int pos) const {
	return _has_changed[pos];
}

bool NY98_TransProbMosaic::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// thetaT
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// gamma
			return(dynamic_cast<ContinuousMosaicVariable*>(parameter));
		case 3:	// pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("NY98_TransProbMosaic::check_parameter_type(): parameter not found");
	}
	return false;
}

void NY98_TransProbMosaic::set_thetaT(ContinuousVariable* thetaT) {
	set_parameter(0,(Variable*)thetaT);
}

void NY98_TransProbMosaic::set_kappa(ContinuousVariable* kappa) {
	set_parameter(1,(Variable*)kappa);
}

void NY98_TransProbMosaic::set_gamma(ContinuousMosaicVariable* gamma) {
	set_parameter(2,(Variable*)gamma);
}

void NY98_TransProbMosaic::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(3,(Variable*)pi);
}

const ContinuousVariable* NY98_TransProbMosaic::get_thetaT() const {
	return (const ContinuousVariable*)get_parameter(0);
}

const ContinuousVariable* NY98_TransProbMosaic::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousMosaicVariable* NY98_TransProbMosaic::get_gamma() const {
	return (const ContinuousMosaicVariable*)get_parameter(2);
}

const ContinuousVectorVariable* NY98_TransProbMosaic::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(3);
}

void NY98_TransProbMosaic::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_thetaT_changed = _kappa_changed = _gamma_changed = _pi_changed = false;
	}
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_recalculate = true;
	}
	int pos;
	if(v==(const Value*)get_thetaT()) {
		if(sgl==Variable::_SET) {
			_thetaT_changed = true;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_thetaT_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(_is_block_start[pos]) {
					_previous_P[pos] = _P[pos];
				}
			}
			_previous_is_block_start = _is_block_start;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_ACCEPT) {
			_has_changed = Vector<bool>(_n,false);
		}
		else if(sgl==Variable::_REVERT) {
			for(pos=0;pos<_n;pos++) {
				if(_previous_is_block_start[pos]) {
					_P[pos] = _previous_P[pos];
				}
			}
			_is_block_start = _previous_is_block_start;
			_has_changed = Vector<bool>(_n,false);
		}
	}
	if(v==(const Value*)get_kappa()) {
		if(sgl==Variable::_SET) {
			_kappa_changed = true;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_kappa_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(_is_block_start[pos]) {
					_previous_Eigenvec[pos] = _Eigenvec[pos];
					_previous_Eigenval[pos] = _Eigenval[pos];
					_previous_P[pos] = _P[pos];
				}
			}
			_previous_is_block_start = _is_block_start;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_ACCEPT) {
			_has_changed = Vector<bool>(_n,false);
		}
		else if(sgl==Variable::_REVERT) {
			for(pos=0;pos<_n;pos++) {
				if(_previous_is_block_start[pos]) {
					_Eigenvec[pos] = _previous_Eigenvec[pos];
					_Eigenval[pos] = _previous_Eigenval[pos];
					_P[pos] = _previous_P[pos];
				}
			}
			_is_block_start = _previous_is_block_start;
			_has_changed = Vector<bool>(_n,false);
		}
	}
	if(v==(const Value*)get_gamma()) {
		if(sgl==Variable::_SET) {
			_gamma_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(get_gamma()->has_changed(pos)) _has_changed[pos]=true;
			}
		}
		else if(sgl==Variable::_PROPOSE) {
			_gamma_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(_is_block_start[pos]) {
					_previous_Eigenvec[pos] = _Eigenvec[pos];
					_previous_Eigenval[pos] = _Eigenval[pos];
					_previous_P[pos] = _P[pos];
				}
			}
			_previous_is_block_start = _is_block_start;
			for(pos=0;pos<_n;pos++) {
				if(get_gamma()->has_changed(pos)) _has_changed[pos]=true;
			}
		}
		else if(sgl==Variable::_ACCEPT) {
			_has_changed = Vector<bool>(_n,false);
		}
		else if(sgl==Variable::_REVERT) {
			for(pos=0;pos<_n;pos++) {
				if(_previous_is_block_start[pos]) {
					_Eigenvec[pos] = _previous_Eigenvec[pos];
					_Eigenval[pos] = _previous_Eigenval[pos];
					_P[pos] = _previous_P[pos];
				}
			}
			_is_block_start = _previous_is_block_start;
			_has_changed = Vector<bool>(_n,false);
		}
	}
	if(v==(const Value*)get_pi()) {
		if(sgl==Variable::_SET) {
			_pi_changed = true;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_pi_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(_is_block_start[pos]) {
					_previous_Eigenvec[pos] = _Eigenvec[pos];
					_previous_Eigenval[pos] = _Eigenval[pos];
					_previous_P[pos] = _P[pos];
				}
			}
			_previous_is_block_start = _is_block_start;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_ACCEPT) {
			_has_changed = Vector<bool>(_n,false);
		}
		else if(sgl==Variable::_REVERT) {
			for(pos=0;pos<_n;pos++) {
				if(_previous_is_block_start[pos]) {
					_Eigenvec[pos] = _previous_Eigenvec[pos];
					_Eigenval[pos] = _previous_Eigenval[pos];
					_P[pos] = _previous_P[pos];
				}
			}
			_is_block_start = _previous_is_block_start;
			_has_changed = Vector<bool>(_n,false);
		}
	}
	Transformation::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
