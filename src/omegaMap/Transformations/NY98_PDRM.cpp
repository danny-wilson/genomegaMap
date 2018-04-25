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
 *  NY98_PDRM.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/NY98_PDRM.h>
#include <omegaMap/Utilities/mutation.h>
#include <limits>

using namespace gcat;

namespace gcat_omegaMap {
	
const string NY98_ParentDependentRateMatrixParameterNames[4] = {"theta","kappa","omega","pi"};

NY98_ParentDependentRateMatrix::NY98_ParentDependentRateMatrix(const int n, string name, DAG* dag) : DAGcomponent(name,dag,"NY98_ParentDependentRateMatrix"), Transformation(NY98_ParentDependentRateMatrixParameterNames,4), _n(n), _theta_changed(true), _kappa_changed(true), _omega_changed(true), _pi_changed(true), _has_changed(Vector<bool>(n,true)), _Eigenvec(n), _previous_Eigenvec(n), _Eigenval(n), _previous_Eigenval(n), _meanrate(n), _previous_meanrate(n), _P(n), _previous_P(n), _is_block_start(n,false), _previous_is_block_start(n,false) {
}

NY98_ParentDependentRateMatrix::NY98_ParentDependentRateMatrix(const NY98_ParentDependentRateMatrix& x) : DAGcomponent(x), Transformation(x), _n(x._n), _theta_changed(x._theta_changed), _kappa_changed(x._kappa_changed), _omega_changed(x._omega_changed), _pi_changed(x._pi_changed), _has_changed(x._has_changed), _Eigenvec(x._Eigenvec), _previous_Eigenvec(x._previous_Eigenvec), _Eigenval(x._Eigenval), _previous_Eigenval(x._previous_Eigenval), _meanrate(x._meanrate), _previous_meanrate(x._previous_meanrate), _P(x._P), _previous_P(x._previous_P), _is_block_start(x._is_block_start), _previous_is_block_start(x._previous_is_block_start) {
}

int NY98_ParentDependentRateMatrix::length() const {
	return _n;
}

int NY98_ParentDependentRateMatrix::nrows() const {
	return 61;
}

int NY98_ParentDependentRateMatrix::ncols() const {
	return 61;
}

void NY98_ParentDependentRateMatrix::recalculate() const {
	const ContinuousMosaicVariable& omega = *get_omega();
#pragma omp critical
	if(_kappa_changed || _omega_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const vector<double>& pi = get_pi()->get_doubles();
		int pos;
		for(pos=0;pos<_n;pos++) {
			if(omega.is_block_start(pos)) {
				_is_block_start[pos] = true;
				if((_kappa_changed || _pi_changed || omega.has_changed(pos))) {
					const double omega_pos = omega.get_double(pos);
					// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
					NY98_61::diagonalize_symmetric(_Eigenvec[pos],_meanrate[pos],_Eigenval[pos],kappa,omega_pos,pi);
                    //cout << "NY98_61::diagonalize_symmetric() with kappa " << kappa;
                    //cout << " omega " << omega_pos << " pi " << pi[0] << " ";
                    //cout << pi[1] << " " << pi[2] << " " << pi[3] << " ... returned ";
                    //cout << _meanrate[pos] << endl;
				}
			}
			else {
				_is_block_start[pos] = false;
			}
		}
	}
#pragma omp critical
	if(_theta_changed || _kappa_changed || _omega_changed || _pi_changed) {
		const double kappa = get_kappa()->get_double();
		const double theta = get_theta()->get_double();
		Vector<double> pi(61);
		int ii;
		for(ii=0;ii<61;ii++) pi[ii] = get_pi()->get_double(ii);
		int pos;
		for(pos=0;pos<_n;pos++) {
			if(omega.is_block_start(pos) && (_theta_changed || _kappa_changed || _pi_changed || omega.has_changed(pos))) {
				// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
				// Therefore thetaT is defined as twice the expected number of mutations under neutrality
				//NY98_61::build_Rt(_P[pos],_Eigenvec[pos],_Eigenval[pos],pi,theta*61.0/6.0/(2.0+kappa),1.0);
				// No longer using this equal codon usage approximation, but using the computed mean rate:
				NY98_61::build_Rt(_P[pos],_Eigenvec[pos],_Eigenval[pos],pi,(theta/2.0)/_meanrate[pos],1.0);
                //cout << "NY98_61::build_Rt() called with scale " << (theta/2.0)/_meanrate[pos] << " lambda 1.0 ";
                //cout << " pi " << pi[0] << " ";
                //cout << pi[1] << " " << pi[2] << " " << pi[3] << " ...";
                //cout << _meanrate[pos] << endl;
                // Sweep through to build final approximate parent-independent mutation rate matrix
				int i,j;
				for(i=0;i<61;i++) {
					double rii = _P[pos][i][i];
					for(j=0;j<61;j++) {
						if(j!=i) _P[pos][i][j]/=rii;
					}
					_P[pos][i][i] = 1.0e-6; //numeric_limits<double>::min();
				}
                //cout << "P[0][0..4] = " << _P[pos][0][0] << " " << _P[pos][0][1] << " ";
                //cout << _P[pos][0][2] << " " << _P[pos][0][3] << " " << _P[pos][0][4];
                //cout << endl;
			}
		}
	}
	_theta_changed = _kappa_changed = _omega_changed = _pi_changed = false;
	_recalculate = false;
}

double NY98_ParentDependentRateMatrix::get_double(const int i, const int j, const int k) const {
	if(_recalculate) {
		recalculate();
	}
	int pos = get_omega()->block_start(i);
	return _P[pos][j][k];
}

bool NY98_ParentDependentRateMatrix::has_changed(const int pos) const {
	return _has_changed[pos];
}

bool NY98_ParentDependentRateMatrix::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// theta
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// omega
			return(dynamic_cast<ContinuousMosaicVariable*>(parameter));
		case 3:	// pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("NY98_ParentDependentRateMatrix::check_parameter_type(): parameter not found");
	}
	return false;
}

void NY98_ParentDependentRateMatrix::set_theta(ContinuousVariable* theta) {
	set_parameter(0,(Variable*)theta);
}

void NY98_ParentDependentRateMatrix::set_kappa(ContinuousVariable* kappa) {
	set_parameter(1,(Variable*)kappa);
}

void NY98_ParentDependentRateMatrix::set_omega(ContinuousMosaicVariable* omega) {
	set_parameter(2,(Variable*)omega);
}

void NY98_ParentDependentRateMatrix::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(3,(Variable*)pi);
}

const ContinuousVariable* NY98_ParentDependentRateMatrix::get_theta() const {
	return (const ContinuousVariable*)get_parameter(0);
}

const ContinuousVariable* NY98_ParentDependentRateMatrix::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousMosaicVariable* NY98_ParentDependentRateMatrix::get_omega() const {
	return (const ContinuousMosaicVariable*)get_parameter(2);
}

const ContinuousVectorVariable* NY98_ParentDependentRateMatrix::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(3);
}

void NY98_ParentDependentRateMatrix::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_theta_changed = _kappa_changed = _omega_changed = _pi_changed = false;
	}
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_recalculate = true;
	}
	int pos;
	if(v==(const Value*)get_theta()) {
		if(sgl==Variable::_SET) {
			_theta_changed = true;
			_has_changed = Vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_theta_changed = true;
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
					_previous_meanrate[pos] = _meanrate[pos];
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
					_meanrate[pos] = _previous_meanrate[pos];
					_P[pos] = _previous_P[pos];
				}
			}
			_is_block_start = _previous_is_block_start;
			_has_changed = Vector<bool>(_n,false);
		}
	}
	if(v==(const Value*)get_omega()) {
		if(sgl==Variable::_SET) {
			_omega_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(get_omega()->has_changed(pos)) _has_changed[pos]=true;
			}
		}
		else if(sgl==Variable::_PROPOSE) {
			_omega_changed = true;
			for(pos=0;pos<_n;pos++) {
				if(_is_block_start[pos]) {
					_previous_Eigenvec[pos] = _Eigenvec[pos];
					_previous_Eigenval[pos] = _Eigenval[pos];
					_previous_meanrate[pos] = _meanrate[pos];
					_previous_P[pos] = _P[pos];
				}
			}
			_previous_is_block_start = _is_block_start;
			for(pos=0;pos<_n;pos++) {
				if(get_omega()->has_changed(pos)) _has_changed[pos]=true;
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
					_meanrate[pos] = _previous_meanrate[pos];
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
					_previous_meanrate[pos] = _meanrate[pos];
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
					_meanrate[pos] = _previous_meanrate[pos];
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
