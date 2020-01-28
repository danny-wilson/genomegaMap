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
 *  omega2gammaVector.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 03/03/2010.
 *
 */
#include <genomegaMap/Transformations/Omega2GammaVector.h>
#include <genomegaMap/Utilities/genomegaMapUtils.h>

using namespace gcat;

namespace genomegaMap {
	
const string Omega2GammaVectorTransformParameterNames[1] = {"omega"};

Omega2GammaVectorTransform::Omega2GammaVectorTransform(const int n, string name, DAG* dag) : DAGcomponent(name,dag,"Omega2GammaVectorTransform"), Transformation(Omega2GammaVectorTransformParameterNames,1), _n(n), _omega_changed(true), _x(n), _x_prev(n), _has_changed(n,true), _bad(n), _bad_prev(n) {
}

Omega2GammaVectorTransform::Omega2GammaVectorTransform(const Omega2GammaVectorTransform& x) : DAGcomponent(x), Transformation(x), _n(x._n), _omega_changed(x._omega_changed), _x(x._x), _x_prev(x._x_prev), _has_changed(x._has_changed), _bad(x._bad), _bad_prev(x._bad_prev)  {
}

int Omega2GammaVectorTransform::length() const {
	return _n;
}

double Omega2GammaVectorTransform::get_double(const int i) const {
	if(_recalculate) recalculate();
	// Only throw if the value is requested
	if(_bad[i]) throw BadValueException(to_Value(),"omega out of range");
	return _x[i];
}

void Omega2GammaVectorTransform::recalculate() const {
	int pos;
	double oprev;
	for(pos=0;pos<_n;pos++) {
		double omega = get_omega()->get_double(pos);
		if(pos>0 && omega==oprev) {
			_x[pos] = _x[pos-1];
			_bad[pos] = _bad[pos-1];
		}
		else {
			if(!(omega>0)) _bad[pos] = true;
			else {
				_x[pos] = toGamma(omega);
				_bad[pos] = false;
			}
			oprev = omega;
		}
	}
	_omega_changed = false;
	_recalculate = false;
}

vector<double> Omega2GammaVectorTransform::get_doubles() const {
	vector<double> ret(_n);
	int i;
	for(i=0;i<_n;i++) ret[i] = get_double(i);
	return ret;
}

bool Omega2GammaVectorTransform::has_changed(const int i) const {
	return _has_changed[i];
}

vector<bool> Omega2GammaVectorTransform::has_changed() const {
	return _has_changed;
}

bool Omega2GammaVectorTransform::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// omega
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("Omega2GammaVectorTransform::check_parameter_type(): parameter not found");
	}
	return false;
}

void Omega2GammaVectorTransform::set_omega(ContinuousVectorVariable* omega) {
	set_parameter(0,(Variable*)omega);
}

ContinuousVectorVariable const* Omega2GammaVectorTransform::get_omega() const {
	return (ContinuousVectorVariable const*)get_parameter(0);
}

void Omega2GammaVectorTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_omega_changed = false;
		_has_changed = vector<bool>(_n,false);
	}
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_recalculate = true;
	}
	if(v==(const Value*)get_omega()) {
		if(sgl==Variable::_SET) {
			_omega_changed = true;
			_has_changed = get_omega()->has_changed();
		}
		else if(sgl==Variable::_PROPOSE) {
			_omega_changed = true;
			_has_changed = get_omega()->has_changed();
			_x_prev = _x;
			_bad_prev = _bad;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_x = _x_prev;
			_bad = _bad_prev;
		}
	}
	// Call default implementation, which is to call Variable::send_signal_to_children(sgl)
	Transformation::receive_signal_from_parent(v,sgl);
}

double Omega2GammaVectorTransform::toGamma(const double omega) const {
	return omega + ((omega>1) ? genomegaMapUtils::LambertW(-omega*exp(-omega)) : genomegaMapUtils::LambertW1(-omega*exp(-omega)));
}
	
} // namespace genomegaMap

