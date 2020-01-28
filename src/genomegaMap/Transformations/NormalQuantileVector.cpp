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
 *  NormalQuantileVector.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 16/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/NormalQuantileVector.h>
#include <gsl/gsl_cdf.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string NormalQuantileVectorTransformParameterNames[3] = {"mean","sd","quantile"};

NormalQuantileVectorTransform::NormalQuantileVectorTransform(const int n, string name, DAG* dag) : DAGcomponent(name,dag,"NormalQuantileVectorTransform"), Transformation(NormalQuantileVectorTransformParameterNames,3), _n(n), _mean_changed(true), _sd_changed(true), _quantile_changed(true), _x(n), _x_prev(n), _has_changed(n,true), _bad(n), _bad_prev(n) {
}

NormalQuantileVectorTransform::NormalQuantileVectorTransform(const NormalQuantileVectorTransform& x) : DAGcomponent(x), Transformation(x), _n(x._n), _mean_changed(x._mean_changed), _sd_changed(x._sd_changed), _quantile_changed(x._quantile_changed), _x(x._x), _x_prev(x._x_prev), _has_changed(x._has_changed), _bad(x._bad), _bad_prev(x._bad_prev)  {
}

int NormalQuantileVectorTransform::length() const {
	return _n;
}

double NormalQuantileVectorTransform::get_double(const int i) const {
	if(_recalculate) recalculate();
	// Only throw if the value is requested
	if(_bad[i]) throw BadValueException(to_Value(),"Standard deviation or quantile out of range");
	return _x[i];
}

void NormalQuantileVectorTransform::recalculate() const {
	//	if(_mean_changed || _sd_changed || _quantile_changed) {
	// Recalculate
	double mean = get_mean()->get_double();
	double sd = get_sd()->get_double();
	if(!sd>0) _bad = vector<bool>(_n,true);
	else {
		int pos;
		for(pos=0;pos<_n;pos++) {
			double quantile = get_quantile()->get_double(pos);
			if(!(quantile>0 & quantile<1)) _bad[pos] = true;
			else {
				_x[pos] = mean+sd*standard_normal_quantile_function(quantile);
				_bad[pos] = false;
			}
		}
	}
	_mean_changed = _sd_changed = _quantile_changed = false;
	//	}
	_recalculate = false;
}

vector<double> NormalQuantileVectorTransform::get_doubles() const {
	vector<double> ret(_n);
	int i;
	for(i=0;i<_n;i++) ret[i] = get_double(i);
	return ret;
}

bool NormalQuantileVectorTransform::has_changed(const int i) const {
	return _has_changed[i];
}

vector<bool> NormalQuantileVectorTransform::has_changed() const {
	return _has_changed;
}

bool NormalQuantileVectorTransform::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// mean
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// sd
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// quantile
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("NormalQuantileVectorTransform::check_parameter_type(): parameter not found");
	}
	return false;
}

void NormalQuantileVectorTransform::set_mean(ContinuousVariable* mean) {
	set_parameter(0,(Variable*)mean);
}

void NormalQuantileVectorTransform::set_sd(ContinuousVariable* sd) {
	set_parameter(1,(Variable*)sd);
}

void NormalQuantileVectorTransform::set_quantile(ContinuousVectorVariable* quantile) {
	set_parameter(2,(Variable*)quantile);
}

ContinuousVariable const* NormalQuantileVectorTransform::get_mean() const {
	return (ContinuousVariable const*)get_parameter(0);
}

ContinuousVariable const* NormalQuantileVectorTransform::get_sd() const {
	return (ContinuousVariable const*)get_parameter(1);
}

ContinuousVectorVariable const* NormalQuantileVectorTransform::get_quantile() const {
	return (ContinuousVectorVariable const*)get_parameter(2);
}

double NormalQuantileVectorTransform::standard_normal_quantile_function(const double q) const {
	return gsl_cdf_ugaussian_Pinv(q);
}

void NormalQuantileVectorTransform::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_mean_changed = _sd_changed = _quantile_changed = false;
		_has_changed = vector<bool>(_n,false);
	}
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_recalculate = true;
	}
	if(v==(const Value*)get_mean()) {
		if(sgl==Variable::_SET) {
			_mean_changed = true;
			_has_changed = vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_mean_changed = true;
			_has_changed = vector<bool>(_n,true);
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
	else if(v==(const Value*)get_sd()) {
		if(sgl==Variable::_SET) {
			_sd_changed = true;
			_has_changed = vector<bool>(_n,true);
		}
		else if(sgl==Variable::_PROPOSE) {
			_sd_changed = true;
			_has_changed = vector<bool>(_n,true);
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
	else if(v==(const Value*)get_quantile()) {
		if(sgl==Variable::_SET) {
			_quantile_changed = true;
			_has_changed = get_quantile()->has_changed();
		}
		else if(sgl==Variable::_PROPOSE) {
			_quantile_changed = true;
			_has_changed = get_quantile()->has_changed();
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
	
} // namespace gcat_omegaMap

