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
 *  LogLogNormalQuantile.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 1/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Transformations/LogNormalQuantile.h>
#include <gsl/gsl_cdf.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string LogNormalQuantileTransformParameterNames[3] = {"mean","sd","quantile"};

LogNormalQuantileTransform::LogNormalQuantileTransform(string name, DAG* dag) : DAGcomponent(name,dag,"LogNormalQuantileTransform"), Transformation(LogNormalQuantileTransformParameterNames,3), _prev_q(0.5), _prev_z(0) {
}

LogNormalQuantileTransform::LogNormalQuantileTransform(const LogNormalQuantileTransform& x) : DAGcomponent(x), Transformation(x) {
}

double LogNormalQuantileTransform::get_double() const {
	double mean = get_mean()->get_double();
	double sd = get_sd()->get_double();
	double quantile = get_quantile()->get_double();
	
	if(!sd>0) throw BadValueException(to_Value(),"Standard deviation must be positive");
	if(quantile!=_prev_q) {
		if(!(quantile>0 & quantile<1)) throw BadValueException(to_Value(),"Quantile out of range");
		_prev_q = quantile;
		_prev_z = standard_normal_quantile_function(_prev_q);
	}
	return exp(mean+sd*_prev_z);
}

bool LogNormalQuantileTransform::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	// mean
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1:	// sd
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:	// quantile
			return(dynamic_cast<ContinuousVariable*>(parameter));
		default:
			error("LogNormalQuantileTransform::check_parameter_type(): parameter not found");
	}
	return false;
}

void LogNormalQuantileTransform::set_mean(ContinuousVariable* mean) {
	set_parameter(0,(Variable*)mean);
}

void LogNormalQuantileTransform::set_sd(ContinuousVariable* sd) {
	set_parameter(1,(Variable*)sd);
}

void LogNormalQuantileTransform::set_quantile(ContinuousVariable* quantile) {
	set_parameter(2,(Variable*)quantile);
}

ContinuousVariable const* LogNormalQuantileTransform::get_mean() const {
	return (ContinuousVariable const*)get_parameter(0);
}

ContinuousVariable const* LogNormalQuantileTransform::get_sd() const {
	return (ContinuousVariable const*)get_parameter(1);
}

ContinuousVariable const* LogNormalQuantileTransform::get_quantile() const {
	return (ContinuousVariable const*)get_parameter(2);
}

double LogNormalQuantileTransform::standard_normal_quantile_function(const double q) const {
	return gsl_cdf_ugaussian_Pinv(q);
}
	
} // namespace gcat_omegaMap
