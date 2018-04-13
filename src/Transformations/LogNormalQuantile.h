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
 *  LogNormalQuantile.h
 *  gcat
 *
 *  Created by Daniel Wilson on 1/6/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _LOG_NORMAL_QUANTILE_TRANSFORM_H_
#define _LOG_NORMAL_QUANTILE_TRANSFORM_H_
#include <Variables/Continuous.h>
#include <DAG/Transformation.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class LogNormalQuantileTransform : public ContinuousVariable, public Transformation {
private:
	mutable double _prev_q, _prev_z;
public:
	// Constructor
	LogNormalQuantileTransform(string name="", DAG* dag=0);
	// Copy constructor
	LogNormalQuantileTransform(const LogNormalQuantileTransform& x);
	
	// Implementation of virtual functions inherited from base classes
	double get_double() const;
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_mean(ContinuousVariable* mean);
	void set_sd(ContinuousVariable* sd);
	void set_quantile(ContinuousVariable* quantile);
	ContinuousVariable const* get_mean() const;	
	ContinuousVariable const* get_sd() const;	
	ContinuousVariable const* get_quantile() const;
	
protected:
	double standard_normal_quantile_function(const double q) const;
};
	
} // namespace gcat_omegaMap

#endif // _LOG_NORMAL_QUANTILE_TRANSFORM_H_

