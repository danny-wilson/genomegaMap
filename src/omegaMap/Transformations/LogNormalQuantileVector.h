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
 *  LogNormalQuantileVector.h
 *  gcat
 *
 *  Created by Daniel Wilson on 24/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _LOG_NORMAL_QUANTILE_VECTOR_TRANSFORM_H_
#define _LOG_NORMAL_QUANTILE_VECTOR_TRANSFORM_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <DAG/Transformation.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class LogNormalQuantileVectorTransform : public ContinuousVectorVariable, public Transformation {
private:
	int _n;
	mutable bool _mean_changed, _sd_changed, _quantile_changed;
	mutable vector<double> _x, _x_prev;
	mutable vector<bool> _bad, _bad_prev;
	vector<bool> _has_changed;
public:
	// Constructor
	LogNormalQuantileVectorTransform(const int n, string name="", DAG* dag=0);
	// Copy constructor
	LogNormalQuantileVectorTransform(const LogNormalQuantileVectorTransform& x);
	
	// Implementation of virtual functions inherited from base classes
	// Get length of the variable
	int length() const;
	// Get value at position i
	double get_double(const int i) const;
	// Get vector of values
	vector<double> get_doubles() const;
	// Has the value changed at position i?
	bool has_changed(const int i) const;
	// Has the value changed at each position?
	vector<bool> has_changed() const;
	// Type-checking for parameter(s)
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_mean(ContinuousVariable* mean);
	void set_sd(ContinuousVariable* sd);
	void set_quantile(ContinuousVectorVariable* quantile);
	ContinuousVariable const* get_mean() const;	
	ContinuousVariable const* get_sd() const;	
	ContinuousVectorVariable const* get_quantile() const;
	
	// Overload method inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	void recalculate() const;
	
protected:
	double standard_normal_quantile_function(const double q) const;
};
	
} // namespace gcat_omegaMap

#endif // _LOG_NORMAL_QUANTILE_VECTOR_TRANSFORM_H_



