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
 *  NormalQuantileMosaic.h
 *  gcat
 *
 *  Created by Daniel Wilson on 1/8/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _NORMAL_QUANTILE_MOSAIC_TRANSFORM_H_
#define _NORMAL_QUANTILE_MOSAIC_TRANSFORM_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousMosaic.h>
#include <DAG/Transformation.h>

/*	TO DO:
 
	Need to separate _has_changed from flags to indicate that recalculation is necessary.
	Need to update all elements if recalculation is necessary, because it's possible that
	only some might be queried, and then subsequently uncalculated values may be used.
	(Could this also be a problem for the similar Transformations in omegaMap?)
 
 */

using namespace gcat;

namespace gcat_omegaMap {
	
class NormalQuantileMosaicTransform : public ContinuousMosaicVariable, public Transformation {
private:
	int _n;
	mutable bool _mean_changed, _sd_changed, _quantile_changed;
	mutable vector<double> _x, _x_prev;
	mutable vector<bool> _bad, _bad_prev;
	vector<bool> _has_changed;
public:
	// Constructor
	NormalQuantileMosaicTransform(const int n, string name="", DAG* dag=0);
	// Copy constructor
	NormalQuantileMosaicTransform(const NormalQuantileMosaicTransform& x);
	
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
	// Get the number of breakpoints
	int nblocks() const;
	// Is there a left breakpoint at position i?
	bool is_block_start(const int i) const;
	// Is there a right breakpoint at position i?
	bool is_block_end(const int i) const;
	// Where is the start of the current block?
	int block_start(const int i) const;
	// Where is the end of the current block?
	int block_end(const int i) const;
	// Type-checking for parameter(s)
	bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_mean(ContinuousVariable* mean);
	void set_sd(ContinuousVariable* sd);
	void set_quantile(ContinuousMosaicVariable* quantile);
	ContinuousVariable const* get_mean() const;	
	ContinuousVariable const* get_sd() const;	
	ContinuousMosaicVariable const* get_quantile() const;
	
	// Overload method inherited from Transformation
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	void recalculate() const;
	
protected:
	double standard_normal_quantile_function(const double q) const;
};
	
} // namespace gcat_omegaMap

#endif // _NORMAL_QUANTILE_MOSAIC_TRANSFORM_H_

