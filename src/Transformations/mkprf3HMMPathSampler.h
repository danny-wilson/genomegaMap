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
 *  mkprf3HMMPathSampler.h
 *  gcat
 *
 *  Created by Daniel Wilson on 17/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _MKPRF3_HMM_PATH_SAMPLER_H_
#define _MKPRF3_HMM_PATH_SAMPLER_H_
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <DAG/Transformation.h>
#include <omegaMap/Distributions/mkprf3HMM.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class mkprf3HMMPathSampler : public ContinuousVectorVariable, public Transformation {
private:
	string _distribution_name;		// Initially store just the name, until validate() is called
	mutable mkprf3HMM* _hmm;		// Distribution
	RandomVariable* _rv;			// RV index wrt the distribution
	int _sp;						// Species of interest
	mutable int _L;					// Length

public:
	// Constructor
	mkprf3HMMPathSampler(const int species, string rv_name, string distribution_name, string name="", DAG* dag=0);
	// Copy constructor
	mkprf3HMMPathSampler(const mkprf3HMMPathSampler& x);
	
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
	
	// Overload method inherited from Transformation
	//void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	//void recalculate() const;

	// Overload method inherited from ContinuousVectorVariable
	void print(ostream& out, string sep);

protected:
	// Overload method inherited from Component
	string validate() const;
};
	
} // namespace gcat_omegaMap

#endif // _MKPRF3_HMM_PATH_SAMPLER_H_


