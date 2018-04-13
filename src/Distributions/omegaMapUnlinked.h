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
 *  omegaMap_unlinked.h
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _OMEGAMAP_UNLINKED_H_
#define _OMEGAMAP_UNLINKED_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_PDRM.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class omegaMapUnlinked : public Distribution {
private:
	// Indicates changes to parameter
	bool _mut_changed;
	// Storage for RV- and site- specific likelihoods
	Matrix<mydouble> _likpos, _previous_likpos;
	// Mapping of random variables to their index in _likpos and _previous_likpos
	map< const RandomVariable*, int > _rv_index;	
	
public:
	// Constructor
	omegaMapUnlinked(string name="", DAG* dag=0);
	// Copy constructor
	omegaMapUnlinked(const omegaMapUnlinked& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_mut(NY98_ParentDependentRateMatrix* mut);
	const NY98_ParentDependentRateMatrix* get_mut() const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
};
	
} // namespace gcat_omegaMap

#endif // _OMEGAMAP_UNLINKED_H_
