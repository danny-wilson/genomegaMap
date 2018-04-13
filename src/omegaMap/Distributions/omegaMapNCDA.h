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
 *  omegaMapNCDA.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _OMEGAMAP_NCDA_H_
#define _OMEGAMAP_NCDA_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class omegaMapNCDA : public Distribution {
private:
	// Indicates changes to parameter
	bool _mut_changed, _sel_changed;

public:
	// Constructor
	omegaMapNCDA(string name="", DAG* dag=0);
	// Copy constructor
	omegaMapNCDA(const omegaMapNCDA& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	void set_mut(ParentDependentRateMatrix* mut);
	void set_sel(ContinuousVectorVariable* sel);
	const ParentDependentRateMatrix* get_mut() const;
	const ContinuousVectorVariable* get_sel() const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
};
	
} // namespace gcat_omegaMap

#endif // _OMEGAMAP_NCDA_H_
