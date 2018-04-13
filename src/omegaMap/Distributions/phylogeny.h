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
 *  phylogeny.h
 *  gcat
 *
 *  Created by Daniel Wilson on 11/25/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _PHYLOGENY_H_
#define _PHYLOGENY_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class Phylogeny : public Distribution {
private:
	// Mapping of random variables to their index in _likpos and _previous_likpos
	map< const RandomVariable*, int > _rv_index;
	
public:
	// Constructor
	Phylogeny(string name="", DAG* dag=0);
	// Copy constructor
	Phylogeny(const Phylogeny& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions: indexed by species
	void set_phylo(NY98_TransProbMosaic* phylo);
	const NY98_TransProbMosaic* get_phylo() const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
};
	
} // namespace gcat_omegaMap

#endif // _PHYLOGENY_H_
