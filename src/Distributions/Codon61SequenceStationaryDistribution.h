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
 *  Codon61SequenceStationaryDistribution.h
 *  gcat
 *
 *  Created by Daniel Wilson on 29/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _CODON_61_SEQUENCE_STATIONARY_DISTRIBUTION_H_
#define _CODON_61_SEQUENCE_STATIONARY_DISTRIBUTION_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class Codon61SequenceStationaryDistribution : public Distribution {
protected:
public:
	// Constructor
	Codon61SequenceStationaryDistribution(string name="", DAG* dag=0);
	// Copy constructor
	Codon61SequenceStationaryDistribution(const Codon61SequenceStationaryDistribution& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	bool check_parameter_type(const int i, Variable* parameter);
	void set_pi(ContinuousVectorVariable* pi);
	ContinuousVectorVariable const* get_pi() const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
};
	
} // namespace gcat_omegaMap

#endif // _CODON_61_SEQUENCE_STATIONARY_DISTRIBUTION_H_

