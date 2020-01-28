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
 *  Codon61SequenceStationaryDistribution.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 29/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/Codon61SequenceStationaryDistribution.h>
#include <omegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string Codon61SequenceStationaryDistributionParameterNames[1] = {"pi"};

Codon61SequenceStationaryDistribution::Codon61SequenceStationaryDistribution(string name, DAG* dag) : DAGcomponent(name,dag,"Codon61SequenceStationaryDistribution"), Distribution(Codon61SequenceStationaryDistributionParameterNames,1) {
}

Codon61SequenceStationaryDistribution::Codon61SequenceStationaryDistribution(const Codon61SequenceStationaryDistribution &x) : DAGcomponent(x), Distribution(x) {
}

bool Codon61SequenceStationaryDistribution::check_random_variable_type(RandomVariable* random_variable) {
	return(dynamic_cast<Codon61SequenceRV*>(random_variable));
	return false;
}

bool Codon61SequenceStationaryDistribution::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	{	//	pi
			ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
			if(!(rv==0) && rv->length()!=61) error("Codon61SequenceStationaryDistribution::check_parameter_type(): pi must have length 61");
			return(rv);
		}
		default:
			error("Codon61SequenceStationaryDistribution::check_parameter_type(): parameter not found");
	}
	return false;
}

void Codon61SequenceStationaryDistribution::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(0,(Variable*)pi);
}

ContinuousVectorVariable const*  Codon61SequenceStationaryDistribution::get_pi() const {
	return (ContinuousVectorVariable const*)get_parameter(0);
}

mydouble Codon61SequenceStationaryDistribution::likelihood(const RandomVariable* rv, const Value* val) {
	if(val==0) error("Codon61SequenceStationaryDistribution::log_likelihood(): variable not found");
	Codon61SequenceRV& x = *((Codon61SequenceRV*)val);
	
	vector< mydouble > pi(61);
	int i;
	for(i=0;i<61;i++) {
		pi[i] = get_pi()->get_double(i);
	}
	
	mydouble ret(1.0);
	for(i=0;i<x.length();i++) {
		ret *= pi[x[i]];
	}
	return ret;
}
	
} // namespace gcat_omegaMap

