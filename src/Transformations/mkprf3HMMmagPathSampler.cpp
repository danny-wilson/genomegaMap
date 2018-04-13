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
 *  mkprf3HMMmagmagPathSampler.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 24/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <DAG/DAG.h>
#include <omegaMap/Transformations/mkprf3HMMmagPathSampler.h>
#include <gsl/gsl_cdf.h>
#include <DAG/RandomVariable.h>
#include <omegaMap/RandomVariables/Codon61GroupCount.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string mkprf3HMMmagPathSamplerParameterNames[0];

mkprf3HMMmagPathSampler::mkprf3HMMmagPathSampler(const int species, string rv_name, string distribution_name, string name, DAG* dag) : DAGcomponent(name,dag,"mkprf3HMMmagPathSampler"), Transformation(mkprf3HMMmagPathSamplerParameterNames,0), _distribution_name(distribution_name), _sp(species) {
	_rv = getDAG()->get_random_variable(rv_name);
}

mkprf3HMMmagPathSampler::mkprf3HMMmagPathSampler(const mkprf3HMMmagPathSampler& x) : DAGcomponent(x), Transformation(x), _distribution_name(x._distribution_name), _hmm(x._hmm), _rv(x._rv), _sp(x._sp), _L(x._L) {
}

string mkprf3HMMmagPathSampler::validate() const {
	// For starters, do not allow this to parameterize anything else
	if(n_child_distributions()+n_child_transformations()>0) {
		string errTxt = "mkprf3HMMmagPathSampler: object " + name() + " may not parameterize other objects";
		return errTxt;
	}
	// Secondly, finish identifying the distribution to which it relates (Distributions
	// are initialized after Transformations, so this could not be done in the constructor)
	_hmm = dynamic_cast<mkprf3HMMmag*>(getDAG()->get_distribution(_distribution_name));
	if(!_hmm) {
		string errTxt = "mkprf3HMMmagPathSampler: object " + name() + " cannot find mkprf3HMMmag object " + _distribution_name;
		return errTxt;
	}
	const Codon61GroupCount* ct = dynamic_cast<const Codon61GroupCount*>(_rv);
	//	const Codon61GroupCount* ct = dynamic_cast<const Codon61GroupCount*>(_rv->to_Value());
	if(!ct) {
		string errTxt = "mkprf3HMMmagPathSampler: " + _rv->name() + " is not of type Codon61GroupCount";
		return errTxt;
	}
	_L = ct->length();
	return "";
}

int mkprf3HMMmagPathSampler::length() const {
	return _L;
}

// NB: Repeatedly calling get_double() gives marginal draws, whereas get_doubles() draws a complete path
// For that reason, the print() method is overwritten
double mkprf3HMMmagPathSampler::get_double(const int i) const {
	if(i<0) error("mkprf3HMMmagPathSampler::get_double(): index cannot be negative");
	vector<double> x = get_doubles();
	if(i>=_L) error("mkprf3HMMmagPathSampler::get_double(): index too large");
	return x[i];
}

vector<double> mkprf3HMMmagPathSampler::get_doubles() const {
	return _hmm->path_sampler(_rv,_sp);
}

bool mkprf3HMMmagPathSampler::has_changed(const int i) const {
	return true;
}

vector<bool> mkprf3HMMmagPathSampler::has_changed() const {
	return vector<bool>(length(),true);
}

bool mkprf3HMMmagPathSampler::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		default:
			error("mkprf3HMMmagPathSampler::check_parameter_type(): parameter not found");
	}
	return false;
}

void mkprf3HMMmagPathSampler::print(ostream& out, string sep) {
	int i;
	vector<double> x = get_doubles();
	for(i=0;i<length();i++) {
		if(i>0) out << sep;
		out << x[i];
	}
}
	
} // namespace gcat_omegaMap
