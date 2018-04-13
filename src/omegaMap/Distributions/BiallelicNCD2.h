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
 *  BiallelicNCD2.h
 *  gcat
 *
 *  Created by Daniel Wilson on 23/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _BIALLELIC_NCD2_H_
#define _BIALLELIC_NCD2_H_
#include <DAG/Distribution.h>
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <ostream>

using std::ostream;
using std::endl;

using namespace gcat;

namespace gcat_omegaMap {
	
// Forward declaration
//class BiallelicCodonGroupCount;

class BiallelicNCD2 : public Distribution {
private:
	// Indicates changes to parameter: 1 for each species
	vector<bool> _thetaS_changed, _thetaR_changed, _gamma_changed, _tau_changed;
	// Mapping of random variables to their index in _likpos and _previous_likpos, etc
	map< const RandomVariable*, int > _rv_index;
	// Storage for the likelihoods for each RV, site and species
	vector< vector< vector< mydouble > > > _likpos, _previous_likpos;
	
public:
	// Constructor
	BiallelicNCD2(string name="", DAG* dag=0);
	// Copy constructor
	BiallelicNCD2(const BiallelicNCD2& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions: indexed by species
	void set_thetaS(ContinuousVariable* thetaS, const int sp);
	void set_thetaR(ContinuousVariable* thetaR, const int sp);
	void set_gamma(ContinuousVectorVariable* gamma, const int sp);
	void set_tau(ContinuousVariable* tau, const int sp);
	const ContinuousVariable* get_thetaS(const int sp) const;
	const ContinuousVariable* get_thetaR(const int sp) const;
	const ContinuousVectorVariable* get_gamma(const int sp) const;
	const ContinuousVariable* get_tau(const int sp) const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);

	// Site-specific likelihood
	mydouble site_likelihood(const double n[], const double x[], const double theta[], const double tau[], const double gamma[]) const;
	// Ancestrally conditioned sampling formula
	double acsf(const double n, const double x, const int anc, const double theta, const double gamma) const;
	// Ancestrally conditioned sampling formula, pooling segregated classes into one
	double acsf_pool(const double n, const double x, const int anc, const double theta, const double gamma) const;
	// Probability a site is segregating in a population
	double pseg(const double n, const double theta, const double gamma) const;
	// Probability a site is fixed in a population
	double pfxd(const double n, bool fixed_for_anc, const double theta, const double gamma) const;
	
private:
	double hypergeometric1F1(const double a, const double b, const double c) const;
};
	
} // namespace gcat_omegaMap

#endif // _BIALLELIC_NCD2_H_
