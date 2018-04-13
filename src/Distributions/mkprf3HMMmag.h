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
 *  mkprf3HMMmagmag.h
 *  gcat
 *
 *  Created by Daniel Wilson on 24/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	Like mkprf3HMM but parameterized by the magnitude of gamma, and
 *	q, the probability that a site is under positive selection.
 *
 *	While variation in gamma is modelled using the HMM, each site
 *	has an independent probability of positive or negative selection,
 *	and the resulting selection coefficient is +gamma or -gamma
 *	respectively.
 *
 */
#ifndef _MKPRF_3_HMM_MAG_H_
#define _MKPRF_3_HMM_MAG_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>
#include <ostream>
#include <random.h>

using std::ostream;
using std::endl;
using myutils::Random;

using namespace gcat;

// External global variable (see DAG.cpp)
namespace gcat {
	extern Random _ran;
}

namespace gcat_omegaMap {
	
// Forward declaration
class Codon61GroupCount;

class mkprf3HMMmag : public Distribution {
private:
	// Indicates changes to parameter: 1 for each species
	vector<bool> _thetaS_changed, _thetaR_changed, _gamma_changed, _tau_changed, _p_changed, _q_changed;
	// Mapping of random variables to their index in _likpos and _previous_likpos, etc
	map< const RandomVariable*, int > _rv_index;
	// Storage for the likelihoods for each RV, site and species
	vector< vector< vector< mydouble > > > _likpos, _previous_likpos;
	// Coding of the data, for each RV, site and species
	vector< vector< vector<int> > > _site;
	// Indicates whether the data are valid, for each RV and site
	vector< vector<bool> > _bad;
	
	// Storage for likelihood calculations
	Matrix< Matrix<mydouble> > _alpha, _previous_alpha;
	// Storage for probability of negative selection
	Matrix< Matrix<mydouble> > _pneg, _previous_pneg;
	
public:
	// Constructor
	mkprf3HMMmag(string name="", DAG* dag=0);
	// Copy constructor
	mkprf3HMMmag(const mkprf3HMMmag& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions: indexed by species
	void set_thetaS(ContinuousVariable* thetaS, const int sp);
	void set_thetaR(ContinuousVariable* thetaR, const int sp);
	void set_gamma(ContinuousVectorVariable* gamma, const int sp);
	void set_tau(ContinuousVariable* tau, const int sp);
	void set_p(ContinuousVariable* p, const int sp);
	void set_q(ContinuousVariable* q, const int sp);
	const ContinuousVariable* get_thetaS(const int sp) const;
	const ContinuousVariable* get_thetaR(const int sp) const;
	const ContinuousVectorVariable* get_gamma(const int sp) const;
	const ContinuousVariable* get_tau(const int sp) const;
	const ContinuousVariable* get_p(const int sp) const;
	const ContinuousVariable* get_q(const int sp) const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
	// Functions used by the likelihood
	double F(const double n, const double theta, const double gamma) const;
	double G(const double n, const double theta, const double gamma) const;
	
	void summarize(ostream& out, Codon61GroupCount* rv, const int ix);
	
	vector<double> path_sampler(const RandomVariable* const rv, const int sp);
	
private:
	void codify_sites(Codon61GroupCount* rv, const int ix);
	double hypergeometric1F1(const double a, const double b, const double c) const;
};
	
} // namespace gcat_omegaMap

#endif // _MKPRF_3_HMM_MAG_H_


