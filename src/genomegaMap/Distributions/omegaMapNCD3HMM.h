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
 *  omegaMapNCD3HMM.h
 *  gcat
 *
 *  Created by Daniel Wilson on 18/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _OMEGAMAP_NCD3_HMM_H_
#define _OMEGAMAP_NCD3_HMM_H_
#include <DAG/Distribution.h>
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>
#include <omegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class omegaMapNCD3HMM : public Distribution {
private:
	// Indicates change to parameter
	bool _anc_changed;
	// Indicates changes to parameter: 1 for each species
	vector<bool> _mut_changed, _sel_changed, _phylo_changed, _p_changed;
	// Storage for species-, RV- and site- specific likelihoods
	// I.e. indexed by species (0 pop, 1 pop, 2 pop, 0 branch, 1 branch, 2 branch), then RV, then site, then codon state
	vector< vector< vector< vector<mydouble> > > > _prune_lik, _previous_prune_lik;
	// Similarly, storage for the full likelihoods
	Matrix<mydouble> _likpos, _previous_likpos;
	// Mapping of random variables to their index in _likpos and _previous_likpos
	map< const RandomVariable*, int > _rv_index;
	// Storage for likelihood calculations
	Matrix< Matrix<mydouble> > _alpha, _previous_alpha;
	
public:
	// Constructor
	omegaMapNCD3HMM(string name="", DAG* dag=0);
	// Copy constructor
	omegaMapNCD3HMM(const omegaMapNCD3HMM& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions: indexed by species in most cases
	void set_anc(Codon61SequenceRV* anc);
	void set_mut(ParentDependentRateMatrix* mut, const int sp);
	void set_sel(ContinuousVectorVariable* sel, const int sp);
	void set_phylo(NY98_TransProbMosaic* phylo, const int sp);
	void set_p(ContinuousVariable* p, const int sp);
	const Codon61SequenceRV* get_anc() const;
	const ParentDependentRateMatrix* get_mut(const int sp) const;
	const ContinuousVectorVariable* get_sel(const int sp) const;
	const NY98_TransProbMosaic* get_phylo(const int sp) const;
	const ContinuousVariable* get_p(const int sp) const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
};

class omegaMap3HMMPathSampler : public ContinuousVectorVariable, public Transformation {
private:
	string _distribution_name;		// Initially store just the name, until validate() is called
	mutable omegaMapNCD3HMM* _hmm;	// Distribution
	RandomVariable* _rv;			// RV index wrt the distribution
	int _sp;						// Species of interest
	mutable int _L;					// Length
	
public:
	// Constructor
	omegaMap3HMMPathSampler(const int species, string rv_name, string distribution_name, string name="", DAG* dag=0);
	// Copy constructor
	omegaMap3HMMPathSampler(const omegaMap3HMMPathSampler& x);
	
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

#endif // _OMEGAMAP_NCD3_HMM_H_

