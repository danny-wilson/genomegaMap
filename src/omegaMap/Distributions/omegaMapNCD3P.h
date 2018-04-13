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
 *  omegaMapNCD3Pparsimony.h
 *  gcat
 *
 *  Created by Daniel Wilson on 1/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	An approximation to omegaMapNCD3P in which (originally) ancestral states are inferred parsimoniously.
 *	Now the code has been adapted to restrict the permissible ancestors to those observed in at least
 *	one sample. The permissible ancestors are calculated independently for each codon.
 *
 */
#ifndef _OMEGAMAP_NCD3_PARSIMONY_H_
#define _OMEGAMAP_NCD3_PARSIMONY_H_
#include <DAG/Distribution.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>

using namespace gcat;

namespace gcat_omegaMap {
	
// Forward declaration
class Codon61GroupCount;

class omegaMapNCD3P : public Distribution {
private:
	// Indicates changes to parameter: 1 for each species
	vector<bool> _mut_changed, _sel_changed, _phylo_changed;
	// Storage for species-, RV- and site- specific likelihoods
	// I.e. indexed by species (0 pop, 1 pop, 2 pop, 0 branch, 1 branch, 2 branch), then RV, then site, then codon state
	vector< vector< vector< vector<mydouble> > > > _prune_lik, _previous_prune_lik;
	// Similarly, storage for the full likelihoods
	Matrix<mydouble> _likpos, _previous_likpos;
	// Mapping of random variables to their index in _likpos and _previous_likpos
	map< const RandomVariable*, int > _rv_index;
	// Indicates the ancestral codon for every internal node
	//vector< vector< vector<int> > > _ancestral_codon;
	// Indicates whether a codon is a permissible ancestor, which requires that it was observed in at least one sample
	vector< vector< vector<bool> > > _permissible_ancestor;
	
public:
	// Constructor
	omegaMapNCD3P(string name="", DAG* dag=0);
	// Copy constructor
	omegaMapNCD3P(const omegaMapNCD3P& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions: indexed by species
	void set_mut(ParentDependentRateMatrix* mut, const int sp);
	void set_sel(ContinuousVectorVariable* sel, const int sp);
	void set_phylo(NY98_TransProbMosaic* phylo, const int sp);
	const ParentDependentRateMatrix* get_mut(const int sp) const;
	const ContinuousVectorVariable* get_sel(const int sp) const;
	const NY98_TransProbMosaic* get_phylo(const int sp) const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
private:
	// Compute the _ancestral_codon matrix by parsimony
	//void infer_ancestral_codons(Codon61GroupCount* rv, const int ix);
	// Identify the permissible ancestors for each site
	void permissible_ancestors(Codon61GroupCount* rv, const int ix);
};
	
} // namespace gcat_omegaMap

#endif // _OMEGAMAP_NCD3_PARSIMONY_H_


