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
 *  omegaMapNCDHMMHybridHybrid.h
 *  gcat
 *
 *  Created by Daniel Wilson on 02/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _OMEGAMAP_NCD_HMM_HYBRID_H_
#define _OMEGAMAP_NCD_HMM_HYBRID_H_
#include <DAG/Distribution.h>
#include <Variables/Continuous.h>
#include <Variables/ContinuousVector.h>
#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>
#include <omegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

// External global variable (see DAG.cpp)
namespace gcat {
	extern Random _ran;
}

namespace gcat_omegaMap {

class omegaMapNCDHMMHybrid : public Distribution {
private:
	// Sequence length
	int _L;
	// Number of classes of gamma (the latent variable) in sliding window model
	int _ng1;
	// Number of classes of gamma (the latent variable) in sitewise model
	int _ng2;
	// Indicates change to parameters. Only allow receive_signal_from_parent to changes these.
	bool _anc_changed, _theta_changed, _kappa_changed, _gamma_changed, _T_changed, _p_changed, _pi_changed, _gamma_wt_changed;
	// Indicates changes are required to internal variables. Set by receive_signal_from_parent but can be unset by other functions.
	bool _update_phylo_rate_matrix, _update_popgen_rate_matrix, _update_pclik, _update_pemiss, _update_alphabeta;
	// Mapping of random variables to their index in _likpos, etc
	map< const RandomVariable*, int > _rv_index;
	// Storage for RV-specific variables
	// Position at which to calculate the HMM likelihood. Indexed by rv.
	vector< int > _likpos, _previous_likpos;
	// Emission probabilities. Indexed by rv, site and latent variable.
	vector< Matrix< mydouble> > _pemiss, _previous_pemiss;
	// Forward/backward probabilities. Indexed by rv, site and latent variable.
	vector< Matrix< mydouble > > _alphabeta, _previous_alphabeta;
	// Ancestor-dependent (i.e. phylogenetically conditioned) likelihoods. Indexed by rv, site, latent variable and ancestral allele.
	vector< vector< Matrix< mydouble > > > _pclik, _previous_pclik;
	// Storage for non-RV-specific variables ************* NEED TO ADD TO BLAH BLAH
	// Phylogenetic transition probabilities. Indexed by latent variable, ancestral and descendant allele.
	Matrix<double> _phylo_temp;
	vector< Matrix< mydouble > > _phylo, _previous_phylo;
	vector< Matrix<double> > _phylo_Eigenvec, _previous_phylo_Eigenvec;
	vector< Vector<double> > _phylo_Eigenval, _previous_phylo_Eigenval;
	// Parent dependent rate matrix. Indexed by ancestral and descendant allele.
	Matrix<double> _pdrm, _previous_pdrm;
	Matrix<double> _pdrm_Eigenvec, _previous_pdrm_Eigenvec;
	Vector<double> _pdrm_Eigenval, _previous_pdrm_Eigenval;
	
public:
	// Constructor
	omegaMapNCDHMMHybrid(const int L, const int ng1, const int ng2, string name="", DAG* dag=0);
	// Copy constructor
	//omegaMapNCDHMMHybrid(const omegaMapNCDHMMHybrid& x);
	// Implementation of virtual function inherited from base class Distribution
	bool check_random_variable_type(RandomVariable* random_variable);
	// Implementation of virtual function inherited from base class DependentVariable
	virtual bool check_parameter_type(const int i, Variable* parameter);
	
	// Convenience functions
	int seqlen() const;
	int ngamma1() const;
	int ngamma2() const;
	void set_anc(Codon61SequenceRV* anc);
	void set_theta(ContinuousVariable* theta);
	void set_kappa(ContinuousVariable* kappa);
	void set_gamma1(ContinuousVectorVariable* gamma1);
	void set_gamma2(ContinuousVectorVariable* gamma2);
	void set_T(ContinuousVariable* T);
	void set_p(ContinuousVariable* p);
	void set_pi(ContinuousVectorVariable* pi);
	void set_gamma1_wt(ContinuousVectorVariable* gamma1_wt);
	void set_gamma2_wt(ContinuousVectorVariable* gamma2_wt);
	const Codon61SequenceRV* get_anc() const;
	const ContinuousVariable* get_theta() const;
	const ContinuousVariable* get_kappa() const;
	const ContinuousVectorVariable* get_gamma1() const;
	const ContinuousVectorVariable* get_gamma2() const;
	const ContinuousVariable* get_T() const;
	const ContinuousVariable* get_p() const;
	const ContinuousVectorVariable* get_pi() const;
	const ContinuousVectorVariable* get_gamma1_wt() const;
	const ContinuousVectorVariable* get_gamma2_wt() const;
	
	// Compute likelihood
	mydouble likelihood(const RandomVariable* rv, const Value* val);
	
	// Overload initialization functions to allow readying for efficient likelihood calculation
	void add_random_variable(RandomVariable* random_variable);
	void remove_random_variable(RandomVariable* random_variable);
	// Overload signalling function inherited from Distribution
	void receive_signal_from_parent(const Value* v, const Variable::Signal sgl);
	
	vector<double> path_sampler(const RandomVariable* const rv);
	
private:
	// Various intermediate steps in the likelihood calculation, modularised for computational efficiency
	void update_phylo_rate_matrix();
	void update_popgen_rate_matrix();
	void update_pclik(const Value* val, const int ix);
	void update_pemiss(const int ix);
	void update_alphabeta(const int ix);
	// Recentre likpos
	int recentre() const;
	// When _likpos == i, (0<= i < L), then alphabeta[,0:i] contain forward probabilities for
	// sites 0:i and alphabeta[,(i+1):L] contain backward probabilities for sites i:(L-1).
	// Update forward probabilities for sites left:_likpos (0<= left < L)
	void forward(const int ix, const int left);
	// Update backward probabilities for sites _likpos:right (0<= right < L)
	void backward(const int ix, const int right);
	// HMM likelihood
	mydouble HMMlik(const int ix);
	// Part of posterior decoding, required for the hybrid model
	double final_draw(const int pos, const int draw, const double U_in, const mydouble m, const vector<double> &gamma1, const vector<double> &gamma2, const vector<mydouble> &pg, vector<mydouble> qdraw, const Matrix<mydouble> &pemiss);
};

class omegaMapNCDHMMHybridPathSampler : public ContinuousVectorVariable, public Transformation {
private:
	string _distribution_name;		// Initially store just the name, until validate() is called
	mutable omegaMapNCDHMMHybrid* _hmm;	// Distribution
	RandomVariable* _rv;			// RV index wrt the distribution
	mutable int _L;					// Length
	
public:
	// Constructor
	omegaMapNCDHMMHybridPathSampler(string rv_name, string distribution_name, string name="", DAG* dag=0);
	// Copy constructor
	omegaMapNCDHMMHybridPathSampler(const omegaMapNCDHMMHybridPathSampler& x);
	
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

#endif // _OMEGAMAP_NCD_HMM_HYBRID_H_
