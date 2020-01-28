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
 *  omegaMapNCD2.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/16/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/omegaMapNCD2.h>
#include <omegaMap/RandomVariables/Codon61GroupCount.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>
//#include <omp.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string omegaMapNCD2ParameterNames[6] = {"mut0","sel0","phylo0","mut1","sel1","phylo1"};

omegaMapNCD2::omegaMapNCD2(string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapNCD2"), Distribution(omegaMapNCD2ParameterNames,6),
_mut_changed(vector<bool>(2,true)), _sel_changed(vector<bool>(2,true)), _phylo_changed(vector<bool>(2,true)),
_prune_lik(vector< vector< vector< vector<mydouble> > > >(4,vector< vector< vector<mydouble> > >(0))), _previous_prune_lik(vector< vector< vector< vector<mydouble> > > >(4,vector< vector< vector<mydouble> > >(0))) {
}

omegaMapNCD2::omegaMapNCD2(const omegaMapNCD2& x) : DAGcomponent(x), Distribution(x), 
_mut_changed(x._mut_changed), _sel_changed(x._sel_changed), _phylo_changed(x._phylo_changed),
_likpos(x._likpos), _previous_likpos(x._previous_likpos), _prune_lik(x._prune_lik), _previous_prune_lik(x._previous_prune_lik) {
}

bool omegaMapNCD2::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61GroupCount*>(random_variable));
	return false;
}

bool omegaMapNCD2::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0: case 3:	//	mut
			// Insist on the particular Transformation (not just MatrixVariable) as it also specifies the size, encoding, etc.
			return(dynamic_cast<ParentDependentRateMatrix*>(parameter));
		case 1: case 4:	//	sel
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		case 2:	case 5: //	phylo
			return(dynamic_cast<NY98_TransProbMosaic*>(parameter));
		default:
			error("omegaMapNCD2::check_parameter_type(): parameter not found");
	}
	return false;
}

void omegaMapNCD2::set_mut(ParentDependentRateMatrix* mut, const int sp) {
	set_parameter(0+3*sp,(Variable*)mut);
}

void omegaMapNCD2::set_sel(ContinuousVectorVariable* sel, const int sp) {
	set_parameter(1+3*sp,(Variable*)sel);
}

void omegaMapNCD2::set_phylo(NY98_TransProbMosaic* phylo, const int sp) {
	set_parameter(2+3*sp,(Variable*)phylo);
}

const ParentDependentRateMatrix* omegaMapNCD2::get_mut(const int sp) const {
	return (const ParentDependentRateMatrix*)get_parameter(0+3*sp);
}

const ContinuousVectorVariable* omegaMapNCD2::get_sel(const int sp) const {
	return (const ContinuousVectorVariable*)get_parameter(1+3*sp);
}

const NY98_TransProbMosaic* omegaMapNCD2::get_phylo(const int sp) const {
	return (const NY98_TransProbMosaic*)get_parameter(2+3*sp);
}

mydouble omegaMapNCD2::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("omegaMapNCD2::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
	const Codon61GroupCount& ct = *((const Codon61GroupCount*)val);
	const ParentDependentRateMatrix* mut[2] = {get_mut(0), get_mut(1)};
	const ContinuousVectorVariable* sel[2] = {get_sel(0),get_sel(1)};
	const NY98_TransProbMosaic* phylo[2] = {get_phylo(0),get_phylo(1)};
	const double N[2] = {(double)(ct.n(0)),(double)(ct.n(1))};
	
	mydouble lik(1.0);
	int pos;
	const int L = ct.length();
//#pragma omp parallel for
	for(pos=0;pos<L;pos++) {
		// First, re-calculate species-specific conditional likelihoods, if necessary
		int sp;
		vector<bool> sp_changed(4,false);
		for(sp=0;sp<2;sp++) {
			/*if(true) {*/if(_mut_changed[sp] || (_sel_changed[sp] && sel[sp]->has_changed(pos))) {
				sp_changed[sp] = true;
				// Precalculate common factors in the likelihood at this site (i.e. irrespective of ancestral codon)
				double phi = lgamma(N[sp]+1.);
				int i;
				for(i=0;i<61;i++) {
					const int cti = ct[sp][pos][i];
					phi -= lgamma((double)(cti+1));
				}
				
				int anc;
				// Sum over the population ancestral amino acid
				for(anc=0;anc<61;anc++) {
					// Likelihood conditional on ancestor
					double likposanc = 0.0;
					// Ancestral amino acid
					const char aAA = omegaMapUtils::oneLetterCodes[anc];
					// Enumerate the partition
					double pt1 = 0;			// Number of observed individuals in that partition
					double ptheta = 0;		// Total mutation rate of codons in that partition
					double ttheta = 0;		// Total mutation rate
					for(i=0;i<61;i++) {
						const double cti = (double)ct[sp][pos][i];
						const double muti = mut[sp]->get_double(anc,i);
						ttheta += muti;
						if(omegaMapUtils::oneLetterCodes[i]==aAA) {
							pt1 += cti;
							ptheta += muti;
						}
						if(i==anc) likposanc += lgamma(cti+muti+1.0) - lgamma(muti+1.0);
						else likposanc += lgamma(cti+muti) - lgamma(muti);
					}
					likposanc += lgamma(ttheta) - lgamma(N[sp]+ttheta) + log(ptheta) - log(pt1+ptheta);					
					// Selection coefficient
					const double sigma = sel[sp]->get_double(pos);
					if(fabs(sigma)<1e-3) {
						// Ratio in the limit that sigma equals zero
						likposanc += log(pt1+ptheta) - log(N[sp]+ttheta) - log(ptheta) + log(ttheta);
					}
					else {
						const double hpnum = (omegaMapUtils::hypergeometric1F1(pt1+ptheta,N[sp]+ttheta,-sigma)-1.0);
						const double hpden = (omegaMapUtils::hypergeometric1F1(ptheta,ttheta,-sigma)-1.0);
						if(hpden==0.0 || hpnum/hpden<=0.0) {
							likposanc += log(pt1+ptheta) - log(N[sp]+ttheta) - log(ptheta) + log(ttheta);
						}
						else {
							likposanc += log(hpnum/hpden);
						}
					}
					// Store the within-population conditional likelihood (conditional on ancestral codon)
					_prune_lik[sp][ix][pos][anc].setlog(phi+likposanc);
				}	// anc
			}	// if needs updating
		}	// sp
		// Second, re-calculate the single-branch likelihoods conditional on the common ancestral state
		for(sp=0;sp<2;sp++) {
			if/*(true) {*/(sp_changed[sp] || (_phylo_changed[sp] && phylo[sp]->has_changed(pos))) {
				sp_changed[2+sp] = true;
				int anc;
				for(anc=0;anc<61;anc++) {
					_prune_lik[2+sp][ix][pos][anc] = 0.0;
					int dec;
					for(dec=0;dec<61;dec++) {
						_prune_lik[2+sp][ix][pos][anc] += _prune_lik[sp][ix][pos][dec] * mydouble(phylo[sp]->get_double(pos,anc,dec));
					}
				}
			}
		}
		// Third, re-calculate the common ancestral likelihood, summing over (i.e. marginal to) the ancestral species' state
		// Assume equal codon usage
		mydouble pifreq = 1./61.;
		if(sp_changed[2] || sp_changed[3]) {
			_likpos[ix][pos] = 0.0;
			int anc;
			for(anc=0;anc<61;anc++) {
				// Assume equal codon usage
				_likpos[ix][pos] += pifreq * _prune_lik[2][ix][pos][anc] * _prune_lik[3][ix][pos][anc];
			}
		}
		lik *= _likpos[ix][pos];
//		cout << "\t" << _likpos[ix][pos].LOG();
	}	//cout << endl;
	return lik;
}

void omegaMapNCD2::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61GroupCount* rv = dynamic_cast<Codon61GroupCount*>(random_variable);
	int sp;
	for(sp=0;sp<4;sp++) {
		_prune_lik[sp].push_back(vector< vector<mydouble> >(rv->length(),vector<mydouble>(61)));
		_previous_prune_lik[sp].push_back(vector< vector<mydouble> >(rv->length(),vector<mydouble>(61)));
	}
	_likpos.resize(n_random_variables(),rv->length());
	_previous_likpos.resize(n_random_variables(),rv->length());
}

void omegaMapNCD2::remove_random_variable(RandomVariable* random_variable) {
	error("omegaMapNCD2::remove_random_variable(): not permitted");
}

void omegaMapNCD2::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_mut_changed = vector<bool>(2,false);
		_sel_changed = vector<bool>(2,false);
		_phylo_changed = vector<bool>(2,false);
	}
	int sp;
	for(sp=0;sp<2;sp++) {
		if(v==(const Value*)get_mut(sp)) {
			if(sgl==Variable::_SET) {
				_mut_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_mut_changed[sp] = true;
				_previous_prune_lik[sp] = _prune_lik[sp];
				_previous_prune_lik[2+sp] = _prune_lik[2+sp];
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
//				_mut_changed[sp] = false;
			}
			else if(sgl==Variable::_REVERT) {
//				_mut_changed[sp] = false;
				_prune_lik[sp] = _previous_prune_lik[sp];
				_prune_lik[2+sp] = _previous_prune_lik[2+sp];
				_likpos = _previous_likpos;
			}
		}
		if(v==(const Value*)get_sel(sp)) {
			if(sgl==Variable::_SET) {
				_sel_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_sel_changed[sp] = true;
				_previous_prune_lik[sp] = _prune_lik[sp];
				_previous_prune_lik[2+sp] = _prune_lik[2+sp];
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
//				_sel_changed[sp] = false;
			}
			else if(sgl==Variable::_REVERT) {
//				_sel_changed[sp] = false;
				_prune_lik[sp] = _previous_prune_lik[sp];
				_prune_lik[2+sp] = _previous_prune_lik[2+sp];
				_likpos = _previous_likpos;
			}
		}
		if(v==(const Value*)get_phylo(sp)) {
			if(sgl==Variable::_SET) {
				_phylo_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_phylo_changed[sp] = true;
				_previous_prune_lik[2+sp] = _prune_lik[2+sp];
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
//				_phylo_changed[sp] = false;
			}
			else if(sgl==Variable::_REVERT) {
//				_phylo_changed[sp] = false;
				_prune_lik[2+sp] = _previous_prune_lik[2+sp];
				_likpos = _previous_likpos;
			}
		}
	}
	Distribution::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
