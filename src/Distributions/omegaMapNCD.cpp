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
 *  omegaMapNCD.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 03/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/omegaMapNCD.h>
#include <omegaMap/RandomVariables/Codon61Count.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string omegaMapNCDParameterNames[2] = {"mut","sel"};

omegaMapNCD::omegaMapNCD(string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapNCD"), Distribution(omegaMapNCDParameterNames,2), _mut_changed(true), _sel_changed(true) {
}

omegaMapNCD::omegaMapNCD(const omegaMapNCD& x) : DAGcomponent(x), Distribution(x), _mut_changed(x._mut_changed), _sel_changed(x._sel_changed), _likpos(x._likpos), _previous_likpos(x._previous_likpos), _rv_index(x._rv_index) {
}

bool omegaMapNCD::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61Count*>(random_variable));
	return false;
}

bool omegaMapNCD::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	//	mut
			// Insist on the particular Transformation (not just MatrixVariable) as it also specifies the size, encoding, etc.
			return(dynamic_cast<ParentDependentRateMatrix*>(parameter));
		case 1:	//	sel
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("omegaMapNCD::check_parameter_type(): parameter not found");
	}
	return false;
}

void omegaMapNCD::set_mut(ParentDependentRateMatrix* mut) {
	set_parameter(0,(Variable*)mut);
}

void omegaMapNCD::set_sel(ContinuousVectorVariable* sel) {
	set_parameter(1,(Variable*)sel);
}

const ParentDependentRateMatrix* omegaMapNCD::get_mut() const {
	return (const ParentDependentRateMatrix*)get_parameter(0);
}

const ContinuousVectorVariable* omegaMapNCD::get_sel() const {
	return (const ContinuousVectorVariable*)get_parameter(1);
}

mydouble omegaMapNCD::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("omegaMapNCD::likelihood(): RandomVariable not recognised");
	const int ix = it->second;

	const Codon61Count& ct = *((const Codon61Count*)val);
	const ParentDependentRateMatrix& mut = *get_mut();
	const ContinuousVectorVariable& sel = *get_sel();
	const double N = (double)(ct.n());
	
	mydouble lik(1.0);
	int pos;
	for(pos=0;pos<ct.length();pos++) {
		if(_mut_changed || (_sel_changed && sel.has_changed(pos))) {
			// Precalculate common factors in the likelihood at this site (i.e. irrespective of ancestral codon)
			double phi = lgamma(N+1.);
			int i;
			for(i=0;i<61;i++) {
				const int cti = ct[pos][i];
				phi -= lgamma((double)(cti+1));
			}
			
			_likpos[ix][pos] = mydouble(0.0);
			int anc;
			// Sum over the population ancestral amino acid
			for(anc=0;anc<61;anc++) {
				// Likelihood conditional on ancestor (assume equal codon usage, for now!)
				double likposanc = -log(61.0);
				// Ancestral amino acid
				const char aAA = omegaMapUtils::oneLetterCodes[anc];
				// Enumerate the partition
				double pt1 = 0;			// Number of observed individuals in that partition
				double ptheta = 0;		// Total mutation rate of codons in that partition
				double ttheta = 0;		// Total mutation rate
				for(i=0;i<61;i++) {
					const double cti = (double)ct[pos][i];
					const double muti = mut.get_double(anc,i);
					ttheta += muti;
					if(omegaMapUtils::oneLetterCodes[i]==aAA) {
						pt1 += cti;
						ptheta += muti;
					}
					if(i==anc) likposanc += lgamma(cti+muti+1.0) - lgamma(muti+1.0);
					else likposanc += lgamma(cti+muti) - lgamma(muti);
				}
				likposanc += lgamma(ttheta) - lgamma(N+ttheta) + log(ptheta) - log(pt1+ptheta);
				// Selection coefficient
				const double sigma = sel.get_double(pos);
				if(fabs(sigma)<1e-3) {
					// Ratio in the limit that sigma equals zero
					likposanc += log(pt1+ptheta) - log(N+ttheta) - log(ptheta) + log(ttheta);
				}
				else {
					//30/3/15: from omegaMapNCD3.cpp
					const double hpnum = (omegaMapUtils::hypergeometric1F1(pt1+ptheta,N+ttheta,-sigma)-1.0);
					const double hpden = (omegaMapUtils::hypergeometric1F1(ptheta,ttheta,-sigma)-1.0);
//					const double hpnum = (gsl_sf_hyperg_1F1(pt1+ptheta,N+ttheta,-sigma)-1.0);
//					const double hpden = (gsl_sf_hyperg_1F1(ptheta,ttheta,-sigma)-1.0);
					if(hpden==0.0 || hpnum/hpden<=0.0) {
						likposanc += log(pt1+ptheta) - log(N+ttheta) - log(ptheta) + log(ttheta);
					}
					else {
						likposanc += log(hpnum/hpden);
					}
				}
				mydouble mylikposanc;
				mylikposanc.setlog(likposanc);
				_likpos[ix][pos] += mylikposanc;
			}
			// Multiply position-specific likelihood by common factors
			mydouble myphi;
			myphi.setlog(phi);
			_likpos[ix][pos] *= myphi;
		}
		// Multiply likelihood by position-specific likelihood
		lik *= _likpos[ix][pos];
	}
	return lik;
}

void omegaMapNCD::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61Count* rv = dynamic_cast<Codon61Count*>(random_variable);
	_likpos.resize(n_random_variables(),rv->length());
	_previous_likpos.resize(n_random_variables(),rv->length());
}

void omegaMapNCD::remove_random_variable(RandomVariable* random_variable) {
	error("omegaMapNCD::remove_random_variable(): not permitted");
}

void omegaMapNCD::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
//	if(v==(const Value*)get_mut()) _mut_changed = true;
//	if(v==(const Value*)get_sel()) _sel_changed = true;
	if(v==(const Value*)get_mut()) {
		if(sgl==Variable::_SET) {
			_mut_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_mut_changed = true;
			_previous_likpos = _likpos;
		}
		else if(sgl==Variable::_ACCEPT) {
			_mut_changed = false;
		}
		else if(sgl==Variable::_REVERT) {
			_mut_changed = false;
			_likpos = _previous_likpos;
		}
	}
	if(v==(const Value*)get_sel()) {
		if(sgl==Variable::_SET) {
			_sel_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_sel_changed = true;
			_previous_likpos = _likpos;
		}
		else if(sgl==Variable::_ACCEPT) {
			_sel_changed = false;
		}
		else if(sgl==Variable::_REVERT) {
			_sel_changed = false;
			_likpos = _previous_likpos;
		}
	}
	Distribution::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
