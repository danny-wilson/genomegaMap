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
 *  omegaMap_unlinked.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 06/03/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/omegaMapUnlinked.h>
#include <omegaMap/RandomVariables/Codon61Count.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string omegaMapUnlinkedParameterNames[1] = {"mut"};

omegaMapUnlinked::omegaMapUnlinked(string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapUnlinked"), Distribution(omegaMapUnlinkedParameterNames,1), _mut_changed(true) {
}

omegaMapUnlinked::omegaMapUnlinked(const omegaMapUnlinked& x) : DAGcomponent(x), Distribution(x), _mut_changed(x._mut_changed), _likpos(x._likpos), _previous_likpos(x._previous_likpos), _rv_index(x._rv_index) {
}

bool omegaMapUnlinked::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61Count*>(random_variable));
	return false;
}

bool omegaMapUnlinked::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	//	mut
			// Insist on the particular Transformation (not just MatrixVariable) as it also specifies the size, encoding, etc.
			return(dynamic_cast<NY98_ParentDependentRateMatrix*>(parameter));
		default:
			error("omegaMapUnlinked::check_parameter_type(): parameter not found");
	}
	return false;
}

void omegaMapUnlinked::set_mut(NY98_ParentDependentRateMatrix* mut) {
	set_parameter(0,(Variable*)mut);
}

const NY98_ParentDependentRateMatrix* omegaMapUnlinked::get_mut() const {
	return (const NY98_ParentDependentRateMatrix*)get_parameter(0);
}

mydouble omegaMapUnlinked::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("omegaMapUnlinked::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
	const Codon61Count& ct = *((const Codon61Count*)val);
	const NY98_ParentDependentRateMatrix& mut = *get_mut();
	vector<double> pi = mut.get_pi()->get_doubles();
	
	mydouble lik(1.0);
	int pos;
    //cout << "pi =";
    //for(int i=0;i<61;i++) cout << " " << pi[i];
    //cout << endl;
    //cout << "theta = " << mut.get_theta()->get_double() << " ";
    //cout << "kappa = " << mut.get_kappa()->get_double() << " ";
    //cout << endl;
	for(pos=0;pos<ct.length();pos++) {
		if(_mut_changed && mut.has_changed(pos)) {
			// Precalculate common factors in the likelihood at this site (i.e. irrespective of ancestral codon)
			int i;
			double phi = 0., N = 0.;
            //if(pos==0) cout << "cti =";
			for(i=0;i<61;i++) {
				const int cti = ct[pos][i];
				phi -= lgamma((double)(cti+1));
				N += (double)cti;
                //if(pos==0) cout << " " << cti;
			}
			phi += lgamma(N+1.);
            //if(pos==0) cout << endl << "phi = " << phi << endl << "likposanc =";
			
			_likpos[ix][pos] = mydouble(0.0);
			int anc;
			// Sum over the population ancestral amino acid
			for(anc=0;anc<61;anc++) {
				// Likelihood conditional on ancestor, multiplied by ancestral probability
				double likposanc = log(pi[anc]);
				// Ancestral amino acid
				const char aAA = omegaMapUtils::oneLetterCodes[anc];
				// Enumerate the partition
				double pt1 = 0;			// Number of observed individuals in that partition
				double ptheta = 0;		// Total mutation rate of codons in that partition
				double ttheta = 0;		// Total mutation rate
				for(i=0;i<61;i++) {
					const double cti = (double)ct[pos][i];
					const double muti = mut.get_double(pos,anc,i);
					ttheta += muti;
					if(omegaMapUtils::oneLetterCodes[i]==aAA) {
						pt1 += cti;
						ptheta += muti;
					}
					if(i==anc) likposanc += lgamma(cti+muti+1.0) - lgamma(muti+1.0);
					else likposanc += lgamma(cti+muti) - lgamma(muti);
				}
				likposanc += lgamma(ttheta) - lgamma(N+ttheta) - log(N+ttheta) + log(ttheta);
                //if(pos==0) {
                    //cout << "anc " << anc << ": alpha =";
                    //for(i=0;i<61;i++) cout << " " << mut.get_double(pos,anc,i);
                    //cout << endl;
                    //cout << " " << likposanc;
                //}

				mydouble mylikposanc;
				mylikposanc.setlog(likposanc);
				_likpos[ix][pos] += mylikposanc;
			}
			// Multiply position-specific likelihood by common factors
			mydouble myphi;
			myphi.setlog(phi);
			_likpos[ix][pos] *= myphi;
            //if(pos==0) cout << endl << "myphi = " << myphi.LOG() << endl;
		}
		// Multiply likelihood by position-specific likelihood
		lik *= _likpos[ix][pos];
        //cout << "pos " << pos << ": omega = " << mut.get_omega()->get_double(pos) << " loglik = " << _likpos[ix][pos].LOG() << endl;
    }
	return lik;
}

void omegaMapUnlinked::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61Count* rv = dynamic_cast<Codon61Count*>(random_variable);
	_likpos.resize(n_random_variables(),rv->length());
	_previous_likpos.resize(n_random_variables(),rv->length());
}

void omegaMapUnlinked::remove_random_variable(RandomVariable* random_variable) {
	error("omegaMapUnlinked::remove_random_variable(): not permitted");
}

void omegaMapUnlinked::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
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
	Distribution::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
