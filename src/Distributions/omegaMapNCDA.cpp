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
 *  omegaMapNCDA.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/omegaMapNCDA.h>
#include <omegaMap/RandomVariables/Codon61AncestralCount.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string omegaMapNCDAParameterNames[2] = {"mut","sel"};

omegaMapNCDA::omegaMapNCDA(string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapNCDA"), Distribution(omegaMapNCDAParameterNames,2), _mut_changed(true), _sel_changed(true) {
}

omegaMapNCDA::omegaMapNCDA(const omegaMapNCDA& x) : DAGcomponent(x), Distribution(x), _mut_changed(x._mut_changed), _sel_changed(x._sel_changed) {
}

bool omegaMapNCDA::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61AncestralCount*>(random_variable));
	return false;
}

bool omegaMapNCDA::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	//	mut
			// Insist on the particular Transformation (not just MatrixVariable) as it also specifies the size, encoding, etc.
			return(dynamic_cast<ParentDependentRateMatrix*>(parameter));
		case 1:	//	sel
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		default:
			error("omegaMapNCDA::check_parameter_type(): parameter not found");
	}
	return false;
}

void omegaMapNCDA::set_mut(ParentDependentRateMatrix* mut) {
	set_parameter(0,(Variable*)mut);
}

void omegaMapNCDA::set_sel(ContinuousVectorVariable* sel) {
	set_parameter(1,(Variable*)sel);
}

const ParentDependentRateMatrix* omegaMapNCDA::get_mut() const {
	return (const ParentDependentRateMatrix*)get_parameter(0);
}

const ContinuousVectorVariable* omegaMapNCDA::get_sel() const {
	return (const ContinuousVectorVariable*)get_parameter(1);
}

mydouble omegaMapNCDA::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	const Codon61AncestralCount& ct = *((const Codon61AncestralCount*)val);
	const ParentDependentRateMatrix& mut = *get_mut();
	const ContinuousVectorVariable& sel = *get_sel();
	const double N = (double)(ct.n());
	
	double lik = 0.0;
	int pos;
	for(pos=0;pos<ct.length();pos++) {
		// Ancestral codon and its amino acid
		const int anc = ct.ancestor(pos);
		const char aAA = omegaMapUtils::oneLetterCodes[anc];
		
		// Enumerate the partition
		double pt1 = 0;			// Number of observed individuals in that partition
		double ptheta = 0;		// Total mutation rate of codons in that partition
		double ttheta = 0;		// Total mutation rate

		double phi = lgamma(N+1.);
		double likpos = 0.0;
		int i;
		for(i=0;i<61;i++) {
			const double cti = (double)ct[pos][i];
			const double muti = mut.get_double(anc,i);
			ttheta += muti;
			if(omegaMapUtils::oneLetterCodes[i]==aAA) {
				pt1 += cti;
				ptheta += muti;
			}
			phi -= lgamma(cti+1.0);
			if(i==anc) likpos += lgamma(cti+muti+1.0) - lgamma(muti+1.0);
			else likpos += lgamma(cti+muti) - lgamma(muti);
		}
		likpos += lgamma(ttheta) - lgamma(N+ttheta) + log(ptheta) - log(pt1+ptheta);

		// Selection coefficient
		const double sigma = sel.get_double(pos);
		if(fabs(sigma)<1e-3) {
			// Ratio in the limit that sigma equals zero
			likpos += log(pt1+ptheta) - log(N+ttheta) - log(ptheta) + log(ttheta);
		}
		else {
			const double hpnum = (gsl_sf_hyperg_1F1(pt1+ptheta,N+ttheta,-sigma)-1.0);
			const double hpden = (gsl_sf_hyperg_1F1(ptheta,ttheta,-sigma)-1.0);
			if(hpden==0.0 || hpnum/hpden<=0.0) {
				likpos += log(pt1+ptheta) - log(N+ttheta) - log(ptheta) + log(ttheta);
			}
			else {
				likpos += log(hpnum/hpden);
			}
		}
		lik += phi+likpos;
	}
	mydouble ret;
	ret.setlog(lik);
	return ret;
}

void omegaMapNCDA::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(v==(const Value*)get_mut()) _mut_changed = true;
	if(v==(const Value*)get_sel()) _sel_changed = true;
	Distribution::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
