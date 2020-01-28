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
 *  BiallelicNCD2.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 23/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/BiallelicNCD2.h>
#include <omegaMap/RandomVariables/BiallelicCodingGroupCount.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const int NSP = 2;

const string BiallelicNCD2ParameterNames[4*NSP] = {	"thetaS0","thetaR0","gamma0","tau0",
													"thetaS1","thetaR1","gamma1","tau1"};

BiallelicNCD2::BiallelicNCD2(string name, DAG* dag) : DAGcomponent(name,dag,"BiallelicNCD2"), Distribution(BiallelicNCD2ParameterNames,4*NSP),
_thetaS_changed(vector<bool>(NSP,true)), _thetaR_changed(vector<bool>(NSP,true)), _gamma_changed(vector<bool>(NSP,true)),
_tau_changed(vector<bool>(NSP,true)) {
}

BiallelicNCD2::BiallelicNCD2(const BiallelicNCD2& x) : DAGcomponent(x), Distribution(x), 
_thetaS_changed(x._thetaS_changed), _thetaR_changed(x._thetaR_changed), _gamma_changed(x._gamma_changed),
_tau_changed(x._tau_changed), _likpos(x._likpos), _previous_likpos(x._previous_likpos) {
}

bool BiallelicNCD2::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<BiallelicCodingGroupCount*>(random_variable));
	return false;
}

bool BiallelicNCD2::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0: case 4:		//	thetaS
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1: case 5:		//	thetaR
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2: case 6:		//	gamma
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		case 3: case 7:		//	tau
			return(dynamic_cast<ContinuousVariable*>(parameter));
		default:
			error("BiallelicNCD2::check_parameter_type(): parameter not found");
	}
	return false;
}

void BiallelicNCD2::set_thetaS(ContinuousVariable* thetaS, const int sp) {
	set_parameter(0+4*sp,(Variable*)thetaS);
}

void BiallelicNCD2::set_thetaR(ContinuousVariable* thetaR, const int sp) {
	set_parameter(1+4*sp,(Variable*)thetaR);
}

void BiallelicNCD2::set_gamma(ContinuousVectorVariable* gamma, const int sp) {
	set_parameter(2+4*sp,(Variable*)gamma);
}

void BiallelicNCD2::set_tau(ContinuousVariable* tau, const int sp) {
	set_parameter(3+4*sp,(Variable*)tau);
}

const ContinuousVariable* BiallelicNCD2::get_thetaS(const int sp) const {
	return (const ContinuousVariable*)get_parameter(0+4*sp);
}

const ContinuousVariable* BiallelicNCD2::get_thetaR(const int sp) const {
	return (const ContinuousVariable*)get_parameter(1+4*sp);
}

const ContinuousVectorVariable* BiallelicNCD2::get_gamma(const int sp) const {
	return (const ContinuousVectorVariable*)get_parameter(2+4*sp);
}

const ContinuousVariable* BiallelicNCD2::get_tau(const int sp) const {
	return (const ContinuousVariable*)get_parameter(3+4*sp);
}

mydouble BiallelicNCD2::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("BiallelicNCD2::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
	const BiallelicCodingGroupCount& ct = *((const BiallelicCodingGroupCount*)val);
	// Conversion factor to equate MK thetas to site-specific thetas
	// The correction factor includes L/3 = ct.length() to convert from per-gene to per-codon, and
	// a factor of 1/2 which reflects that each codon can only be (with equal probability)
	// synonymous OR non-synonymous (PRF sort of assumes each site can be either, because
	// events are so rare).
	const double k = (double)(ct.length())/2.0;
	const double thetaS[NSP] = {get_thetaS(0)->get_double()/k, get_thetaS(1)->get_double()/k};
	const double thetaR[NSP] = {get_thetaR(0)->get_double()/k, get_thetaR(1)->get_double()/k};
	const double tau[NSP] = {get_tau(0)->get_double(), get_tau(1)->get_double()};
	const double N[NSP] = {(double)(ct.n(0)),(double)(ct.n(1))};
	const double gamma_neut[NSP] = {0.0, 0.0};
	
	bool thetaS_changed = _thetaS_changed[0] || _thetaS_changed[1];
	bool thetaR_changed = _thetaR_changed[0] || _thetaR_changed[1];
	bool tau_changed = _tau_changed[0] || _tau_changed[1];
	bool gamma_changed = _gamma_changed[0] || _gamma_changed[1];
	mydouble psyn  = 0.5; //mydouble((thetaS[0]+thetaS[1])/(thetaS[0]+thetaS[1]+thetaR[0]+thetaR[1]));
	mydouble pnsyn = 0.5; //mydouble((thetaR[0]+thetaR[1])/(thetaS[0]+thetaS[1]+thetaR[0]+thetaR[1]));
	int pos;
	const int L = ct.length();
	mydouble ret(1);
	for(pos=0;pos<L;pos++) {
		mydouble& likpos = _likpos[ix][pos][0];
		// Is the site valid?
		bool valid = ct.is_valid(pos);
		// Is the site variable?
		bool vari = ct.codon0(pos)!=ct.codon1(pos) && valid;
		// Is the site non-synonymous?
		bool nsyn = ct.aminoacid0(pos)!=ct.aminoacid1(pos) && vari;
		// Is the site synonymous?
		bool syn = ct.aminoacid0(pos)==ct.aminoacid1(pos) && vari;
		// Do I need to recalculate at this position?
//		if(!vari || tau_changed  || thetaS_changed || thetaR_changed || (gamma_changed && (get_gamma(0)->has_changed(pos) || get_gamma(1)->has_changed(pos)))) {
		if(!vari || tau_changed  || (syn && thetaS_changed) || (nsyn && thetaR_changed) || (nsyn && gamma_changed && (get_gamma(0)->has_changed(pos) || get_gamma(1)->has_changed(pos)))) {
//		if(true) {
			if(!valid) {
				likpos = mydouble(1);
			}
			else {
				mydouble likpos_nsyn = 0;
				mydouble likpos_syn = 0;
				const double X[2] = {ct[0][pos][0], ct[1][pos][0]};
				if(nsyn || !vari) {
					const double gamma[NSP] = {get_gamma(0)->get_double(pos), get_gamma(1)->get_double(pos)};
					likpos_nsyn = site_likelihood(N,X,thetaR,tau,gamma);
				}
				if(syn || !vari) {
					likpos_syn = site_likelihood(N,X,thetaS,tau,gamma_neut);
				}
				likpos = psyn*likpos_syn+pnsyn*likpos_nsyn;
			}
		}
		ret *= likpos;
	}
	
	return ret;
}

mydouble BiallelicNCD2::site_likelihood(const double N[], const double x[], const double theta[], const double tau[], const double gamma[]) const {
	// Population-specific likelihood indexed by population then ancestral allele
	Matrix<mydouble> poplik(2,2);
	poplik[0][0] = acsf(N[0],x[0],0,theta[0],gamma[0]);
	poplik[0][1] = acsf(N[0],x[0],1,theta[0],gamma[0]);
	poplik[1][0] = acsf(N[1],x[1],0,theta[1],gamma[1]);
	poplik[1][1] = acsf(N[1],x[1],1,theta[1],gamma[1]);
	// DN/DS ratio
	double omega[2];
	omega[0] = (fabs(gamma[0])<1.0e-3) ? 1.0 : gamma[0]/(1-exp(-gamma[0]));
	omega[1] = (fabs(gamma[1])<1.0e-3) ? 1.0 : gamma[1]/(1-exp(-gamma[1]));
	// Phylogenetic substitution probabilities, indexed by lineage
	mydouble psame[2] = {mydouble(0.5*(1+exp(-theta[0]*omega[0]*tau[0]))), mydouble(0.5*(1+exp(-theta[1]*omega[1]*tau[1])))};
	mydouble pdiff[2] = {mydouble(0.5*(1-exp(-theta[0]*omega[0]*tau[0]))), mydouble(0.5*(1-exp(-theta[1]*omega[1]*tau[1])))};
	// Marginal likelihood conditional on MRCA allele 0
	mydouble likanc0 = (psame[0]*poplik[0][0]+pdiff[0]*poplik[0][1])*(psame[1]*poplik[1][0]+pdiff[1]*poplik[1][1]);
	// Marginal likelihood conditional on MRCA allele 1
	mydouble likanc1 = (psame[0]*poplik[0][1]+pdiff[0]*poplik[0][0])*(psame[1]*poplik[1][1]+pdiff[1]*poplik[1][0]);
	// Site likelihood marginal to ancestor
	return mydouble(0.5)*(likanc0+likanc1);
}

void BiallelicNCD2::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	BiallelicCodingGroupCount* rv = dynamic_cast<BiallelicCodingGroupCount*>(random_variable);
	const int L = rv->length();
	
	_likpos.push_back(vector< vector<mydouble> >(L,vector<mydouble>(NSP)));
	_previous_likpos.push_back(vector< vector<mydouble> >(L,vector<mydouble>(NSP)));
}

void BiallelicNCD2::remove_random_variable(RandomVariable* random_variable) {
	error("BiallelicNCD2::remove_random_variable(): not permitted");
}

void BiallelicNCD2::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_thetaS_changed = vector<bool>(NSP,false);
		_thetaR_changed = vector<bool>(NSP,false);
		_gamma_changed = vector<bool>(NSP,false);
		_tau_changed = vector<bool>(NSP,false);
	}
	int sp;
	for(sp=0;sp<NSP;sp++) {
		if(v==(const Value*)get_thetaS(sp)) {
			if(sgl==Variable::_SET) {
				_thetaS_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_thetaS_changed[sp] = true;
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
			}
		}
		if(v==(const Value*)get_thetaR(sp)) {
			if(sgl==Variable::_SET) {
				_thetaR_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_thetaR_changed[sp] = true;
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
			}
		}
		if(v==(const Value*)get_gamma(sp)) {
			if(sgl==Variable::_SET) {
				_gamma_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_gamma_changed[sp] = true;
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
			}
		}
		if(v==(const Value*)get_tau(sp)) {
			if(sgl==Variable::_SET) {
				_tau_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_tau_changed[sp] = true;
				_previous_likpos = _likpos;
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
			}
		}
	}
	Distribution::receive_signal_from_parent(v,sgl);
}

double BiallelicNCD2::hypergeometric1F1(const double a, const double b, const double c) const {
	gsl_sf_result res;
	gsl_sf_hyperg_1F1_e(a, b, c, &res);
	return res.val;
}

// Ancestrally conditioned sampling formula
//	n		Number of sequences
//	x		Number of copies of allele 0
//	anc		Ancestral allele (0 or 1)
//	theta	Population-scaled mutation rate
//	gamma	Population-scaled selection coefficient
double BiallelicNCD2::acsf(const double n, const double x, const int anc, const double theta, const double gamma) const {
	double ret;
	if(anc==0) {
		if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(x+theta+1,n-x+theta)/gsl_sf_beta(theta+1,theta);
		else ret = gsl_sf_beta(x+theta,n-x+theta)*(1-hypergeometric1F1(x+theta,n+2*theta,-gamma))/
			gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	}
	else {
		if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(n-x+theta+1,x+theta)/gsl_sf_beta(theta+1,theta);
		else ret = gsl_sf_beta(n-x+theta,x+theta)*(1-hypergeometric1F1(n-x+theta,n+2*theta,-gamma))/
			gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	}
	return ret;
}

// Ancestrally conditioned sampling formula, pooling polymorphic categories into one
//	n		Number of sequences
//	x		Number of copies of allele 0
//	anc		Ancestral allele (0 or 1)
//	theta	Population-scaled mutation rate
//	gamma	Population-scaled selection coefficient
double BiallelicNCD2::acsf_pool(const double n, const double x, const int anc, const double theta, const double gamma) const {
	const bool FIXED_FOR_ANC = true, NOT_FIXED_FOR_ANC = false;
	if(anc==0) {
		if(x==0) return pfxd(n,NOT_FIXED_FOR_ANC,theta,gamma);
		if(x==n) return pfxd(n,FIXED_FOR_ANC,theta,gamma);
		return pseg(n,theta,gamma);
	}
	else {
		if(x==0) return pfxd(n,FIXED_FOR_ANC,theta,gamma);
		if(x==n) return pfxd(n,NOT_FIXED_FOR_ANC,theta,gamma);
		return pseg(n,theta,gamma);
	}
}

// Probability that a site is segregating in a population
double BiallelicNCD2::pseg(const double n, const double theta, const double gamma) const {
	double ret;
	if(fabs(gamma)<1.0e-3) ret = 1-(gsl_sf_beta(theta,n+theta+1)+gsl_sf_beta(n+theta,theta+1))/gsl_sf_beta(theta,theta+1);
	else ret = 1-(gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))+
				  gsl_sf_beta(n+theta,theta)*(1-hypergeometric1F1(n+theta,n+2*theta,-gamma)))/
		gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	return ret;
}

// Probability that a site is fixed in a population
double BiallelicNCD2::pfxd(const double n, const bool fixed_for_anc, const double theta, const double gamma) const {
	double ret;
	if(fixed_for_anc) {
		if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(theta,n+theta+1)/gsl_sf_beta(theta,theta+1);
		else ret =	gsl_sf_beta(n+theta,theta)*(1-hypergeometric1F1(n+theta,n+2*theta,-gamma))/
					gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	}
	else {
		if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(n+theta,theta+1)/gsl_sf_beta(theta,theta+1);
		else ret =	gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))/
					gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	}
	return ret;
}

} // namespace gcat_omegaMap
