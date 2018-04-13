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
 *  mkprf3.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 09/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/mkprf3.h>
#include <omegaMap/RandomVariables/Codon61GroupCount.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string mkprf3ParameterNames[12] = {"thetaS0","thetaR0","gamma0","tau0",
										 "thetaS1","thetaR1","gamma1","tau1",
										 "thetaS2","thetaR2","gamma2","tau2"};

mkprf3::mkprf3(string name, DAG* dag) : DAGcomponent(name,dag,"mkprf3"), Distribution(mkprf3ParameterNames,12),
_thetaS_changed(vector<bool>(3,true)), _thetaR_changed(vector<bool>(3,true)), _gamma_changed(vector<bool>(3,true)),
_tau_changed(vector<bool>(3,true)) {
}

mkprf3::mkprf3(const mkprf3& x) : DAGcomponent(x), Distribution(x), 
_thetaS_changed(x._thetaS_changed), _thetaR_changed(x._thetaR_changed), _gamma_changed(x._gamma_changed),
_tau_changed(x._tau_changed), _likpos(x._likpos), _previous_likpos(x._previous_likpos), 
_site(x._site), _bad(x._bad) {
}

bool mkprf3::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61GroupCount*>(random_variable));
	return false;
}

bool mkprf3::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0: case 4: case 8:		//	thetaS
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1: case 5: case 9:		//	thetaR
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2: case 6: case 10:	//	gamma
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		case 3: case 7: case 11:	//	tau
			return(dynamic_cast<ContinuousVariable*>(parameter));
		default:
			error("mkprf3::check_parameter_type(): parameter not found");
	}
	return false;
}

void mkprf3::set_thetaS(ContinuousVariable* thetaS, const int sp) {
	set_parameter(0+4*sp,(Variable*)thetaS);
}

void mkprf3::set_thetaR(ContinuousVariable* thetaR, const int sp) {
	set_parameter(1+4*sp,(Variable*)thetaR);
}

void mkprf3::set_gamma(ContinuousVectorVariable* gamma, const int sp) {
	set_parameter(2+4*sp,(Variable*)gamma);
}

void mkprf3::set_tau(ContinuousVariable* tau, const int sp) {
	set_parameter(3+4*sp,(Variable*)tau);
}

const ContinuousVariable* mkprf3::get_thetaS(const int sp) const {
	return (const ContinuousVariable*)get_parameter(0+4*sp);
}

const ContinuousVariable* mkprf3::get_thetaR(const int sp) const {
	return (const ContinuousVariable*)get_parameter(1+4*sp);
}

const ContinuousVectorVariable* mkprf3::get_gamma(const int sp) const {
	return (const ContinuousVectorVariable*)get_parameter(2+4*sp);
}

const ContinuousVariable* mkprf3::get_tau(const int sp) const {
	return (const ContinuousVariable*)get_parameter(3+4*sp);
}

mydouble mkprf3::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("omegaMapNCD3::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
	const Codon61GroupCount& ct = *((const Codon61GroupCount*)val);
	const ContinuousVariable* thetaS[3] = {get_thetaS(0), get_thetaS(1), get_thetaS(2)};
	const ContinuousVariable* thetaR[3] = {get_thetaR(0), get_thetaR(1), get_thetaR(2)};
	const ContinuousVectorVariable* gamma[3] = {get_gamma(0), get_gamma(1), get_gamma(2)};
	const ContinuousVariable* tau[3] = {get_tau(0), get_tau(1), get_tau(2)};
	const double N[3] = {(double)(ct.n(0)),(double)(ct.n(1)),(double)(ct.n(2))};
	
	int pos,sp;
	const int L = ct.length();
	// Calculations!!!!!
	Matrix<double> lambda(3,5,0);
	for(sp=0;sp<3;sp++) {
		const double THETAS = thetaS[sp]->get_double();
		const double TAU = tau[sp]->get_double();
		lambda[sp][1] = F(N[sp],THETAS,0);
		lambda[sp][3] = THETAS*TAU/2.0+G(N[sp],THETAS,0);
	}
	mydouble lik(1);
	int preveval[3] = {-10,-10,-10};
	for(pos=0;pos<L;pos++) {
		for(sp=0;sp<3;sp++) {
			mydouble& likpossp = _likpos[ix][pos][sp];
			if(_thetaS_changed[sp] || _thetaR_changed[sp] || (_gamma_changed[sp] && gamma[sp]->has_changed(pos)) || _tau_changed[sp]) {
				const double GAMMA = gamma[sp]->get_double(pos);
				if(pos!=preveval[sp]+1 || GAMMA!=gamma[sp]->get_double(pos-1)) {
					const double OMEGA = (fabs(GAMMA)<1.0e-3) ? 1.0 : GAMMA/(1-exp(-GAMMA));
					const double THETAR = thetaR[sp]->get_double();
					const double TAU = tau[sp]->get_double();
					lambda[sp][2] = F(N[sp],THETAR,GAMMA);
					lambda[sp][4] = THETAR*TAU/2.0*OMEGA+G(N[sp],THETAR,GAMMA);
					preveval[sp] = pos;
				}
				if(!_bad[ix][pos]) {
					const int X = _site[ix][pos][sp];
					const double LAMBDA = lambda[sp][1]+lambda[sp][2]+lambda[sp][3]+lambda[sp][4];					
					if(X==0) likpossp.setlog(-LAMBDA);
					else likpossp = lambda[sp][X]/LAMBDA*(1-exp(-LAMBDA));
				}
				else likpossp = 1;
			}
			lik *= likpossp;
		}
	}
	return lik;
}

void mkprf3::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61GroupCount* rv = dynamic_cast<Codon61GroupCount*>(random_variable);
	const int L = rv->length();
	
	_likpos.push_back(vector< vector<mydouble> >(L,vector<mydouble>(3)));
	_previous_likpos.push_back(vector< vector<mydouble> >(L,vector<mydouble>(3)));
	_site.push_back(vector< vector<int> >(L,vector<int>(3)));
	_bad.push_back(vector<bool>(L));
	codify_sites(rv,ix);
	summarize(cout,rv,ix);
}

void mkprf3::remove_random_variable(RandomVariable* random_variable) {
	error("mkprf3::remove_random_variable(): not permitted");
}

void mkprf3::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_thetaS_changed = vector<bool>(3,false);
		_thetaR_changed = vector<bool>(3,false);
		_gamma_changed = vector<bool>(3,false);
		_tau_changed = vector<bool>(3,false);
	}
	int sp;
	for(sp=0;sp<3;sp++) {
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

void mkprf3::codify_sites(Codon61GroupCount* rv, const int ix) {
	const Codon61GroupCount& ct = *rv;
	const int L = ct.length();
	
	int cod, pos;
	// At each site, identify the ancestral and derived codon
	// under a strict definition in which the derived allele must be
	// unique to one species.
	for(pos=0;pos<L;pos++) {
		_site[ix][pos][0] = _site[ix][pos][1] = _site[ix][pos][2] = 0;
		// First determine that the site is no more than di-morphic
		for(cod=0;cod<61;cod++) if(ct[0][pos][cod]+ct[1][pos][cod]+ct[2][pos][cod]>0) break;
		const int al0 = cod;
		for(++cod;cod<61;cod++) if(ct[0][pos][cod]+ct[1][pos][cod]+ct[2][pos][cod]>0) break;
		const int al1 = cod;
		if(cod==61) {
			// Monomorphic (same as ancestor in each species)
			_bad[ix][pos] = false;
		}
		else {
			for(++cod;cod<61;cod++) if(ct[0][pos][cod]+ct[1][pos][cod]+ct[2][pos][cod]>0) break;
			if(cod!=61) {
				// More tha dimorphic
				_bad[ix][pos] = true;
			}
			else {
				// Synonymous or non-synonymous?
				bool syn = omegaMapUtils::oneLetterCodes[al0] == omegaMapUtils::oneLetterCodes[al1];
				// Polymorphic or fixed?
				bool poly;
				// Derived/fixed in which species?
				int sp;
				// Identify which allele is derived
				if((int)(ct[0][pos][al0]>0) + (int)(ct[1][pos][al0]>0) + (int)(ct[2][pos][al0]>0)>1) {
					// al0 is shared between species
					if((int)(ct[0][pos][al1]>0) + (int)(ct[1][pos][al1]>0) + (int)(ct[2][pos][al1]>0)>1) {
						// al1 is also shared between species
						_bad[ix][pos] = true;
					}
					else {
						// al1 is unique to one species, hence derived
						_bad[ix][pos] = false;
						sp = 1*(int)(ct[1][pos][al1]>0) + 2*(int)(ct[2][pos][al1]>0);
						poly = ct[sp][pos][al0]>0;
					}
				}
				else {
					// al0 is unique to one species, hence derived
					_bad[ix][pos] = false;
					sp = 1*(int)(ct[1][pos][al0]>0) + 2*(int)(ct[2][pos][al0]>0);
					poly = ct[sp][pos][al1]>0;
				}
				if(!_bad[ix][pos]) {
					if(poly & syn) _site[ix][pos][sp] = 1;
					else if(poly & !syn) _site[ix][pos][sp] = 2;
					else if(!poly & syn) _site[ix][pos][sp] = 3;
					else if(!poly & !syn) _site[ix][pos][sp] = 4;
				}
			}
		}
	}
}

double mkprf3::hypergeometric1F1(const double a, const double b, const double c) const {
	gsl_sf_result res;
	gsl_sf_hyperg_1F1_e(a, b, c, &res);
	return res.val;
}

double mkprf3::F(const double n, const double theta, const double gamma) const {
	double ret;
	if(fabs(gamma)<1.0e-3) ret = 1-(gsl_sf_beta(theta,n+theta+1)+gsl_sf_beta(n+theta,theta+1))/gsl_sf_beta(theta,theta+1);
	else ret = 1-(gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))+
					gsl_sf_beta(n+theta,theta)*(1-hypergeometric1F1(n+theta,n+2*theta,-gamma)))/
					gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	return ret;
}

double mkprf3::G(const double n, const double theta, const double gamma) const {
	double ret;
	if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(theta+1,n+theta)/gsl_sf_beta(theta+1,theta);
	else ret = gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))/
					gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	return ret;
}

void mkprf3::summarize(ostream& out, Codon61GroupCount* rv, const int ix) {
	int pos, nNA=0;
	vector<int> n1(3,0), n2(3,0), n3(3,0), n4(3,0);
	const int L = rv->length();
	out << "mkprf3 summary for rv " << rv->name() << endl;
	out << endl;
	out << "Uninterpretable sites:";
	for(pos=0;pos<L;pos++) {
		if(_bad[ix][pos]) {
			out << " " << pos;
			++nNA;
		}
	}	out << endl;
	out << "Synonymous polymorphisms (species):";
	for(pos=0;pos<L;pos++) {
		if(!_bad[ix][pos]) {
			if(_site[ix][pos][0]==1) {out << " " << pos << "(0)"; ++n1[0];}
			if(_site[ix][pos][1]==1) {out << " " << pos << "(1)"; ++n1[1];}
			if(_site[ix][pos][2]==1) {out << " " << pos << "(2)"; ++n1[2];}
		}
	}	out << endl;
	out << "Non-synonymous polymorphisms (species):";
	for(pos=0;pos<L;pos++) {
		if(!_bad[ix][pos]) {
			if(_site[ix][pos][0]==2) {out << " " << pos << "(0)"; ++n2[0];}
			if(_site[ix][pos][1]==2) {out << " " << pos << "(1)"; ++n2[1];}
			if(_site[ix][pos][2]==2) {out << " " << pos << "(2)"; ++n2[2];}
		}
	}	out << endl;
	out << "Synonymous substitutions (species):";
	for(pos=0;pos<L;pos++) {
		if(!_bad[ix][pos]) {
			if(_site[ix][pos][0]==3) {out << " " << pos << "(0)"; ++n3[0];}
			if(_site[ix][pos][1]==3) {out << " " << pos << "(1)"; ++n3[1];}
			if(_site[ix][pos][2]==3) {out << " " << pos << "(2)"; ++n3[2];}
		}
	}	out << endl;
	out << "Non-synonymous substitutions (species):";
	for(pos=0;pos<L;pos++) {
		if(!_bad[ix][pos]) {
			if(_site[ix][pos][0]==4) {out << " " << pos << "(0)"; ++n4[0];}
			if(_site[ix][pos][1]==4) {out << " " << pos << "(1)"; ++n4[1];}
			if(_site[ix][pos][2]==4) {out << " " << pos << "(2)"; ++n4[2];}
		}
	}	out << endl;
	char tab = '\t';
	out << endl;
	out << "Totals:" << endl;
	out << "0. N0 = " << L-nNA-n1[0]-n2[0]-n3[0]-n4[0] << tab << L-nNA-n1[1]-n2[1]-n3[1]-n4[1] << tab << L-nNA-n1[2]-n2[2]-n3[2]-n4[2] << endl;
	out << "1. PS = " << n1[0] << tab << n1[1] << tab << n1[2] << endl;
	out << "2. PN = " << n2[0] << tab << n2[1] << tab << n2[2] << endl;
	out << "3. DS = " << n3[0] << tab << n3[1] << tab << n3[2] << endl;
	out << "4. DN = " << n4[0] << tab << n4[1] << tab << n4[2] << endl;
	out << "Other = " << nNA << endl;
	out << endl;
}
	
} // namespace gcat_omegaMap

