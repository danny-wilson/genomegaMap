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
 *  mkprf3HMMmagmag.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 24/02/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/mkprf3HMMmag.h>
#include <omegaMap/RandomVariables/Codon61GroupCount.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string mkprf3HMMmagParameterNames[18] = {"thetaS0","thetaR0","gamma0","tau0","p0","q0",
												"thetaS1","thetaR1","gamma1","tau1","p1","q1",
												"thetaS2","thetaR2","gamma2","tau2","p2","q2"};

mkprf3HMMmag::mkprf3HMMmag(string name, DAG* dag) : DAGcomponent(name,dag,"mkprf3HMMmag"), Distribution(mkprf3HMMmagParameterNames,18),
_thetaS_changed(vector<bool>(3,true)), _thetaR_changed(vector<bool>(3,true)), _gamma_changed(vector<bool>(3,true)),
_tau_changed(vector<bool>(3,true)), _p_changed(vector<bool>(3,true)), _q_changed(vector<bool>(3,true)) {
}

mkprf3HMMmag::mkprf3HMMmag(const mkprf3HMMmag& x) : DAGcomponent(x), Distribution(x), 
_thetaS_changed(x._thetaS_changed), _thetaR_changed(x._thetaR_changed), _gamma_changed(x._gamma_changed),
_tau_changed(x._tau_changed), _p_changed(x._p_changed), _q_changed(x._q_changed), _likpos(x._likpos), _previous_likpos(x._previous_likpos), 
_site(x._site), _bad(x._bad), _alpha(x._alpha), _previous_alpha(x._previous_alpha), _pneg(x._pneg), _previous_pneg(x._previous_pneg) {
}

bool mkprf3HMMmag::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61GroupCount*>(random_variable));
	return false;
}

bool mkprf3HMMmag::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0: case 6: case 12:	//	thetaS
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 1: case 7: case 13:	//	thetaR
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2: case 8: case 14:	//	gamma
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		case 3: case 9: case 15:	//	tau
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 4: case 10: case 16:	//	p
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 5: case 11: case 17:	//	q
			return(dynamic_cast<ContinuousVariable*>(parameter));
		default:
			error("mkprf3HMMmag::check_parameter_type(): parameter not found");
	}
	return false;
}

void mkprf3HMMmag::set_thetaS(ContinuousVariable* thetaS, const int sp) {
	set_parameter(0+6*sp,(Variable*)thetaS);
}

void mkprf3HMMmag::set_thetaR(ContinuousVariable* thetaR, const int sp) {
	set_parameter(1+6*sp,(Variable*)thetaR);
}

void mkprf3HMMmag::set_gamma(ContinuousVectorVariable* gamma, const int sp) {
	set_parameter(2+6*sp,(Variable*)gamma);
}

void mkprf3HMMmag::set_tau(ContinuousVariable* tau, const int sp) {
	set_parameter(3+6*sp,(Variable*)tau);
}

void mkprf3HMMmag::set_p(ContinuousVariable* p, const int sp) {
	set_parameter(4+6*sp,(Variable*)p);
}

void mkprf3HMMmag::set_q(ContinuousVariable* q, const int sp) {
	set_parameter(5+6*sp,(Variable*)q);
}

const ContinuousVariable* mkprf3HMMmag::get_thetaS(const int sp) const {
	return (const ContinuousVariable*)get_parameter(0+6*sp);
}

const ContinuousVariable* mkprf3HMMmag::get_thetaR(const int sp) const {
	return (const ContinuousVariable*)get_parameter(1+6*sp);
}

const ContinuousVectorVariable* mkprf3HMMmag::get_gamma(const int sp) const {
	return (const ContinuousVectorVariable*)get_parameter(2+6*sp);
}

const ContinuousVariable* mkprf3HMMmag::get_tau(const int sp) const {
	return (const ContinuousVariable*)get_parameter(3+6*sp);
}

const ContinuousVariable* mkprf3HMMmag::get_p(const int sp) const {
	return (const ContinuousVariable*)get_parameter(4+6*sp);
}

const ContinuousVariable* mkprf3HMMmag::get_q(const int sp) const {
	return (const ContinuousVariable*)get_parameter(5+6*sp);
}

mydouble mkprf3HMMmag::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("mkprf3HMMmag::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
	const Codon61GroupCount& ct = *((const Codon61GroupCount*)val);
	const int L = ct.length();
	
	int pos,sp,g;
	Vector<double> lambda(5);
	Matrix<mydouble> pemiss;
	mydouble lik(1);
	for(sp=0;sp<3;sp++) {
		mydouble& liksp = _likpos[ix][0][sp];
		if(_thetaS_changed[sp] || _thetaR_changed[sp] || _gamma_changed[sp] || _tau_changed[sp] || _p_changed[sp] || _q_changed[sp]) {
			// Pre-calculations!!!!!
			const double N = (double)(ct.n(sp));
			const int ng = get_gamma(sp)->length();
			Matrix<mydouble>& alpha = _alpha[ix][sp];
			alpha = Matrix<mydouble>(L,ng);
			Matrix<mydouble>& pneg = _pneg[ix][sp];
			pneg = Matrix<mydouble>(L,ng);
			const double THETAS = get_thetaS(sp)->get_double();
			const double THETAR = get_thetaR(sp)->get_double();
			const double TAU = get_tau(sp)->get_double();
			mydouble Q = get_q(sp)->get_double();
			mydouble OneMinusQ = 1.0-get_q(sp)->get_double();
			lambda[1] = MAX(0.0,F(N,THETAS,0));
			lambda[3] = MAX(0.0,THETAS*TAU/2.0+G(N,THETAS,0));
			pemiss = Matrix<mydouble>(5,2*ng);
			for(g=0;g<ng;g++) {
				// When gamma > 0
				double GAMMA = get_gamma(sp)->get_double(g);
				double OMEGA = (fabs(GAMMA)<1.0e-3) ? 1.0 : GAMMA/(1-exp(-GAMMA));
				lambda[2] = MAX(0.0,F(N,THETAR,GAMMA));
				lambda[4] = MAX(0.0,THETAR*TAU/2.0*OMEGA+G(N,THETAR,GAMMA));
				lambda[0] = lambda[1]+lambda[2]+lambda[3]+lambda[4];
				int x;
				pemiss[0][g].setlog(-lambda[0]);
				for(x=1;x<5;x++) {
					pemiss[x][g] = lambda[x]/lambda[0]*(1-exp(-lambda[0]));
				}
				// When gamma < 0
				GAMMA = -GAMMA;
				OMEGA = (fabs(GAMMA)<1.0e-3) ? 1.0 : GAMMA/(1-exp(-GAMMA));
				lambda[2] = MAX(0.0,F(N,THETAR,GAMMA));
				lambda[4] = MAX(0.0,THETAR*TAU/2.0*OMEGA+G(N,THETAR,GAMMA));
				lambda[0] = lambda[1]+lambda[2]+lambda[3]+lambda[4];
				pemiss[0][ng+g].setlog(-lambda[0]);
				for(x=1;x<5;x++) {
					pemiss[x][ng+g] = lambda[x]/lambda[0]*(1-exp(-lambda[0]));
				}
				// Now let pemiss[x][g] be the unconditional probability
				// and pemiss[x][ng+g] the conditional probaility of negative selection
				for(x=0;x<5;x++) {
					pemiss[x][g] = Q*pemiss[x][g]+OneMinusQ*pemiss[x][ng+g];
					pemiss[x][ng+g] = (pemiss[x][g].iszero()) ? OneMinusQ : OneMinusQ*pemiss[x][ng+g]/pemiss[x][g];
				}
			}
			mydouble pswitch = get_p(sp)->get_double();
			mydouble pnoswitch = 1-pswitch;
			mydouble pi = 1.0/(double)(ng);
			mydouble sumalpha;
			// Forward algorithm
			{	// pos=0
				if(!_bad[ix][0]) {
					const int X = _site[ix][0][sp];
					sumalpha = 0;
					for(g=0;g<ng;g++) {
						alpha[0][g] = pemiss[X][g]*pi;
						sumalpha += alpha[0][g];
						pneg[0][g] = pemiss[X][ng+g];
					}
				}
				else {
					sumalpha = 0;
					for(g=0;g<ng;g++) {
						alpha[0][g] = pi;
						sumalpha += alpha[0][g];
						pneg[0][g] = OneMinusQ;
					}
				}
			}
			for(pos=1;pos<L;pos++) {
				if(!_bad[ix][pos]) {
					const int X = _site[ix][pos][sp];
					mydouble newsumalpha = 0;
					for(g=0;g<ng;g++) {
						alpha[pos][g] = (pswitch*sumalpha*pi+pnoswitch*alpha[pos-1][g])*pemiss[X][g];
						newsumalpha += alpha[pos][g];
						pneg[pos][g] = pemiss[X][ng+g];
					}
					sumalpha = newsumalpha;
				}
				else {
					mydouble newsumalpha = 0;
					for(g=0;g<ng;g++) {
						alpha[pos][g] = (pswitch*sumalpha*pi+pnoswitch*alpha[pos-1][g]);
						newsumalpha += alpha[pos][g];
						pneg[pos][g] = OneMinusQ;
					}
					sumalpha = newsumalpha;
				}
			}
			liksp = sumalpha;
		}
		lik *= liksp;
	}
	return lik;
}

void mkprf3HMMmag::add_random_variable(RandomVariable* random_variable) {
	if(n_random_variables()>0) {
		string errTxt = "mkprf3HMMmag::add_random_variable(): only one child RV allowed for object " + name();
		error(errTxt.c_str());
	}
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61GroupCount* rv = dynamic_cast<Codon61GroupCount*>(random_variable);
	const int L = rv->length();
	
	_likpos.push_back(vector< vector<mydouble> >(1,vector<mydouble>(3)));
	_previous_likpos.push_back(vector< vector<mydouble> >(1,vector<mydouble>(3)));
	_site.push_back(vector< vector<int> >(L,vector<int>(3)));
	_bad.push_back(vector<bool>(L));
	codify_sites(rv,ix);
	summarize(cout,rv,ix);
	
	_alpha.resize(n_random_variables(), 3);
	_previous_alpha.resize(n_random_variables(), 3);
	_pneg.resize(n_random_variables(), 3);
	_previous_pneg.resize(n_random_variables(), 3);
}

void mkprf3HMMmag::remove_random_variable(RandomVariable* random_variable) {
	error("mkprf3HMMmag::remove_random_variable(): not permitted");
}

void mkprf3HMMmag::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_thetaS_changed = vector<bool>(3,false);
		_thetaR_changed = vector<bool>(3,false);
		_gamma_changed = vector<bool>(3,false);
		_tau_changed = vector<bool>(3,false);
		_p_changed = vector<bool>(3,false);
		_q_changed = vector<bool>(3,false);
	}
	int sp,ix;
	for(sp=0;sp<3;sp++) {
		if(v==(const Value*)get_thetaS(sp)) {
			if(sgl==Variable::_SET) {
				_thetaS_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_thetaS_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
		if(v==(const Value*)get_thetaR(sp)) {
			if(sgl==Variable::_SET) {
				_thetaR_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_thetaR_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
		if(v==(const Value*)get_gamma(sp)) {
			if(sgl==Variable::_SET) {
				_gamma_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_gamma_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
		if(v==(const Value*)get_tau(sp)) {
			if(sgl==Variable::_SET) {
				_tau_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_tau_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
		if(v==(const Value*)get_p(sp)) {
			if(sgl==Variable::_SET) {
				_p_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_p_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
		if(v==(const Value*)get_q(sp)) {
			if(sgl==Variable::_SET) {
				_q_changed[sp] = true;
			}
			else if(sgl==Variable::_PROPOSE) {
				_q_changed[sp] = true;
				_previous_likpos = _likpos;
				for(ix=0;ix<n_random_variables();ix++) _previous_alpha[ix][sp] = _alpha[ix][sp];
			}
			else if(sgl==Variable::_ACCEPT) {
			}
			else if(sgl==Variable::_REVERT) {
				_likpos = _previous_likpos;
				for(ix=0;ix<n_random_variables();ix++) _alpha[ix][sp] = _previous_alpha[ix][sp];
			}
		}
	}
	Distribution::receive_signal_from_parent(v,sgl);
}

void mkprf3HMMmag::codify_sites(Codon61GroupCount* rv, const int ix) {
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

double mkprf3HMMmag::hypergeometric1F1(const double a, const double b, const double c) const {
	gsl_sf_result res;
	gsl_sf_hyperg_1F1_e(a, b, c, &res);
	return res.val;
}

double mkprf3HMMmag::F(const double n, const double theta, const double gamma) const {
	double ret;
	if(fabs(gamma)<1.0e-3) ret = 1-(gsl_sf_beta(theta,n+theta+1)+gsl_sf_beta(n+theta,theta+1))/gsl_sf_beta(theta,theta+1);
	else ret = 1-(gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))+
				  gsl_sf_beta(n+theta,theta)*(1-hypergeometric1F1(n+theta,n+2*theta,-gamma)))/
		gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	return ret;
}

double mkprf3HMMmag::G(const double n, const double theta, const double gamma) const {
	double ret;
	if(fabs(gamma)<1.0e-3) ret = gsl_sf_beta(theta+1,n+theta)/gsl_sf_beta(theta+1,theta);
	else ret = gsl_sf_beta(theta,n+theta)*(1-hypergeometric1F1(theta,n+2*theta,-gamma))/
		gsl_sf_beta(theta,theta)/(1-hypergeometric1F1(theta,2*theta,-gamma));
	return ret;
}

void mkprf3HMMmag::summarize(ostream& out, Codon61GroupCount* rv, const int ix) {
	int pos, nNA=0;
	vector<int> n1(3,0), n2(3,0), n3(3,0), n4(3,0);
	const int L = rv->length();
	out << "mkprf3HMMmag summary for rv " << rv->name() << endl;
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

vector<double> mkprf3HMMmag::path_sampler(const RandomVariable* const rv, const int sp) {
	// Remember: there may be multiple child RVs!
	// COULD SCREW UP: if likelihood for one but not another child RV is requested,
	// and the path_sampler for that RV is then requested, it will return nonsense.
	// HACKED SOLUTION: cause error if multiple children are added (see add_random_variable)
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("mkprf3HMMmag::path_sampler(): RandomVariable not recognised");
	const int ix = it->second;
	
	// Use the alpha matrix
	Matrix<mydouble>& alpha = _alpha[ix][sp];
	Matrix<mydouble>& pneg = _pneg[ix][sp];
	const int L = alpha.nrows();
	const int ng = alpha.ncols();
	mydouble pi = 1.0/(double)ng;
	mydouble pswitch = get_p(sp)->get_double();
	mydouble pnoswitch = 1-pswitch;
	vector<double> gamma = get_gamma(sp)->get_doubles();
	vector<mydouble> pdraw(ng);
	vector<int> draw(L);
	vector<double> ret(L);
	int pos,g;
	{	// pos = L-1
		mydouble sumalpha = 0;
		for(g=0;g<ng;g++) {
			pdraw[g] = alpha[L-1][g];
			sumalpha += pdraw[g];
		}
		double U = _ran.U();
		for(g=0;g<ng;g++) {
			if((U -= (pdraw[g]/sumalpha).todouble())<=0.0) break;
		}
		if(g==ng) error("mkprf3HMMmag::path_sampler(): did not draw path correctly");
		draw[L-1] = g;
		// Negative or positive selection?
		const bool neg = _ran.bernoulliTF(pneg[L-1][g].todouble());
		ret[L-1] = (neg) ? -gamma[g] : gamma[g];
	}
	for(pos=L-2;pos>=0;pos--) {
		mydouble sumalpha = 0;
		for(g=0;g<ng;g++) {
			pdraw[g] = alpha[pos][g]*pswitch*pi;
			if(g==draw[pos+1]) pdraw[g] += alpha[pos][g]*pnoswitch;
			sumalpha += pdraw[g];
		}
		double U = _ran.U();
		for(g=0;g<ng;g++) {
			if((U -= (pdraw[g]/sumalpha).todouble())<=0.0) break;
		}
		if(g==ng) error("mkprf3HMMmag::path_sampler(): did not draw path correctly");
		draw[pos] = g;
		// Negative or positive selection?
		const bool neg = _ran.bernoulliTF(pneg[pos][g].todouble());
		ret[pos] = (neg) ? -gamma[g] : gamma[g];
	}
	
	return ret;
}
	
} // namespace gcat_omegaMap

