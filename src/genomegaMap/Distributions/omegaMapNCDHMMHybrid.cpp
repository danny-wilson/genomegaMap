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
 *  omegaMapNCDHMMHybridHybrid.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 02/08/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <DAG/DAG.h>
#include <omegaMap/Distributions/omegaMapNCDHMMHybrid.h>
#include <omegaMap/RandomVariables/Codon61Count.h>
#include <omegaMap/Utilities/omegaMapUtils.h>
#include <gsl/gsl_sf_hyperg.h>
#include <omegaMap/Utilities/mutation.h>
#include <gsl/gsl_cdf.h>
#include <DAG/RandomVariable.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string omegaMapNCDHMMHybridParameterNames[10] = {"anc","theta","kappa","gamma1","gamma2","T","p","pi","gamma1_wt","gamma2_wt"};

omegaMapNCDHMMHybrid::omegaMapNCDHMMHybrid(const int L, const int ng1, const int ng2, string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapNCDHMMHybrid"), Distribution(omegaMapNCDHMMHybridParameterNames,10), _L(L), _ng1(ng1), _ng2(ng2), _anc_changed(true),
_theta_changed(true), _kappa_changed(true), _gamma_changed(true), _T_changed(true), _p_changed(true), _pi_changed(true), _gamma_wt_changed(true), _update_phylo_rate_matrix(true), _update_popgen_rate_matrix(true), _update_pclik(true),
_update_pemiss(true), _update_alphabeta(true), _phylo_temp(Matrix<double>(61,61)), _phylo(vector< Matrix<mydouble> >(ng1+ng2,Matrix<mydouble>(61,61))), _phylo_Eigenvec(vector< Matrix<double> >(ng1+ng2,Matrix<double>(61,61))),
_phylo_Eigenval(vector< Vector<double> >(ng1+ng2,Vector<double>(61))), _previous_phylo(vector< Matrix<mydouble> >(ng1+ng2,Matrix<mydouble>(61,61))), _previous_phylo_Eigenvec(vector< Matrix<double> >(ng1+ng2,Matrix<double>(61,61))),
_previous_phylo_Eigenval(vector< Vector<double> >(ng1+ng2,Vector<double>(61))), _pdrm(Matrix<double>(61,61)), _pdrm_Eigenvec(Matrix<double>(61,61)), _pdrm_Eigenval(Vector<double>(61)) {
}

/*omegaMapNCDHMMHybrid::omegaMapNCDHMMHybrid(const omegaMapNCDHMMHybrid& x) : DAGcomponent(x), Distribution(x), _anc_changed(x._anc_changed),
 _theta_changed(x._theta_changed), _kappa_changed(x._kappa_changed), _gamma_changed(x._gamma_changed), _T_changed(x._T_changed), _p_changed(x._p_changed), _pi_changed(x._pi_changed), _gamma_wt_changed(x._gamma_wt_changed),
 _likpos(x._likpos), _previous_likpos(x._previous_likpos), _pemiss(x._pemiss), _previous_pemiss(x._previous_pemiss), _alphabeta(x._alphabeta), _previous_alphabeta(x._previous_alphabeta),
 _pclik(x._pclik), _previous_pclik(x._previous_pclik) {
 }*/

bool omegaMapNCDHMMHybrid::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61Count*>(random_variable));
	return false;
}

bool omegaMapNCDHMMHybrid::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:	{	//	anc
			// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
			Codon61SequenceRV* rv = dynamic_cast<Codon61SequenceRV*>(parameter);
			// Check the ancestral sequence length is compatible
			if(!(rv==0) && rv->length()!=seqlen()) {
				stringstream errMsg;
				errMsg << "omegaMapNCDHMMHybrid::check_parameter_type(): Codon61SequenceRV object " << rv->name() << " has sequence length " << rv->length();
				errMsg << " whereas omegaMapNCDHMMHybrid object " << name() << " has sequence length " << seqlen();
				error(errMsg.str().c_str());
			}
			return rv;
		}
		case 1:		//	theta
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 2:		//	kappa
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 3:	{	//	gamma1
			ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
			// Check the ancestral sequence length is compatible
			if(!(rv==0) && rv->length()!=ngamma1()) {
				stringstream errMsg;
				errMsg << "omegaMapNCDHMMHybrid::check_parameter_type(): ContinuousVectorVariable object " << rv->name() << " has " << rv->length();
				errMsg << " entries whereas omegaMapNCDHMMHybrid object " << name() << " requires " << ngamma1();
				error(errMsg.str().c_str());
			}
			return rv;
		}
		case 4:	{	//	gamma2
			ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
			// Check the ancestral sequence length is compatible
			if(!(rv==0) && rv->length()!=ngamma2()) {
				stringstream errMsg;
				errMsg << "omegaMapNCDHMMHybrid::check_parameter_type(): ContinuousVectorVariable object " << rv->name() << " has " << rv->length();
				errMsg << " entries whereas omegaMapNCDHMMHybrid object " << name() << " requires " << ngamma2();
				error(errMsg.str().c_str());
			}
			return rv;
		}
		case 5:		//	T
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 6:		//	p
			return(dynamic_cast<ContinuousVariable*>(parameter));
		case 7:		//	pi
			return(dynamic_cast<ContinuousVectorVariable*>(parameter));
		case 8:	{	//	gamma1_wt
			ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
			// Check the ancestral sequence length is compatible
			if(!(rv==0) && rv->length()!=ngamma1()) {
				stringstream errMsg;
				errMsg << "omegaMapNCDHMMHybrid::check_parameter_type(): ContinuousVectorVariable object " << rv->name() << " has " << rv->length();
				errMsg << " entries whereas omegaMapNCDHMMHybrid object " << name() << " requires " << ngamma1();
				error(errMsg.str().c_str());
			}
			return rv;
		}
		case 9:	{	//	gamma2_wt
			ContinuousVectorVariable* rv = dynamic_cast<ContinuousVectorVariable*>(parameter);
			// Check the ancestral sequence length is compatible
			if(!(rv==0) && rv->length()!=ngamma2()) {
				stringstream errMsg;
				errMsg << "omegaMapNCDHMMHybrid::check_parameter_type(): ContinuousVectorVariable object " << rv->name() << " has " << rv->length();
				errMsg << " entries whereas omegaMapNCDHMMHybrid object " << name() << " requires " << ngamma2();
				error(errMsg.str().c_str());
			}
			return rv;
		}
		default:
			error("omegaMapNCDHMMHybrid::check_parameter_type(): parameter not found");
	}
	return false;
}

int omegaMapNCDHMMHybrid::seqlen() const {
	return _L;
}

int omegaMapNCDHMMHybrid::ngamma1() const {
	return _ng1;
}

int omegaMapNCDHMMHybrid::ngamma2() const {
	return _ng2;
}

void omegaMapNCDHMMHybrid::set_anc(Codon61SequenceRV* anc) {
	set_parameter(0,(Variable*)anc);
}

void omegaMapNCDHMMHybrid::set_theta(ContinuousVariable* theta) {
	set_parameter(1,(Variable*)theta);
}

void omegaMapNCDHMMHybrid::set_kappa(ContinuousVariable* kappa) {
	set_parameter(2,(Variable*)kappa);
}

void omegaMapNCDHMMHybrid::set_gamma1(ContinuousVectorVariable* gamma) {
	set_parameter(3,(Variable*)gamma);
}

void omegaMapNCDHMMHybrid::set_gamma2(ContinuousVectorVariable* gamma) {
	set_parameter(4,(Variable*)gamma);
}

void omegaMapNCDHMMHybrid::set_T(ContinuousVariable* T) {
	set_parameter(5,(Variable*)T);
}

void omegaMapNCDHMMHybrid::set_p(ContinuousVariable* p) {
	set_parameter(6,(Variable*)p);
}

void omegaMapNCDHMMHybrid::set_pi(ContinuousVectorVariable* pi) {
	set_parameter(7,(Variable*)pi);
}

void omegaMapNCDHMMHybrid::set_gamma1_wt(ContinuousVectorVariable* gamma_wt) {
	set_parameter(8,(Variable*)gamma_wt);
}

void omegaMapNCDHMMHybrid::set_gamma2_wt(ContinuousVectorVariable* gamma_wt) {
	set_parameter(9,(Variable*)gamma_wt);
}

const Codon61SequenceRV* omegaMapNCDHMMHybrid::get_anc() const {
	return (const Codon61SequenceRV*)get_parameter(0);
}

const ContinuousVariable* omegaMapNCDHMMHybrid::get_theta() const {
	return (const ContinuousVariable*)get_parameter(1);
}

const ContinuousVariable* omegaMapNCDHMMHybrid::get_kappa() const {
	return (const ContinuousVariable*)get_parameter(2);
}

const ContinuousVectorVariable* omegaMapNCDHMMHybrid::get_gamma1() const {
	return (const ContinuousVectorVariable*)get_parameter(3);
}

const ContinuousVectorVariable* omegaMapNCDHMMHybrid::get_gamma2() const {
	return (const ContinuousVectorVariable*)get_parameter(4);
}

const ContinuousVariable* omegaMapNCDHMMHybrid::get_T() const {
	return (const ContinuousVariable*)get_parameter(5);
}

const ContinuousVariable* omegaMapNCDHMMHybrid::get_p() const {
	return (const ContinuousVariable*)get_parameter(6);
}

const ContinuousVectorVariable* omegaMapNCDHMMHybrid::get_pi() const {
	return (const ContinuousVectorVariable*)get_parameter(7);
}

const ContinuousVectorVariable* omegaMapNCDHMMHybrid::get_gamma1_wt() const {
	return (const ContinuousVectorVariable*)get_parameter(8);
}

const ContinuousVectorVariable* omegaMapNCDHMMHybrid::get_gamma2_wt() const {
	return (const ContinuousVectorVariable*)get_parameter(9);
}

mydouble omegaMapNCDHMMHybrid::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("omegaMapNCDHMMHybrid::likelihood(): RandomVariable not recognised");
	const int ix = it->second;
	
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
	_update_phylo_rate_matrix = _update_popgen_rate_matrix = _update_pclik = _update_pemiss = _update_alphabeta = true;
#endif
	if(_update_phylo_rate_matrix) update_phylo_rate_matrix();
	if(_update_popgen_rate_matrix) update_popgen_rate_matrix();
	if(_update_pclik) update_pclik(val,ix);
	if(_update_pemiss) update_pemiss(ix);
	if(_update_alphabeta) update_alphabeta(ix);
	// Recalculate the final likelihood
	return HMMlik(ix);
}

void omegaMapNCDHMMHybrid::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
	Codon61Count* rv = dynamic_cast<Codon61Count*>(random_variable);
	// Check the random variable has a compatible sequence length
	if(rv->length()!=seqlen()) {
		stringstream errMsg;
		errMsg << "omegaMapNCDHMMHybrid::add_random_variable(): Codon61Count object " << rv->name() << " has sequence length " << rv->length();
		errMsg << " whereas omegaMapNCDHMMHybrid object " << name() << " has sequence length " << seqlen();
		error(errMsg.str().c_str());
	}
	// Memory allocation
	_likpos.push_back(recentre());
	_pemiss.push_back(Matrix<mydouble>(_L,_ng1+_ng2));
	_alphabeta.push_back(Matrix<mydouble>(_L+1,_ng1));
	_pclik.push_back(vector< Matrix<mydouble> >(_L,Matrix<mydouble>(_ng1+_ng2,61)));
	_previous_likpos.push_back(recentre());
	_previous_pemiss.push_back(Matrix<mydouble>(_L,_ng1+_ng2));
	_previous_alphabeta.push_back(Matrix<mydouble>(_L+1,_ng1));
	_previous_pclik.push_back(vector< Matrix<mydouble> >(_L,Matrix<mydouble>(_ng1+_ng2,61)));
}

void omegaMapNCDHMMHybrid::remove_random_variable(RandomVariable* random_variable) {
	error("omegaMapNCDHMMHybrid::remove_random_variable(): not permitted");
}

void omegaMapNCDHMMHybrid::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	if(sgl==Variable::_ACCEPT || sgl==Variable::_REVERT) {
		_anc_changed = _theta_changed = _kappa_changed = _gamma_changed = _T_changed = _p_changed = _pi_changed = _gamma_wt_changed = false;
		_update_phylo_rate_matrix = _update_popgen_rate_matrix = _update_pclik = _update_pemiss = _update_alphabeta = false;
	}
	if(v==(const Value*)get_anc()) {
		if(sgl==Variable::_SET) {
			_anc_changed = true;
			_previous_likpos = _likpos;
		}
		else if(sgl==Variable::_PROPOSE) {
			_anc_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			//			_previous_pclik = _pclik;
			//			_previous_phylo = _phylo;
			//			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			//			_previous_phylo_Eigenval = _phylo_Eigenval;
			//			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			//			_pclik = _previous_pclik;
			//			_phylo = _previous_phylo;
			//			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			//			_phylo_Eigenval = _previous_phylo_Eigenval;
			//			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_theta()) {
		if(sgl==Variable::_SET) {
			_theta_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_theta_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			_previous_pclik = _pclik;
			_previous_phylo = _phylo;
			//			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			//			_previous_phylo_Eigenval = _phylo_Eigenval;
			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			_pclik = _previous_pclik;
			_phylo = _previous_phylo;
			//			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			//			_phylo_Eigenval = _previous_phylo_Eigenval;
			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_kappa()) {
		if(sgl==Variable::_SET) {
			_kappa_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_kappa_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			_previous_pclik = _pclik;
			_previous_phylo = _phylo;
			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			_previous_phylo_Eigenval = _phylo_Eigenval;
			_previous_pdrm = _pdrm;
			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			_pclik = _previous_pclik;
			_phylo = _previous_phylo;
			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			_phylo_Eigenval = _previous_phylo_Eigenval;
			_pdrm = _previous_pdrm;
			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_gamma1() || v==(const Value*)get_gamma2()) {
		if(sgl==Variable::_SET) {
			_gamma_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_gamma_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			_previous_pclik = _pclik;
			_previous_phylo = _phylo;
			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			_previous_phylo_Eigenval = _phylo_Eigenval;
			//			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			_pclik = _previous_pclik;
			_phylo = _previous_phylo;
			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			_phylo_Eigenval = _previous_phylo_Eigenval;
			//			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_T()) {
		if(sgl==Variable::_SET) {
			_T_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_T_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			//			_previous_pclik = _pclik;
			_previous_phylo = _phylo;
			//			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			//			_previous_phylo_Eigenval = _phylo_Eigenval;
			//			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			//			_pclik = _previous_pclik;
			_phylo = _previous_phylo;
			//			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			//			_phylo_Eigenval = _previous_phylo_Eigenval;
			//			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_p()) {
		if(sgl==Variable::_SET) {
			_p_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_p_changed = true;
			_previous_likpos = _likpos;
			//			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			//			_previous_pclik = _pclik;
			//			_previous_phylo = _phylo;
			//			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			//			_previous_phylo_Eigenval = _phylo_Eigenval;
			//			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			//			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			//			_pclik = _previous_pclik;
			//			_phylo = _previous_phylo;
			//			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			//			_phylo_Eigenval = _previous_phylo_Eigenval;
			//			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_pi()) {
		if(sgl==Variable::_SET) {
			_pi_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_pi_changed = true;
			_previous_likpos = _likpos;
			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			_previous_pclik = _pclik;
			_previous_phylo = _phylo;
			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			_previous_phylo_Eigenval = _phylo_Eigenval;
			_previous_pdrm = _pdrm;
			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			_pclik = _previous_pclik;
			_phylo = _previous_phylo;
			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			_phylo_Eigenval = _previous_phylo_Eigenval;
			_pdrm = _previous_pdrm;
			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	if(v==(const Value*)get_gamma1_wt() || v==(const Value*)get_gamma2_wt()) {
		if(sgl==Variable::_SET) {
			_gamma_wt_changed = true;
		}
		else if(sgl==Variable::_PROPOSE) {
			_gamma_wt_changed = true;
			_previous_likpos = _likpos;
			//			_previous_pemiss = _pemiss;
			_previous_alphabeta = _alphabeta;
			//			_previous_pclik = _pclik;
			//			_previous_phylo = _phylo;
			//			_previous_phylo_Eigenvec = _phylo_Eigenvec;
			//			_previous_phylo_Eigenval = _phylo_Eigenval;
			//			_previous_pdrm = _pdrm;
			//			_previous_pdrm_Eigenvec = _pdrm_Eigenvec;
			//			_previous_pdrm_Eigenval = _pdrm_Eigenval;
		}
		else if(sgl==Variable::_ACCEPT) {
		}
		else if(sgl==Variable::_REVERT) {
			_likpos = _previous_likpos;
			//			_pemiss = _previous_pemiss;
			_alphabeta = _previous_alphabeta;
			//			_pclik = _previous_pclik;
			//			_phylo = _previous_phylo;
			//			_phylo_Eigenvec = _previous_phylo_Eigenvec;
			//			_phylo_Eigenval = _previous_phylo_Eigenval;
			//			_pdrm = _previous_pdrm;
			//			_pdrm_Eigenvec = _previous_pdrm_Eigenvec;
			//			_pdrm_Eigenval = _previous_pdrm_Eigenval;
		}
	}
	// If kappa, gamma or pi have changed, update the phylogenetic rate matrix
	if(_theta_changed || _kappa_changed || _gamma_changed || _pi_changed || _T_changed) _update_phylo_rate_matrix = true;
	// If theta, kappa or pi have changed, update the population genetic parent-dependent rate matrix
	if(_theta_changed|| _kappa_changed || _pi_changed) _update_popgen_rate_matrix = true;
	// If theta, kappa, gamma or pi have changed, update the ancestor-dependent (phylogenetically conditioned) likelihoods
	if(_theta_changed || _kappa_changed || _gamma_changed || _pi_changed) _update_pclik = true;
	// If anc, theta, kappa, gamma, T or pi have changed, update the emission probabilities
	if(_anc_changed || _theta_changed || _kappa_changed || _gamma_changed || _T_changed || _pi_changed) _update_pemiss = true;
	// If anc, theta, kappa, gamma, T, p or pi have changed, update the forward-backward probabilities
	if(_anc_changed || _theta_changed || _kappa_changed || _gamma_changed || _T_changed || _p_changed || _pi_changed || _gamma_wt_changed) _update_alphabeta = true;
	Distribution::receive_signal_from_parent(v,sgl);
}

void omegaMapNCDHMMHybrid::update_phylo_rate_matrix() {
	const ContinuousVectorVariable& _gamma1 = *get_gamma1();
	const ContinuousVectorVariable& _gamma2 = *get_gamma2();
	const double thetaT = get_theta()->get_double()*get_T()->get_double();
	const double kappa = get_kappa()->get_double();
	const vector<double>& pi = get_pi()->get_doubles();
	Vector<double> Pi(61);
	int i;
	for(i=0;i<61;i++) Pi[i] = pi[i];
	int g;
	for(g=0;g<_ng1;g++) {
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
		{
#else
			if(_kappa_changed || _gamma_changed || _pi_changed) {
#endif
				const double gamma = _gamma1.get_double(g);
				const double omega = (fabs(gamma)<1e-3) ? 1.0 : gamma/(1-exp(-gamma));
				// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
				NY98_61::diagonalize_symmetric(_phylo_Eigenvec[g],_phylo_Eigenval[g],kappa,omega,pi);
			}
			// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
			// Therefore thetaT is defined as twice the expected number of mutations under neutrality
			NY98_61::build_Pt(_phylo_temp,_phylo_Eigenvec[g],_phylo_Eigenval[g],Pi,thetaT*61.0/6.0/(2.0+kappa));
			// Convert to mydouble
			for(i=0;i<61;i++) {
				int j;
				for(j=0;j<61;j++) {
					_phylo[g][i][j] = mydouble(_phylo_temp[i][j]);
				}
			}
		}
		int gprime = g;
		for(g=0;g<_ng2;g++,gprime++) {
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
			{
#else
				if(_kappa_changed || _gamma_changed || _pi_changed) {
#endif
					const double gamma = _gamma2.get_double(g);
					const double omega = (fabs(gamma)<1e-3) ? 1.0 : gamma/(1-exp(-gamma));
					// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
					NY98_61::diagonalize_symmetric(_phylo_Eigenvec[gprime],_phylo_Eigenval[gprime],kappa,omega,pi);
				}
				// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
				// Therefore thetaT is defined as twice the expected number of mutations under neutrality
				NY98_61::build_Pt(_phylo_temp,_phylo_Eigenvec[gprime],_phylo_Eigenval[gprime],Pi,thetaT*61.0/6.0/(2.0+kappa));
				// Convert to mydouble
				for(i=0;i<61;i++) {
					int j;
					for(j=0;j<61;j++) {
						_phylo[gprime][i][j] = mydouble(_phylo_temp[i][j]);
					}
				}
			}	
			_update_phylo_rate_matrix = false;
	}
	
	void omegaMapNCDHMMHybrid::update_popgen_rate_matrix() {
		const double theta = get_theta()->get_double();
		const double kappa = get_kappa()->get_double();
		const vector<double>& pi = get_pi()->get_doubles();
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
		{
#else
			if(_kappa_changed || _pi_changed) {
#endif
				// Build an (unscaled) symmetric version of the rate matrix and diagonalize it
				NY98_61::diagonalize_symmetric(_pdrm_Eigenvec,_pdrm_Eigenval,kappa,1.0,pi);
			}
			Vector<double> Pi(61);
			int i,j;
			for(i=0;i<61;i++) Pi[i] = pi[i];
			// Under neutrality, the expected rate from diagonalize_symmetric is 3*(2+kappa)/61
			// Therefore theta is defined as twice the expected rate under neutrality
			NY98_61::build_Rt(_pdrm,_pdrm_Eigenvec,_pdrm_Eigenval,Pi,theta*61.0/6.0/(2.0+kappa),1.0);
			// Sweep through to build final approximate parent-independent mutation rate matrix
			for(i=0;i<61;i++) {
				double rii = _pdrm[i][i];
				for(j=0;j<61;j++) {
					if(j!=i) _pdrm[i][j]/=rii;
				}
				_pdrm[i][i] = 1.0e-6; //numeric_limits<double>::min();
			}		
			_update_popgen_rate_matrix = false;
		}
		
		void omegaMapNCDHMMHybrid::update_pclik(const Value* val, const int ix) {
			// Calculate the phylogenetically conditioned likelihood for each site given
			// the population ancestral allele and the selection coefficient gamma
			bool full_recalc = (_theta_changed || _kappa_changed || _gamma_changed || _pi_changed);
			vector< Matrix<mydouble> > &pclik = _pclik[ix];
			const Codon61Count& ct = *((const Codon61Count*)val);
			const Codon61SequenceRV& _anc = *get_anc();
			const ContinuousVectorVariable& gamma1 = *get_gamma1();
			const ContinuousVectorVariable& gamma2 = *get_gamma2();
			const double N = (double)(ct.n());
			int pos;
			for(pos=0;pos<_L;pos++) {
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
				{
#else
					if(full_recalc || _anc.has_changed(pos)) {
#endif
						// Precalculate common factors in the likelihood at this site (i.e. irrespective of ancestral codon)
						// Remove indels from the sample size (i.e. treat as missing data)
						const double Npos = N-(double)ct[pos][61];
						double phi = lgamma(Npos+1.);
						int i;
						for(i=0;i<61;i++) {
							const int cti = ct[pos][i];
							phi -= lgamma((double)(cti+1));
						}
						
						int anc, g;
						// For each population ancestral allele
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
								const double cti = (double)ct[pos][i];
								const double muti = _pdrm[anc][i];
								ttheta += muti;
								if(omegaMapUtils::oneLetterCodes[i]==aAA) {
									pt1 += cti;
									ptheta += muti;
								}
								if(i==anc) likposanc += lgamma(cti+muti+1.0) - lgamma(muti+1.0);
								else likposanc += lgamma(cti+muti) - lgamma(muti);
							}
							likposanc += lgamma(ttheta) - lgamma(Npos+ttheta) + log(ptheta) - log(pt1+ptheta);
							// Selection coefficient
							for(g=0;g<_ng1;g++) {
								double likposancg = 0;
								const double sigma = gamma1.get_double(g);
								if(fabs(sigma)<1e-3) {
									// Ratio in the limit that sigma equals zero
									likposancg += log(pt1+ptheta) - log(Npos+ttheta) - log(ptheta) + log(ttheta);
								}
								else {
									const double hpnum = (omegaMapUtils::hypergeometric1F1(pt1+ptheta,Npos+ttheta,-sigma)-1.0);
									const double hpden = (omegaMapUtils::hypergeometric1F1(ptheta,ttheta,-sigma)-1.0);
									if(hpden==0.0 || hpnum/hpden<=0.0) {
										likposancg += log(pt1+ptheta) - log(Npos+ttheta) - log(ptheta) + log(ttheta);
									}
									else {
										likposancg += log(hpnum/hpden);
									}
								}
								// Store the within-population conditional likelihood (conditional on ancestral codon)
								pclik[pos][g][anc].setlog(phi+likposanc+likposancg);
							}
							int gprime = g;
							for(g=0;g<_ng2;g++,gprime++) {
								double likposancg = 0;
								const double sigma = gamma2.get_double(g);
								if(fabs(sigma)<1e-3) {
									// Ratio in the limit that sigma equals zero
									likposancg += log(pt1+ptheta) - log(Npos+ttheta) - log(ptheta) + log(ttheta);
								}
								else {
									const double hpnum = (omegaMapUtils::hypergeometric1F1(pt1+ptheta,Npos+ttheta,-sigma)-1.0);
									const double hpden = (omegaMapUtils::hypergeometric1F1(ptheta,ttheta,-sigma)-1.0);
									if(hpden==0.0 || hpnum/hpden<=0.0) {
										likposancg += log(pt1+ptheta) - log(Npos+ttheta) - log(ptheta) + log(ttheta);
									}
									else {
										likposancg += log(hpnum/hpden);
									}
								}
								// Store the within-population conditional likelihood (conditional on ancestral codon)
								pclik[pos][gprime][anc].setlog(phi+likposanc+likposancg);
							}
						}
					}
				}	
			}
			
			void omegaMapNCDHMMHybrid::update_pemiss(const int ix) {
				// Conditional on the ancestor, sum the phylogenetically conditioned
				// likelihoods, weighted by the phylogenetic transition probabilities.
				bool full_recalc = (_theta_changed || _kappa_changed || _gamma_changed || _pi_changed || _T_changed);
				const Codon61SequenceRV& _anc = *get_anc();
				Matrix<mydouble> &pemiss = _pemiss[ix];
				vector< Matrix<mydouble> > &pclik = _pclik[ix];
				int dec, pos, g;
				for(pos=0;pos<_L;pos++) {
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
					{
#else
						if(full_recalc || _anc.has_changed(pos)) {
#endif
							const int anc = _anc[pos];
							for(g=0;g<_ng1+_ng2;g++) {
								pemiss[pos][g] = 0;
								for(dec=0;dec<61;dec++) {
									pemiss[pos][g] += _phylo[g][anc][dec] * pclik[pos][g][dec];
								}
							}
						}
					}
				}
				
				void omegaMapNCDHMMHybrid::update_alphabeta(const int ix) {
					// Assuming pemiss is up-to-date, recalculate the forward-backward probabilities as necessary
					const Codon61SequenceRV& anc = *get_anc();
					int& likpos = _likpos[ix];
					int& previous_likpos = _previous_likpos[ix];
					bool full_recalc = (_theta_changed || _kappa_changed || _gamma_changed || _T_changed || _pi_changed || _p_changed || _gamma_wt_changed);
#ifdef _DEBUG_OMEGAMAP_NCD_HMM
					full_recalc = true;
#endif
					// Test to see whether a single site in _anc has changed
					if(!full_recalc) {
						int wh = -1;
						int i;
						for(i=0;i<anc.length();i++) {
							if(anc.has_changed(i)) {
								if(wh!=-1) {
									full_recalc = true;
									break;
								}
								wh = i;
							}
						}
						if(!full_recalc) {
							if(wh==-1) {
								// Changed to the same thing, so don't change likpos
							}
							else {
								likpos = wh;
							}
						}
					}
					if(full_recalc) {
						// Forward to centre, backward to centre. Set likpos = centre
						likpos = recentre();
						forward(ix,0);
						backward(ix,_L-1);
					}
					else {
						// Only if only _anc_changed and at a single site
						if(previous_likpos<likpos) forward(ix,previous_likpos+1);
						else if(previous_likpos==likpos) forward(ix,likpos);
						else if(likpos<previous_likpos) {
							forward(ix,likpos);
							backward(ix,previous_likpos-1);
						}
					}
				}
				
				int omegaMapNCDHMMHybrid::recentre() const {
					return (_L-1)/2;
				}
				
				void omegaMapNCDHMMHybrid::forward(const int ix, const int left) {
					Matrix<mydouble> &alpha = _alphabeta[ix];
					Matrix<mydouble> &pemiss = _pemiss[ix];
					mydouble pswitch = get_p()->get_double();
					mydouble pnoswitch = 1-pswitch;
					int g, pos = left;
					vector<mydouble> pg1(_ng1);
					vector<mydouble> pg(_ng1+_ng2);
					mydouble pg1tot = 0, pgtot = 0;
					for(g=0;g<_ng1;g++) {
						pg1[g] = mydouble(get_gamma1_wt()->get_double(g));
						pg[g] = pg1[g];
						pg1tot += pg1[g];
						pgtot += pg[g];
					}
					int gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] = mydouble(get_gamma2_wt()->get_double(g));
						pgtot += pg[gprime];
					}
					mydouble sumpg1 = 0;
					for(g=0;g<_ng1;g++) {
						pg1[g] /= pg1tot;
						pg[g] /= pgtot;
						sumpg1 += pg[g];
					}
					mydouble m = sumpg1;
					gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] /= pgtot;
					}
					mydouble sumalpha = 0.0;
					if(pos==0) {
						mydouble c = pg[_ng1]*pemiss[0][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[0][_ng1+g];
						for(g=0;g<_ng1;g++) {
							alpha[0][g] = (m*pemiss[0][g]+c)*pg1[g];
							sumalpha += alpha[0][g];
						}
						sumalpha *= pswitch;
						pos++;
					}
					else {
						for(g=0;g<_ng1;g++) {
							sumalpha += alpha[pos-1][g];
						}
						sumalpha *= pswitch;
					}
					for(;pos<=_likpos[ix];pos++) {
						mydouble newsumalpha = 0.0;
						mydouble c = pg[_ng1]*pemiss[pos][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[pos][_ng1+g];
						for(g=0;g<_ng1;g++) {
							alpha[pos][g] = (sumalpha*pg1[g]+pnoswitch*alpha[pos-1][g])*(m*pemiss[pos][g]+c);
							newsumalpha += alpha[pos][g];
						}
						newsumalpha *= pswitch;
						sumalpha = newsumalpha;
					}
				}
				
				void omegaMapNCDHMMHybrid::backward(const int ix, const int right) {
					Matrix<mydouble> &beta = _alphabeta[ix];
					Matrix<mydouble> &pemiss = _pemiss[ix];
					mydouble pswitch = get_p()->get_double();
					mydouble pnoswitch = 1-pswitch;
					int g, pos = right;
					vector<mydouble> pg1(_ng1);
					vector<mydouble> pg(_ng1+_ng2);
					mydouble pg1tot = 0, pgtot = 0;
					for(g=0;g<_ng1;g++) {
						pg1[g] = mydouble(get_gamma1_wt()->get_double(g));
						pg[g] = pg1[g];
						pg1tot += pg1[g];
						pgtot += pg[g];
					}
					int gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] = mydouble(get_gamma2_wt()->get_double(g));
						pgtot += pg[gprime];
					}
					mydouble sumpg1 = 0;
					for(g=0;g<_ng1;g++) {
						pg1[g] /= pg1tot;
						pg[g] /= pgtot;
						sumpg1 += pg[g];
					}
					mydouble m = sumpg1;
					gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] /= pgtot;
					}
					mydouble sumbeta = 0.0;
					if(pos==_L-1) {
						mydouble c = pg[_ng1]*pemiss[pos][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[pos][_ng1+g];
						for(g=0;g<_ng1;g++) {
							beta[_L][g] = 1.0;
							sumbeta += pg1[g]*(m*pemiss[pos][g]+c);
						}
						sumbeta *= pswitch;
						pos--;
					}
					else {
						mydouble c = pg[_ng1]*pemiss[pos+1][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[pos+1][_ng1+g];
						for(g=0;g<_ng1;g++) {
							sumbeta += pg1[g]*(m*pemiss[pos+1][g]+c)*beta[pos+2][g];
						}
						sumbeta *= pswitch;
					}
					for(;pos>=_likpos[ix];pos--) {
						mydouble newsumbeta = 0.0;
						mydouble c = pg[_ng1]*pemiss[pos+1][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[pos+1][_ng1+g];
						mydouble c2 = pg[_ng1]*pemiss[pos][_ng1];
						for(g=1;g<_ng2;g++) c2 += pg[_ng1+g]*pemiss[pos][_ng1+g];
						for(g=0;g<_ng1;g++) {
							beta[pos+1][g] = sumbeta+pnoswitch*(m*pemiss[pos+1][g]+c)*beta[pos+2][g];
							newsumbeta += pg1[g]*(m*pemiss[pos][g]+c2)*beta[pos+1][g];
						}
						newsumbeta *= pswitch;
						sumbeta = newsumbeta;
					}
				}
				
				mydouble omegaMapNCDHMMHybrid::HMMlik(const int ix) {
					mydouble lik = 0.0;
					int g;
					for(g=0;g<_ng1;g++) {
						lik += _alphabeta[ix][_likpos[ix]][g]*_alphabeta[ix][_likpos[ix]+1][g];
					}
					return lik;
				}
				
				vector<double> omegaMapNCDHMMHybrid::path_sampler(const RandomVariable* const rv) {
					// Remember: a restriction is added in add_random_variable that only 1 rv is allowed
					// Otherwise COULD SCREW UP: if likelihood for one but not another child RV is requested,
					// and the path_sampler for that RV is then requested, it will return nonsense.
					// SOLUTION: aren't all likelihoods called by the DAG to avoid this and other problems?
					map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
					if(it==_rv_index.end()) error("omegaMapNCDHMMHybrid::path_sampler(): RandomVariable not recognised");
					const int ix = it->second;
					
					// Use the forward-backward probabilities
					Matrix<mydouble>& alphabeta = _alphabeta[ix];
					Matrix<mydouble> &pemiss = _pemiss[ix];
					vector<mydouble> pg1(_ng1);
					vector<mydouble> pg(_ng1+_ng2);
					mydouble pg1tot = 0, pgtot = 0;
					int g;
					for(g=0;g<_ng1;g++) {
						pg1[g] = mydouble(get_gamma1_wt()->get_double(g));
						pg[g] = pg1[g];
						pg1tot += pg1[g];
						pgtot += pg[g];
					}
					int gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] = mydouble(get_gamma2_wt()->get_double(g));
						pgtot += pg[gprime];
					}
					mydouble sumpg1 = 0;
					for(g=0;g<_ng1;g++) {
						pg1[g] /= pg1tot;
						pg[g] /= pgtot;
						sumpg1 += pg[g];
					}
					mydouble m = sumpg1;
					gprime = g;
					for(g=0;g<_ng2;g++,gprime++) {
						pg[gprime] /= pgtot;
					}
					mydouble pswitch = get_p()->get_double();
					mydouble pnoswitch = 1-pswitch;
					//	mydouble pswitchpg = pswitch*pg;
					vector<double> gamma1 = get_gamma1()->get_doubles();
					vector<double> gamma2 = get_gamma2()->get_doubles();
					vector<mydouble> pdraw(_ng1), qdraw(_ng2+1);
					vector<int> draw(_L);
					vector<double> ret(_L);
					int pos = _likpos[ix];
					// Start at _likpos
					mydouble denom = 0;
					for(g=0;g<_ng1;g++) {
						pdraw[g] = alphabeta[pos][g]*alphabeta[pos+1][g];
						denom += pdraw[g];
					}
					double U = _ran.U();
					for(g=0;g<_ng1;g++) {
						if((U -= (pdraw[g]/denom).todouble())<=0.0) break;
					}
					if(g==_ng1) error("omegaMapNCDHMMHybrid::path_sampler(): did not draw path correctly");
					draw[pos] = g;
					ret[pos] = final_draw(pos,g,-U*(denom/pdraw[g]).todouble(),m,gamma1,gamma2,pg,qdraw,pemiss);
//					ret[pos] = final_draw(pos,g,_ran.U(),m,gamma1,gamma2,pg,qdraw,pemiss);
					// Work leftwards
					for(pos=_likpos[ix]-1;pos>=0;pos--) {
						denom = 0;
						for(g=0;g<_ng1;g++) {
							pdraw[g] = alphabeta[pos][g]*pswitch*pg1[draw[pos+1]];
							if(g==draw[pos+1]) pdraw[g] += alphabeta[pos][g]*pnoswitch;
							denom += pdraw[g];
						}
						U = _ran.U();
						for(g=0;g<_ng1;g++) {
							if((U -= (pdraw[g]/denom).todouble())<=0.0) break;
						}
						if(g==_ng1) error("omegaMapNCDHMMHybrid::path_sampler(): did not draw path correctly");
						draw[pos] = g;
						ret[pos] = final_draw(pos,g,-U*(denom/pdraw[g]).todouble(),m,gamma1,gamma2,pg,qdraw,pemiss);
//						ret[pos] = final_draw(pos,g,_ran.U(),m,gamma1,gamma2,pg,qdraw,pemiss);
					}
					// Work rightwards
					for(pos=_likpos[ix]+1;pos<_L;pos++) {
						mydouble c = pg[_ng1]*pemiss[pos][_ng1];
						for(g=1;g<_ng2;g++) c += pg[_ng1+g]*pemiss[pos][_ng1+g];
						denom = 0;
						for(g=0;g<_ng1;g++) {
							pdraw[g] = alphabeta[pos+1][g]*pswitch*pg1[g];
							if(g==draw[pos-1]) pdraw[g] += alphabeta[pos+1][g]*pnoswitch;
							// Cf. this next line to working leftwards.
							pdraw[g] *= (m*pemiss[pos][g]+c);
							denom += pdraw[g];
						}
						U = _ran.U();
						for(g=0;g<_ng1;g++) {
							if((U -= (pdraw[g]/denom).todouble())<=0.0) break;
						}
						if(g==_ng1) error("omegaMapNCDHMMHybrid::path_sampler(): did not draw path correctly");
						draw[pos] = g;
						ret[pos] = final_draw(pos,g,-U*(denom/pdraw[g]).todouble(),m,gamma1,gamma2,pg,qdraw,pemiss);
//						ret[pos] = final_draw(pos,g,_ran.U(),m,gamma1,gamma2,pg,qdraw,pemiss);
					}
					return ret;
				}
				
				double omegaMapNCDHMMHybrid::final_draw(const int pos, const int draw, const double U_in, const mydouble m, const vector<double> &gamma1, const vector<double> &gamma2, const vector<mydouble> &pg, vector<mydouble> qdraw, const Matrix<mydouble> &pemiss) {
					// To revert to previous behaviour:
					// return gamma1[draw];
					mydouble denom = 0;
					int h;
					for(h=0;h<_ng2;h++) {
						qdraw[h] = pg[_ng1+h]*pemiss[pos][_ng1+h];
						denom += qdraw[h];
					}
					qdraw[_ng2] = m*pemiss[pos][draw];
					denom += qdraw[_ng2];
					double U = U_in;
					for(h=0;h<_ng2+1;h++) {
						if((U -= (qdraw[h]/denom).todouble())<=0.0) break;
					}
					if(h==_ng2+1) error("omegaMapNCDHMMHybrid::final_draw(): did not draw path correctly");
					if(h==_ng2) return gamma1[draw];
					return gamma2[h];
				}
				
				const string omegaMapNCDHMMHybridPathSamplerParameterNames[0];
				
				omegaMapNCDHMMHybridPathSampler::omegaMapNCDHMMHybridPathSampler(string rv_name, string distribution_name, string name, DAG* dag) : DAGcomponent(name,dag,"omegaMapNCDHMMHybridPathSampler"), Transformation(omegaMapNCDHMMHybridPathSamplerParameterNames,0), _distribution_name(distribution_name) {
					_rv = getDAG()->get_random_variable(rv_name);
				}
				
				omegaMapNCDHMMHybridPathSampler::omegaMapNCDHMMHybridPathSampler(const omegaMapNCDHMMHybridPathSampler& x) : DAGcomponent(x), Transformation(x), _distribution_name(x._distribution_name), _hmm(x._hmm), _rv(x._rv), _L(x._L) {
				}
				
				string omegaMapNCDHMMHybridPathSampler::validate() const {
					// For starters, do not allow this to parameterize anything else
					if(n_child_distributions()+n_child_transformations()>0) {
						string errTxt = "omegaMapNCDHMMHybridPathSampler: object " + name() + " may not parameterize other objects";
						return errTxt;
					}
					// Secondly, finish identifying the distribution to which it relates (Distributions
					// are initialized after Transformations, so this could not be done in the constructor)
					_hmm = dynamic_cast<omegaMapNCDHMMHybrid*>(getDAG()->get_distribution(_distribution_name));
					if(!_hmm) {
						string errTxt = "omegaMapNCDHMMHybridPathSampler: object " + name() + " cannot find omegaMapNCDHMMHybrid object " + _distribution_name;
						return errTxt;
					}
					const Codon61Count* ct = dynamic_cast<const Codon61Count*>(_rv);
					if(!ct) {
						string errTxt = "omegaMapNCDHMMHybridPathSampler: " + _rv->name() + " is not of type Codon61Count";
						return errTxt;
					}
					_L = ct->length();
					return "";
				}
				
				int omegaMapNCDHMMHybridPathSampler::length() const {
					return _L;
				}
				
				// NB: Repeatedly calling get_double() gives marginal draws, whereas get_doubles() draws a complete path
				// For that reason, the print() method is overwritten
				double omegaMapNCDHMMHybridPathSampler::get_double(const int i) const {
					if(i<0) error("omegaMapNCDHMMHybridPathSampler::get_double(): index cannot be negative");
					vector<double> x = get_doubles();
					if(i>=_L) error("omegaMapNCDHMMHybridPathSampler::get_double(): index too large");
					return x[i];
				}
				
				vector<double> omegaMapNCDHMMHybridPathSampler::get_doubles() const {
					return _hmm->path_sampler(_rv);
				}
				
				bool omegaMapNCDHMMHybridPathSampler::has_changed(const int i) const {
					return true;
				}
				
				vector<bool> omegaMapNCDHMMHybridPathSampler::has_changed() const {
					return vector<bool>(length(),true);
				}
				
				bool omegaMapNCDHMMHybridPathSampler::check_parameter_type(const int i, Variable* parameter) {
					switch(i) {
						default:
							error("omegaMapNCDHMMHybridPathSampler::check_parameter_type(): parameter not found");
					}
					return false;
				}
				
				void omegaMapNCDHMMHybridPathSampler::print(ostream& out, string sep) {
					int i;
					vector<double> x = get_doubles();
					for(i=0;i<length();i++) {
						if(i>0) out << sep;
						out << x[i];
					}
				}
				
			} // namespace gcat_omegaMap
