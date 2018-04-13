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
 *  phylogeny.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 11/25/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/Distributions/phylogeny.h>
#include <omegaMap/RandomVariables/Codon61Alignment.h>
#include <omegaMap/Utilities/omegaMapUtils.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const string PhylogenyParameterNames[1] = {"phylo"};

Phylogeny::Phylogeny(string name, DAG* dag) : DAGcomponent(name,dag,"Phylogeny"), Distribution(PhylogenyParameterNames,1) {
}

Phylogeny::Phylogeny(const Phylogeny& x) : DAGcomponent(x), Distribution(x) {
}

bool Phylogeny::check_random_variable_type(RandomVariable* random_variable) {
	// Insist on the particular RandomVariable (not just Variable) as it also specifies the encoding
	return(dynamic_cast<Codon61Alignment*>(random_variable));
	return false;
}

bool Phylogeny::check_parameter_type(const int i, Variable* parameter) {
	switch(i) {
		case 0:
			return(dynamic_cast<NY98_TransProbMosaic*>(parameter));
		default:
			error("Phylogeny::check_parameter_type(): parameter not found");
	}
	return false;
}

void Phylogeny::set_phylo(NY98_TransProbMosaic* phylo) {
	set_parameter(0,(Variable*)phylo);
}

const NY98_TransProbMosaic* Phylogeny::get_phylo() const {
	return (const NY98_TransProbMosaic*)get_parameter(0);
}

mydouble Phylogeny::likelihood(const RandomVariable* rv, const Value* val) {
	// Remember: there may be multiple child RVs!
	map< const RandomVariable*, int >::iterator it = _rv_index.find(rv);
	if(it==_rv_index.end()) error("Phylogeny::likelihood(): RandomVariable not recognised");
	
	const Codon61Alignment& ct = *((const Codon61Alignment*)val);
	if(ct.n()!=2) error("Phylogeny::likelihood(): Alignment must have only 2 sequences");
	const NY98_TransProbMosaic* phylo = get_phylo();

	mydouble lik(1.0);
	int pos;
	for(pos=0;pos<ct.length();pos++) {
		const int anc = ct[0][pos];
		const int dec = ct[1][pos];
		lik *= mydouble(phylo->get_double(pos,anc,dec));
	}
	return lik;
}

void Phylogeny::add_random_variable(RandomVariable* random_variable) {
	Distribution::add_random_variable(random_variable);
	const int ix = n_random_variables()-1;
	_rv_index.insert(pair< RandomVariable*, int>(random_variable,ix));
}

void Phylogeny::remove_random_variable(RandomVariable* random_variable) {
	error("Phylogeny::remove_random_variable(): not permitted");
}

void Phylogeny::receive_signal_from_parent(const Value* v, const Variable::Signal sgl) {
	Distribution::receive_signal_from_parent(v,sgl);
}
	
} // namespace gcat_omegaMap
