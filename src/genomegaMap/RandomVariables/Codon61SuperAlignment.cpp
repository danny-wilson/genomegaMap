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
 *  Codon61SuperAlignment.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/Codon61SuperAlignment.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* Codon61SuperAlignmentEncoding[62] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG","___"};

Codon61SuperAlignment::Codon61SuperAlignment(vector<string> filename, vector<string> groupname, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61SuperAlignment"), RandomVariable(), _filename(filename), _groupname(groupname), _ng(filename.size()) {
	if(_groupname.size()!=_filename.size()) error("Codon61SuperAlignment: filename and groupname vectors must have same length");
	_n.resize(_ng);
	_seq.resize(_ng);
	_label.resize(_ng);
	int g;
	for(g=0;g<_ng;g++) {
		DNA dna(filename[g].c_str());
		tocodon61(dna,_seq[g]);
		_n[g] = _seq[g].size();
		if(g==0) _length = _seq[g][0].size();
		else if(_seq[g][0].size()!=_length) error("Codon61SuperAlignment: all groups must identical sequence lengths");
	}
}

Codon61SuperAlignment::Codon61SuperAlignment(const Codon61SuperAlignment& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _groupname(x._groupname), _ng(x._ng), _n(x._n), _length(x._length), _label(x._label), _seq(x._seq) {
}

Codon61SuperAlignment::~Codon61SuperAlignment() {};

vector<string> Codon61SuperAlignment::encoding() const {
	vector<string> ret(62);
	int i;
	for(i=0;i<ret.size();i++) ret[i] = Codon61SuperAlignmentEncoding[i];
	return ret;
}

int Codon61SuperAlignment::n(const int gp) const {
	return _n[gp];
}

int Codon61SuperAlignment::n_groups() const {
	return _ng;
}

int Codon61SuperAlignment::length() const {
	return _length;
}

const string Codon61SuperAlignment::group_name(const int gp) const {
	return _groupname[gp];
}

const vector<string>& Codon61SuperAlignment::group_names() const {
	return _groupname;
}

string Codon61SuperAlignment::label(const int gp, const int id) const {
	return _label[gp][id];
}

const vector<string>& Codon61SuperAlignment::labels(const int gp) const {
	return _label[gp];
}

const vector< vector<string> >& Codon61SuperAlignment::labels() const {
	return _label;
}

void Codon61SuperAlignment::tocodon61(DNA& dna, vector< vector<int> >& codonsequence, const int offset) {
	if(offset<0) error("Codon61SuperAlignment::tocodon61(): cannot have negative offset");
	if((dna.lseq-offset)%3!=0) error("Codon61SuperAlignment::tocodon61(): DNA length minus offset isn't a multiple of 3");
	const int tlen = (dna.lseq-offset)/3;
	vector<int> blank_sequence(tlen);
	codonsequence = vector< vector<int> >(dna.nseq,blank_sequence);
	int i,j,ctr;
	for(i=offset,ctr=0;i<dna.lseq;i+=3,ctr+=1) {
		for(j=0;j<dna.nseq;j++)	{
			string triplet = dna.sequence[j].substr(i,3);
			codonsequence[j][ctr] = dna.tripletToCodon61(triplet);
		}
	}
}
	
} // namespace gcat_omegaMap
