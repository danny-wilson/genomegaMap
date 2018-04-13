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
 *  Codon61Alignment.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/Codon61Alignment.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* Codon61AlignmentEncoding[62] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG","___"};

Codon61Alignment::Codon61Alignment(string filename, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61Alignment"), RandomVariable(), _filename(filename) {
	DNA dna(filename.c_str());
	tocodon61(dna,_seq);
	_n = _seq.size();
	_length = _seq[0].size();
}

Codon61Alignment::Codon61Alignment(const Codon61Alignment& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _n(x._n), _length(x._length), _label(x._label), _seq(x._seq) {
}

Codon61Alignment::~Codon61Alignment() {};

vector<string> Codon61Alignment::encoding() const {
	vector<string> ret(62);
	int i;
	for(i=0;i<ret.size();i++) ret[i] = Codon61AlignmentEncoding[i];
	return ret;
}

int Codon61Alignment::n() const {
	return _n;
}

int Codon61Alignment::length() const {
	return _length;
}

const vector< vector<int> >& Codon61Alignment::seqs() const {
	return _seq;
}

string Codon61Alignment::label(const int i) const {
	return _label[i];
}

const vector<string>& Codon61Alignment::labels() const {
	return _label;
}

void Codon61Alignment::tocodon61(DNA& dna, vector< vector<int> >& codonsequence, const int offset) {
	if(offset<0) error("Codon61Alignment::tocodon61(): cannot have negative offset");
	if((dna.lseq-offset)%3!=0) error("Codon61Alignment::tocodon61(): DNA length minus offset isn't a multiple of 3");
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
