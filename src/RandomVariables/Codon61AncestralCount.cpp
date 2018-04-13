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
 *  Codon61Count.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/Codon61AncestralCount.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* Codon61AncestralCountEncoding[62] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG","___"};

Codon61AncestralCount::Codon61AncestralCount(string filename, string ancestor, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61Count"), RandomVariable(), _filename(filename) {
	DNA dna(filename.c_str());
	int i;
	int anc = -1;
	for(i=0;i<dna.nseq;i++) {
		if(dna.label[i]==ancestor) {
			if(anc!=-1) error("Codon61AncestralCount: multiple sequence labels match ancestor");
			anc = i;
		}
	}
	if(anc==-1) error("Codon61AncestralCount: no sequence labels match ancestor");
	tally(dna,anc);
	_n = dna.nseq-1;
	_length = _ct.size();
}

Codon61AncestralCount::Codon61AncestralCount(const Codon61AncestralCount& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _n(x._n), _length(x._length), _ct(x._ct), _anc(x._anc) {
}

Codon61AncestralCount::~Codon61AncestralCount() {};

vector<string> Codon61AncestralCount::encoding() const {
	vector<string> ret(62);
	int i;
	for(i=0;i<ret.size();i++) ret[i] = Codon61AncestralCountEncoding[i];
	return ret;
}

int Codon61AncestralCount::n() const {
	return _n;
}

int Codon61AncestralCount::length() const {
	return _length;
}

int Codon61AncestralCount::ancestor(const int site) const {
	return _anc[site];
}

void Codon61AncestralCount::tally(DNA& dna, const int ancestor, const int offset) {
	if(offset<0) error("Codon61AncestralCount::tally(): cannot have negative offset");
	if((dna.lseq-offset)%3!=0) error("Codon61AncestralCount::tally(): DNA length minus offset isn't a multiple of 3");
	const int tlen = (dna.lseq-offset)/3;
	vector<int> empty_count(61,0);
	_ct = vector< vector<int> >(tlen,empty_count);
	_anc = vector<int>(tlen,-1);
	int i,j,ctr;
	for(i=offset,ctr=0;i<dna.lseq;i+=3,ctr+=1) {
		for(j=0;j<dna.nseq;j++)	{
			string triplet = dna.sequence[j].substr(i,3);
			const int cdi = dna.tripletToCodon61(triplet);
			if(cdi<0 || cdi>=61) {
				error("Codon61AncestralCount::tally(): Unrecognized state or indel detected");
			}
			if(j==ancestor) {
				_anc[ctr] = cdi;
			}
			else {
				_ct[ctr][cdi]++;
			}
		}
	}
}
	
} // namespace gcat_omegaMap
