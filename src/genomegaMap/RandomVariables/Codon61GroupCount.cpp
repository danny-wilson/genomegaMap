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
 *  Codon61GroupCount.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/Codon61GroupCount.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* Codon61GroupCountEncoding[62] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG","---"};

Codon61GroupCount::Codon61GroupCount(vector<string> filename, vector<string> groupname, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61GroupCount"), RandomVariable(), _filename(filename), _groupname(groupname), _ng(filename.size()) {
	if(groupname.size()!=_ng) error("Codon61GroupCount: filename and groupname must be vectors of the same length");
	_n.resize(_ng);
	_ct.resize(_ng);
	int i;
	for(i=0;i<_ng;i++) {
		DNA dna(filename[i].c_str());
		tally(dna,i);
		_n[i] = dna.nseq;
		if(i==0) _length = _ct[0].size();
		else if(_ct[i].size()!=_length) error("Codon61GroupCount: all groups must have same sequence lengths");
	}
}

Codon61GroupCount::Codon61GroupCount(const Codon61GroupCount& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _groupname(x._groupname), _ng(x._ng), _n(x._n), _length(x._length), _ct(x._ct) {
}

Codon61GroupCount::~Codon61GroupCount() {};

vector<string> Codon61GroupCount::encoding() const {
	vector<string> ret(62);
	int i;
	for(i=0;i<ret.size();i++) ret[i] = Codon61GroupCountEncoding[i];
	return ret;
}

int Codon61GroupCount::n_groups() const {
	return _ng;
}

int Codon61GroupCount::n(const int gp) const {
	return _n[gp];
}

int Codon61GroupCount::length() const {
	return _length;
}

void Codon61GroupCount::tally(DNA& dna, const int gp, const int offset) {
	if(offset<0) error("Codon61GroupCount::tally(): cannot have negative offset");
	if((dna.lseq-offset)%3!=0) error("Codon61GroupCount::tally(): DNA length minus offset isn't a multiple of 3");
	const int tlen = (dna.lseq-offset)/3;
	vector<int> empty_count(62,0);
	_ct[gp] = vector< vector<int> >(tlen,empty_count);
	int i,j,ctr;
	for(i=offset,ctr=0;i<dna.lseq;i+=3,ctr+=1) {
		for(j=0;j<dna.nseq;j++)	{
			string triplet = dna.sequence[j].substr(i,3);
			int cdi = dna.tripletToCodon61(triplet);
			// Count partial indels as indels
			if(cdi==-1) {
				stringstream wrnTxt;
				wrnTxt << "Codon61GroupCount::tally(): converted partial indel (" << triplet << ") to indel (---) at position " << i << " in sequence " << j;
				warning(wrnTxt.str().c_str());
				cdi=61;
			}
			// Generate an error for unexpected strings
			if(cdi<0 || cdi>=62) {
				stringstream errTxt;
				errTxt << "Codon61GroupCount::tally(): Unrecognized state (" << triplet << ") detected at position " << i << " in sequence " << j;
				error(errTxt.str().c_str());
			}
			_ct[gp][ctr][cdi]++;
		}
	}
}
	
} // namespace gcat_omegaMap
