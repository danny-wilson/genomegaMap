/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the genomegaMap library.
 *
 *  The genomegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The genomegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the genomegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/*
 *  Codon61Count.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *
 */
#include <genomegaMap/RandomVariables/Codon61Count.h>

using namespace gcat;

namespace genomegaMap {
	
const char* Codon61CountEncoding[62] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG","---"};

Codon61Count::Codon61Count(string filename, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61Count"), RandomVariable(), _filename(filename) {
	DNA dna(filename.c_str());
	tally(dna);
	_n = dna.nseq;
	_length = _ct.size();
}

// Alternative constructor
Codon61Count::Codon61Count(vector< vector<int> > _ct_in, string name, DAG* dag) : DAGcomponent(name,dag,"Codon61Count"), RandomVariable(), _filename(""), _ct(_ct_in) {
	// Assumes all rows have the same length (i.e. number of columns)
	if(_ct[0].size()!=62) error("Codon61Count::Codon61Count(): exactly 62 columns expected");
	_length = _ct.size();
	_n = 0;
	for(int i=0;i<62;i++) _n += _ct[0][i];
}

Codon61Count::Codon61Count(const Codon61Count& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _n(x._n), _length(x._length), _ct(x._ct) {
}

Codon61Count::~Codon61Count() {};

vector<string> Codon61Count::encoding() const {
	vector<string> ret(62);
	int i;
	for(i=0;i<ret.size();i++) ret[i] = Codon61CountEncoding[i];
	return ret;
}

int Codon61Count::n() const {
	return _n;
}

int Codon61Count::length() const {
	return _length;
}

void Codon61Count::tally(DNA& dna, const int offset) {
	if(offset<0) error("Codon61Count::tally(): cannot have negative offset");
	if((dna.lseq-offset)%3!=0) error("Codon61Count::tally(): DNA length minus offset isn't a multiple of 3");
	const int tlen = (dna.lseq-offset)/3;
	vector<int> empty_count(62,0);
	_ct = vector< vector<int> >(tlen,empty_count);
	int i,j,ctr;
	for(i=offset,ctr=0;i<dna.lseq;i+=3,ctr+=1) {
		for(j=0;j<dna.nseq;j++)	{
			string triplet = dna.sequence[j].substr(i,3);
			int cdi = dna.tripletToCodon61(triplet);
			// Count partial indels as indels
			if(cdi==-1) {
				stringstream wrnTxt;
				wrnTxt << "Codon61Count::tally(): converted partial indel (" << triplet << ") to indel (---) at position " << i << " in sequence " << j;
				warning(wrnTxt.str().c_str());
				cdi=61;
			}
			// Generate an error for unexpected strings
			if(cdi<0 || cdi>=62) {
				stringstream errTxt;
				errTxt << "Codon61Count::tally(): Unrecognized state (" << triplet << ") detected at position " << i << " in sequence " << j;
				error(errTxt.str().c_str());
			}
			_ct[ctr][cdi]++;
		}
	}
}
	
} // namespace genomegaMap
