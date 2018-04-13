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
 *  BiallelicCodingGroupCount.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 23/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/BiallelicCodingGroupCount.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* SixtyoneCodons[61] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};
const char* SixtyoneAminoAcids[61] = {"Phe","Phe","Leu","Leu","Ser","Ser","Ser","Ser","Tyr","Tyr","Cys","Cys","Trp","Leu","Leu","Leu","Leu","Pro","Pro","Pro","Pro","His","His","Gln","Gln","Arg","Arg","Arg","Arg","Ile","Ile","Ile","Met","Thr","Thr","Thr","Thr","Asn","Asn","Lys","Lys","Ser","Ser","Arg","Arg","Val","Val","Val","Val","Ala","Ala","Ala","Ala","Asp","Asp","Glu","Glu","Gly","Gly","Gly","Gly"};

BiallelicCodingGroupCount::BiallelicCodingGroupCount(vector<string> filename, vector<string> groupname, string name, DAG* dag) : DAGcomponent(name,dag,"BiallelicCodingGroupCount"), RandomVariable(), _filename(filename), _groupname(groupname), _ng(filename.size()) {
	if(groupname.size()!=_ng) error("BiallelicCodingGroupCount: filename and groupname must be vectors of the same length");
	_n.resize(_ng);
	int i;
	vector<DNA> dna(_ng);
	for(i=0;i<_ng;i++) {
		dna[i].readFASTA(filename[i].c_str());
		_n[i] = dna[i].nseq;
		if(i==0) _length = dna[i].lseq;
		else if(dna[i].lseq!=_length) error("BiallelicCodingGroupCount: all groups must have same sequence lengths");
	}
	tally(dna);
	// Adjust to number of codons
	_length = _ct[0].size();
	
	// Summarize sites
	int syn=0, nsyn=0, nvar=0;
	for(i=0;i<_length;i++) {
		if(codon0(i)==codon1(i)) ++nvar;
		else if(aminoacid0(i)==aminoacid1(i)) ++syn;
		else ++nsyn;
	}
	cout << "Found " << syn << " synonymous, " << nsyn << " non-synonymous and " << nvar << " invariant sites\n";
}

BiallelicCodingGroupCount::BiallelicCodingGroupCount(const BiallelicCodingGroupCount& x) : DAGcomponent(x), RandomVariable(x), _filename(x._filename), _groupname(x._groupname), _ng(x._ng), _n(x._n), _length(x._length), _ct(x._ct) {
}

BiallelicCodingGroupCount::~BiallelicCodingGroupCount() {};

vector<string> BiallelicCodingGroupCount::encoding() const {
	vector<string> ret(2,"0");
	ret[1] = "1";
	return ret;
}

int BiallelicCodingGroupCount::n_groups() const {
	return _ng;
}

int BiallelicCodingGroupCount::n(const int gp) const {
	return _n[gp];
}

int BiallelicCodingGroupCount::length() const {
	return _length;
}

void BiallelicCodingGroupCount::tally(vector<DNA>& dna, const int offset) {
	if(offset<0) error("BiallelicCodingGroupCount::tally(): cannot have negative offset");
	// Assumes we've already checked groups have the same length (this happens in the constructor)
	if((dna[0].lseq-offset)%3!=0) error("BiallelicCodingGroupCount::tally(): DNA length minus offset isn't a multiple of 3");
	const int tlen = (dna[0].lseq-offset)/3;
	// Allocate memory and initialize arrays
	int gp;
	vector<int> empty_count(2,0);
	_ct.resize(_ng);
	for(gp=0;gp<_ng;gp++) _ct[gp] = vector< vector<int> >(tlen,empty_count);
	_va = vector<bool>(tlen,true);
	_a0 = vector<int>(tlen,0);
	_a1 = vector<int>(tlen,0);
	// Cycle through sites
	cout << "The following sites were invalid:";
	int nind=0, npol=0;
	int i,j,ctr;
	for(i=offset,ctr=0;i<dna[0].lseq;i+=3,ctr+=1) {
		string triplet = dna[0].sequence[0].substr(i,3);
		const int cd0 = dna[0].tripletToCodon61(triplet);
		_ct[0][ctr][0] = 1;
		int cd1;
		// Number of alleles
		int nalleles = 1;
		bool anyindels = false;
		// Disallow indels
		if(cd0==-1) {
			anyindels = true;
		}
		// Generate an error for unexpected strings
		else if(cd0<0 || cd0>=62) {
			stringstream errTxt;
			errTxt << "BiallelicCodingGroupCount::tally(): Unrecognized state (" << triplet << ") detected in group " << 0 << " at position " << i << " in sequence " << 0;
			error(errTxt.str().c_str());
		}
		// Continue with the remaining sequences
		else {
			for(gp=0;gp<_ng;gp++) {
				const int jmin = (gp==0) ? 1 : 0;
				for(j=jmin;j<_n[gp];j++)	{
					triplet = dna[gp].sequence[j].substr(i,3);
					int cdi = dna[gp].tripletToCodon61(triplet);
					// Disallow indels
					if(cdi==-1) {
						anyindels = true;
						break;
					}
					// Generate an error for unexpected strings
					if(cdi<0 || cdi>=62) {
						stringstream errTxt;
						errTxt << "BiallelicCodingGroupCount::tally(): Unrecognized state (" << triplet << ") detected in group " << gp << " at position " << i << " in sequence " << j;
						error(errTxt.str().c_str());
					}
					// Which allele is it?
					if(cdi==cd0) {
						// Allele 0, so augment the count
						_ct[gp][ctr][0]++;
					}
					else if(nalleles==1) {
						// The second allele, so record its identity
						cd1 = cdi;
						_ct[gp][ctr][1] = 1;
						++nalleles;
					}
					else if(cdi==cd1) {
						// Allele 1, so augment the count
						_ct[gp][ctr][1]++;
					}
					else {
						// A third allele, which invalidates the site
						++nalleles;
						break;
					}
				}
				if(anyindels) break;
				if(nalleles>2) break;
			}
		}
		if(anyindels) {
			for(gp=0;gp<_ng;gp++) {
				_ct[gp][ctr][0] = _n[gp];
				_ct[gp][ctr][1] = 0;
			}
			_va[ctr] = false;
			_a0[ctr] = -1;
			_a1[ctr] = -1;
			cout << " " << ctr << "I" ;
			++nind;
		}
		else if(nalleles>2) {
			for(gp=0;gp<_ng;gp++) {
				_ct[gp][ctr][0] = _n[gp];
				_ct[gp][ctr][1] = 0;
			}
			_va[ctr] = false;
			_a0[ctr] = -2;
			_a1[ctr] = -2;
			cout << " " << ctr << "P" ;
			++npol;
		}
		else if(_va[ctr]) {
			_a0[ctr] = cd0;
			_a1[ctr] = (nalleles==1) ? cd0 : cd1;
		}
	}
	cout << endl;
	cout << "Totals: " << tlen-nind-npol << " valid, " << nind << " with indels, " << npol << " too polymorphic" << endl;
}

bool BiallelicCodingGroupCount::is_valid(const int pos) const {
	return _va[pos];
}

string BiallelicCodingGroupCount::codon0(const int pos) const {
	return (_va[pos]) ? SixtyoneCodons[_a0[pos]] : "???";
}

string BiallelicCodingGroupCount::codon1(const int pos) const {
	return (_va[pos]) ? SixtyoneCodons[_a1[pos]] : "???";
}

string BiallelicCodingGroupCount::aminoacid0(const int pos) const {
	return (_va[pos]) ? SixtyoneAminoAcids[_a0[pos]] : "???";
}

string BiallelicCodingGroupCount::aminoacid1(const int pos) const {
	return (_va[pos]) ? SixtyoneAminoAcids[_a1[pos]] : "???";
}
	
} // namespace gcat_omegaMap

