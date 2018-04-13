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
 *  BiallelicCodingGroupCount.h
 *  gcat
 *
 *  Created by Daniel Wilson on 23/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _BIALLELIC_CODING_GROUP_COUNT_VARIABLE_H_
#define _BIALLELIC_CODING_GROUP_COUNT_VARIABLE_H_
#include <omegaMap/Variables/AlleleGroupCount.h>
#include <DAG/RandomVariable.h>
#include <DNA.h>
#include <vector>
#include <myerror.h>

using std::vector;

using namespace gcat;

namespace gcat_omegaMap {
	
class BiallelicCodingGroupCount : public AlleleGroupCount, public RandomVariable {
private:
	// Number of groups
	int _ng;
	// Filename for each group
	vector<string> _filename;
	// Label for each group
	vector<string> _groupname;
	// Number of sequences for each group
	vector<int> _n;
	// Sequence length
	int _length;
	// Abundance of alleles indexed by groups, then sites, then allele (0 or 1)
	vector< vector< vector<int> > > _ct;
	// Indicates the codon is valid biallelic indexed by sites
	vector<bool> _va;
	// Identity of allele 0 indexed by sites
	vector<int> _a0;
	// Identity of allele 1 indexed by sites
	vector<int> _a1;
	
public:
	// Constructor
	BiallelicCodingGroupCount(vector<string> filename, vector<string> groupname, string name="", DAG* dag=0);
	// Copy constructor
	BiallelicCodingGroupCount(const BiallelicCodingGroupCount &x);
	// Destructor
	virtual ~BiallelicCodingGroupCount();
	
	// Implementation of inherited methods
	// Report encoding
	vector<string> encoding() const;
	// Number of groups
	int n_groups() const;
	// Number of sequences
	int n(const int gp) const;
	// Sequence length: assumed same for all
	int length() const;
	// Return the count of alleles for a particular group.
	inline const vector< vector<int> >& operator[](const int gp) const;
	
	// Derived class-specific methods
	// Is the site valid?
	bool is_valid(const int pos) const;
	// Identity of codon 0
	string codon0(const int pos) const;
	// Identity of codon 1
	string codon1(const int pos) const;
	// Identity of amino acid 0
	string aminoacid0(const int pos) const;
	// Identity of amino acid 1
	string aminoacid1(const int pos) const;
	
protected:
	void tally(vector<DNA>& dna, const int offset=0);
};

const vector< vector<int> >& BiallelicCodingGroupCount::operator[](const int gp) const {
	return _ct[gp];
}	
	
} // namespace gcat_omegaMap

#endif // _BIALLELIC_CODING_GROUP_COUNT_VARIABLE_H_


