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
 *  Codon61GroupCount.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _CODON_61_GROUP_COUNT_VARIABLE_H_
#define _CODON_61_GROUP_COUNT_VARIABLE_H_
#include <omegaMap/Variables/AlleleGroupCount.h>
#include <DAG/RandomVariable.h>
#include <DNA.h>
#include <vector>
#include <myerror.h>

using std::vector;

using namespace gcat;

namespace gcat_omegaMap {
	
class Codon61GroupCount : public AlleleGroupCount, public RandomVariable {
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
	// Abundance of each codon indexed by groups, then sites, then sequences
	vector< vector< vector<int> > > _ct;
	
public:
	// Constructor
	Codon61GroupCount(vector<string> filename, vector<string> groupname, string name="", DAG* dag=0);
	// Copy constructor
	Codon61GroupCount(const Codon61GroupCount &x);
	// Destructor
	virtual ~Codon61GroupCount();
	
	// Implementation of inherited methods
	// Report encoding
	vector<string> encoding() const;
	// Number of groups
	int n_groups() const;
	// Number of sequences
	int n(const int gp) const;
	// Sequence length: assumed same for all
	int length() const;
	// Return the counts for a particular group
	inline const vector< vector<int> >& operator[](const int gp) const;
	
protected:
	void tally(DNA& dna, const int gp, const int offset=0);
};

const vector< vector<int> >& Codon61GroupCount::operator[](const int gp) const {
	return _ct[gp];
}	
	
} // namespace gcat_omegaMap


#endif // _CODON_61_GROUP_COUNT_VARIABLE_H_


