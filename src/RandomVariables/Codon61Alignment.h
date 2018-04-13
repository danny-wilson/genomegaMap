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
 *  Codon61Alignment.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _CODON_61_ALIGNMENT_H_
#define _CODON_61_ALIGNMENT_H_
#include <omegaMap/Variables/Alignment.h>
#include <DAG/RandomVariable.h>
#include <DNA.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class Codon61Alignment : public Alignment, public RandomVariable {
private:
	string _filename;
	int _n;
	int _length;
	vector<string> _label;
	// NB:- rows are sequences, columns are sites
	vector< vector<int> > _seq;
public:
	// Constructor
	Codon61Alignment(string filename, string name="", DAG* dag=0);
	// Copy constructor
	Codon61Alignment(const Codon61Alignment& x);
	// Destructor
	virtual ~Codon61Alignment();
	
	// Implementation of inherited methods
	// Report encoding
	vector<string> encoding() const;
	// Number of sequences
	int n() const;
	// Sequence length
	int length() const;
	// Return a site
	inline int seq(const int i, const int j) const;
	// Return an encoded sequence
	inline const vector<int>& seq(const int i) const;
	// Operator
	inline const vector<int>& operator[](const int i) const;
	// Return all sequences
	const vector< vector<int> >& seqs() const;
	// Return a label
	string label(const int i) const;
	// Return all labels
	const vector<string>& labels() const;
	
protected:
	void tocodon61(DNA& dna, vector< vector<int> >& codonsequence, const int offset=0);
};

int Codon61Alignment::seq(const int i, const int j) const {
	return _seq[i][j];
}

const vector<int>& Codon61Alignment::seq(const int i) const {
	return _seq[i];
}

const vector<int>& Codon61Alignment::operator[](const int i) const {
	return _seq[i];
}	
	
} // namespace gcat_omegaMap

#endif // _CODON_61_ALIGNMENT_H_
