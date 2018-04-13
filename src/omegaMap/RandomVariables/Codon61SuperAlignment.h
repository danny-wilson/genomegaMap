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
 *  Codon61SuperAlignment.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _CODON_61_SUPER_ALIGNMENT_H_
#define _CODON_61_SUPER_ALIGNMENT_H_
#include <omegaMap/Variables/SuperAlignment.h>
#include <DAG/RandomVariable.h>
#include <DNA.h>

using namespace gcat;

namespace gcat_omegaMap {
	
class Codon61SuperAlignment : public SuperAlignment, public RandomVariable {
private:
	vector<string> _filename;
	int _ng;
	vector<int> _n;
	int _length;
	vector< vector<string> > _label;
	vector<string> _groupname;
	// NB:- indexed by groups, then sequences, then sites
	vector< vector< vector<int> > > _seq;
public:
	// Constructor
	Codon61SuperAlignment(vector<string> filename, vector<string> groupname, string name="", DAG* dag=0);
	// Copy constructor
	Codon61SuperAlignment(const Codon61SuperAlignment& x);
	// Destructor
	virtual ~Codon61SuperAlignment();
	
	// Implementation of inherited methods
	// Report encoding
	vector<string> encoding() const;
	// Number of sequences
	int n(const int gp) const;
	// Number of groups
	int n_groups() const;
	// Sequence length
	int length() const;
	// Return a site
	inline int seq(const int gp, const int id, const int pos) const;
	// Return an encoded sequence
	inline const vector<int>& seq(const int gp, const int id) const;
	// Operator
	inline const vector< vector<int> >& operator[](const int gp) const;
	// Return all sequences in a group
	inline const vector< vector<int> >& seqs(const int gp) const;
	// Return all sequences
	inline const vector< vector< vector<int> > >& seqs() const;
	// Return a group name
	const string group_name(const int gp) const;
	// Return all group names
	const vector<string>& group_names() const;
	// Return a label
	string label(const int gp, const int i) const;
	// Return all labels for a group
	const vector<string>& labels(const int gp) const;
	// Return all labels
	const vector< vector<string> >& labels() const;
	
protected:
	void tocodon61(DNA& dna, vector< vector<int> >& codonsequence, const int offset=0);
};

int Codon61SuperAlignment::seq(const int gp, const int id, const int pos) const {
	return _seq[gp][id][pos];
}

const vector<int>& Codon61SuperAlignment::seq(const int gp, const int id) const {
	return _seq[gp][id];
}

const vector< vector<int> >& Codon61SuperAlignment::seqs(const int gp) const {
	return _seq[gp];
}

const vector< vector< vector<int> > >& Codon61SuperAlignment::seqs() const {
	return _seq;
}

const vector< vector<int> >& Codon61SuperAlignment::operator[](const int i) const {
	return _seq[i];
}	
	
} // namespace gcat_omegaMap

#endif // _CODON_61_SUPER_ALIGNMENT_H_
