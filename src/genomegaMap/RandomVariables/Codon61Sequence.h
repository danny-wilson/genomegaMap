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
 *  Codon61Sequence.h
 *  gcat
 *
 *  Created by Daniel Wilson on 25/04/2010.
 *
 *	This concrete class derived from FactorVectorVariable and RandomVariable
 *	guarantees that the levels are the 61 non-stop codons numbered from 0 to 60.
 *
 */
#ifndef _CODON_61_SEQUENCE_RV_H_
#define _CODON_61_SEQUENCE_RV_H_
#include <Variables/FactorVector.h>
#include <DAG/RandomVariable.h>
#include <DNA.h>

using namespace gcat;

namespace genomegaMap {
	
class Codon61SequenceRV : public FactorVectorVariable, public RandomVariable {
private:
	// Length of vector
	int _n;
	// Values
	vector< int > _value, _previous_value;
	// Record whether the values have changed
	vector< bool > _has_changed;
public:
	// Constructor
	Codon61SequenceRV(const int n, string name="", DAG* dag=0, const vector<int> ivalues=vector<int>(0), const vector<string> svalues=vector<string>(0));
	// Copy constructor
	Codon61SequenceRV(const Codon61SequenceRV& x);
	// Destructor
	virtual ~Codon61SequenceRV();
	
	// Implementation of inherited methods
	// Sequence length
	int length() const;
	// Get value at position i
	int get_int(const int i) const;
	inline int operator[](const int i) const;
	// Number of levels
	int nlevels() const;
	// Levels
	vector<string> levels() const;
	// Get vector of values
	vector<int> get_ints() const;
	// Convert string to integer
	int to_int(const string s) const;
	// Convert integer to string
	string to_string(const int i) const;
	// Has the value changed at position i?
	bool has_changed(const int i) const;
	// Has the value changed at each position?
	vector<bool> has_changed() const;
	
	// Signal (set/propose/accept/revert) a change a single value
	void change_value(const int pos, const int value, Variable::Signal sgl);
	// Signal (set/propose/accept/revert) a change to one or more values
	void change_value(vector<int> &value, Variable::Signal sgl);
};

int Codon61SequenceRV::operator[](const int site) const {
	return _value[site];
}	
	
} // namespace genomegaMap

#endif // _CODON_61_SEQUENCE_RV_H_

