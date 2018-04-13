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
 *  SuperAlignment.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/17/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 *	A SuperAlignment is a collection of aligned alignments, which may represent
 *	species, populations or any other arbitrary grouping 
 *
 */
#ifndef _SUPER_ALIGNMENT_VARIABLE_H_
#define _SUPER_ALIGNMENT_VARIABLE_H_
#include <DAG/Value.h>
#include <vector>
#include <myerror.h>
#include <Properties/Length.h>

using std::vector;

using namespace gcat;

namespace gcat_omegaMap {
	
class SuperAlignment : public Value, public LengthProperty {
public:
	// Constructor
	SuperAlignment() {};
	// Copy constructor
	SuperAlignment(const SuperAlignment &x) {};
	// Destructor
	virtual ~SuperAlignment() {};
	
	// Report encoding
	virtual vector<string> encoding() const = 0;
	// Number of groups
	virtual int n_groups() const = 0;
	// Number of sequences in a group
	virtual int n(const int gp) const = 0;
	// Sequence length (assumed same in all groups). Inherited from LengthProperty
	//virtual int length() const = 0;
	// Return a site
	virtual int seq(const int gp, const int id, const int pos) const = 0;
	// Return an encoded sequence
	virtual const vector<int>& seq(const int gp, const int id) const = 0;
	// Operator
	virtual const vector< vector<int> >& operator[](const int gp) const = 0;
	// Return all sequences in a group
	virtual const vector< vector<int> >& seqs(const int gp) const = 0;
	// Return all sequences
	virtual const vector< vector< vector<int> > >& seqs() const = 0;
	// Return a group name
	virtual const string group_name(const int gp) const = 0;
	// Return all group names
	virtual const vector<string>& group_names() const = 0;
	// Return a label
	virtual string label(const int gp, const int id) const = 0;
	// Return all labels for a group
	virtual const vector<string>& labels(const int gp) const = 0;
	// Return all labels
	virtual const vector<string>& label() const = 0;
	
	// Print header (implementation of inherited method)
	virtual void print_header(ostream& out, string sep) {
		myutils::warning("SuperAlignment::print_header(): no print method available");
	}
	// Print value (implementation of inherited method)	
	virtual void print(ostream& out, string sep) {
	}
};
	
} // namespace gcat_omegaMap

#endif // _SUPER_ALIGNMENT_VARIABLE_H_
