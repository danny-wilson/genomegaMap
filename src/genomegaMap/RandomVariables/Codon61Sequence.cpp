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
 *  Codon61Sequence.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 25/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <omegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

namespace gcat_omegaMap {
	
const char* Codon61SequenceRVStrings[61] = {"TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TGT","TGC","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"};

Codon61SequenceRV::Codon61SequenceRV(const int n, string name, DAG* dag, const vector<int> ivalues, const vector<string> svalues) : _n(n), DAGcomponent(name,dag,"Codon61SequenceRV"), RandomVariable(), _value(vector<int>(n)), _previous_value(vector<int>(n)), _has_changed(n,true) {
	if(ivalues.size()==0) {
		if(svalues.size()==0) {
			// If zero is an invalid value, this should cause an error
			string s = to_string(0);
			_value = vector<int>(_n,0);
		}
		else if(svalues.size()==1) {
			// If an invalid value, this should cause an error
			_value = vector<int>(_n,to_int(svalues[0]));
		}
		else if(svalues.size()==_n) {
			int i;
			_value = vector<int>(_n);
			for(i=0;i<_n;i++) {
				_value[i] = to_int(svalues[i]);
			}
		}
		else error("Codon61SequenceRV: initial string vector incompatible size");
	}
	else {
		if(svalues.size()!=0) error("Codon61SequenceRV: cannot specify initial integer and string vector");
		if(ivalues.size()==1) {
			// If an invalid value, this should cause an error
			const int i = ivalues[0];
			string s = to_string(i);
			_value = vector<int>(_n,i);
		}
		else if(ivalues.size()==_n) {
			int i;
			for(i=0;i<_n;i++) {
				// If an invalid value, this should cause an error
				string s = to_string(ivalues[i]);
			}
			_value = ivalues;
		}
		else error("Codon61SequenceRV: initial integer vector incompatible size");
	}
}

Codon61SequenceRV::Codon61SequenceRV(const Codon61SequenceRV& x) : _n(x._n), DAGcomponent(x), RandomVariable(x), _value(x._value), _previous_value(x._previous_value), _has_changed(x._has_changed) {
}

Codon61SequenceRV::~Codon61SequenceRV() {};

int Codon61SequenceRV::length() const {
	return _n;
}

int Codon61SequenceRV::get_int(const int i) const {
	return _value[i];
}

int Codon61SequenceRV::nlevels() const {
	return 61;
}

vector<string> Codon61SequenceRV::levels() const {
	vector<string> s(61);
	int i;
	for(i=0;i<61;i++) s[i] = Codon61SequenceRVStrings[i];
	return s;
}

vector<int> Codon61SequenceRV::get_ints() const {
	return _value;
}

int Codon61SequenceRV::to_int(const string s) const {
	// Not case sensitive
	string S = s;
	int i;
	for(i=0;i<S.length();i++) S[i] = toupper(S[i]);
	for(i=0;i<61;i++) {
		if(S==string(Codon61SequenceRVStrings[i])) break;
	}
	if(i==61) {
		string errMsg = "Codon61SequenceRV::to_int(): ";
		errMsg += s + " not a valid level";
		error(errMsg.c_str());
	}
	return i;
}

string Codon61SequenceRV::to_string(const int i) const {
	if(i<0 || i>=61) {
		stringstream errMsg;
		errMsg << "Codon61SequenceRV::to_string(): " << i << " not a valid level";
		error(errMsg.str().c_str());
	}
	return Codon61SequenceRVStrings[i];
}

bool Codon61SequenceRV::has_changed(const int i) const {
	return _has_changed[i];
}

vector<bool> Codon61SequenceRV::has_changed() const {
	return _has_changed;
}

void Codon61SequenceRV::change_value(const int pos, const int value, Variable::Signal sgl) {
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_previous_value = _value;
		if(value<0 || value>=61) error("ContinuousMosaicRV::change_value(): value out of range");
		_has_changed = vector<bool>(_n,false);
		if(value!=_value[pos]) _has_changed[pos] = true;
		_value[pos] = value;
	}
	else if(sgl==Variable::_ACCEPT) {
		_has_changed[pos] = false;
	}
	else if(sgl==Variable::_REVERT) {
		_value = _previous_value;
		_has_changed[pos] = false;
	}
	else {
		error("ContinuousMosaicRV::change_value(): unexpected Variable signal");
	}
	act_on_signal(sgl);
	send_signal_to_children(sgl);
}

void Codon61SequenceRV::change_value(vector<int> &value, Variable::Signal sgl) {
	if(value.size()!=length()) error("ContinuousMosaicRV::change_value(): value vector incompatible size");
	int i;
	if(sgl==Variable::_SET || sgl==Variable::_PROPOSE) {
		_previous_value = _value;
		_has_changed = vector<bool>(_n,false);
		for(i=0;i<length();i++) {
			if(_value[i]!=value[i]) {
				_value[i] = value[i];
				_has_changed[i] = true;
			}
		}
	}
	else if(sgl==Variable::_ACCEPT) {
		_has_changed = vector<bool>(_n,false);
	}
	else if(sgl==Variable::_REVERT) {
		_value = _previous_value;
		_has_changed = vector<bool>(_n,false);
	}
	else {
		error("ContinuousMosaicRV::change_value(): unexpected Variable signal");
	}
	act_on_signal(sgl);
	send_signal_to_children(sgl);
}
	
} // namespace gcat_omegaMap
