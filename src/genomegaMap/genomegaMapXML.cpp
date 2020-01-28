/*  Copyright 2018 Daniel Wilson.
 *
 *  genomegaMapXML.cpp
 *  Part of the omegaMap library.
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
#include <gsl/gsl_errno.h>
#include <RandomVariables/Continuous.h>
#include <RandomVariables/ContinuousVector.h>
#include <Distributions/DistributionsXML.h>
#include <RandomVariables/RandomVariablesXML.h>
#include <Transformations/TransformationsXML.h>
#include <Inference/InferenceXML.h>
#include <genomegaMap/omegaMapXML.h>
#include <genomegaMap/Distributions/genomegaMap.h>
#include <genomegaMap/RandomVariables/Codon61Count.h>
#include <genomegaMap/RandomVariables/Codon61Sequence.h>
#include <genomegaMap/Transformations/NY98_PDRM.h>
#include <genomegaMap/Transformations/Omega2GammaVector.h>
#include <genomegaMap/Transformations/ParentDependentRateMatrix.h>
#include <genomegaMap/Inference/MCMC/genomegaMapMoves.h>
#include <genomegaMap/omegaMap1.0.xsd.h>
#include <stdexcept>

using namespace gcat;

namespace genomegaMap {

// DISTRIBUTIONS

genomegaMapUnlinked_XMLParser::genomegaMapUnlinked_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapUnlinked_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 2;
	const char* attrNames[nattr] = {"id","mut"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	new omegaMapUnlinked(sattr[0],getDAG());
	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
}

// RANDOM VARIABLES

codon_count_XMLParser::codon_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_count_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 5;
	const char* attrNames[nattr] = {"id","distribution","file","format","encoding"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	if(sattr[3]=="counts") {
		if(sattr[4]=="codon61") {
			// Based on https://github.com/danny-wilson/gcat-lmm/blob/master/src/lmm/lmmXML.cpp continuous_matrix_file_XMLParser::continuous_matrix_file_XMLParser
			// Read values from file
			time_t start = clock();
			vector< vector<int> > mat(0);
			vector<int> val(0);
			ifstream in(sattr[2].c_str());
			if(!in.good()) {
				string errMsg = "codon_count_XMLParser: file " + sattr[2] + " not found";
				error(errMsg.c_str());
			}
			// Read the file one line at a time: impose constraint that the number of elements is the same for all lines
			int intin;
			string line;
			getline(in,line);
			bool firstline = true;
			int nrows=0, ncols;
			while(!in.eof()) {
				if(!in.good()) {
					string errMsg = "codon_count_XMLParser: unexpected problem reading " + sattr[2];
					error(errMsg.c_str());
				}
				// Convert to vector of doubles
				istringstream linein(line);
				int n = 0;
				while(linein >> intin) {
					if(linein.bad()) {
						string errMsg = "codon_count_XMLParser: could not read integer in " + sattr[2];
						error(errMsg.c_str());
					}
					val.push_back(intin);
					++n;
				}
				if(firstline) {
					ncols = n;
					firstline = false;
					if(ncols!=62) {
						stringstream errMsg;
						errMsg << "codon_count_XMLParser: row " << nrows+1 << " contained " << n << " columns, not " << 62 << " as expected";
						error(errMsg.str().c_str());
					}
					mat.push_back(val);
					val = vector<int>(0);
				} else {
					if(ncols!=n) {
						stringstream errMsg;
						errMsg << "codon_count_XMLParser: row " << nrows+1 << " contained " << n << " columns, not " << ncols << " as expected";
						error(errMsg.str().c_str());
					}
					mat.push_back(val);
					val = vector<int>(0);
				}
				++nrows;
				getline(in,line);
			}
			cout << "Read " << nrows << " rows " << ncols << " columns  from " << sattr[2] << " in " << (clock()-start)/CLOCKS_PER_SEC << " s" << endl;
			// Instantiate the variable
			new Codon61Count(mat,sattr[0],getDAG());
		}
		else error("codon_count_XMLParser only codon61 encoding supported");
	} else if(sattr[3]=="fasta") {
		// Instantiate the variable
		if(sattr[4]=="codon61") {
			new Codon61Count(sattr[2],sattr[0],getDAG());
		}
		else error("codon_count_XMLParser only codon61 encoding supported");
	} else {
		error("codon_count_XMLParser: only fasta and counts formats supported");
	}
	if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],attrNames[1],sattr[1]);
	getDAG()->set_constant(sattr[0]);
}

codon_sequence_XMLParser::codon_sequence_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_sequence_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 3;
	const char* attrNames[nattr] = {"id","distribution","encoding"};
	sattr = attributesToStrings(nattr,attrNames,attrs);
	if(sattr[2]!="codon61") error("codon_sequence_XMLParser: currently encoding must equal \"codon61\"");
	// Don't instantiate the variable until the values have been read in
}

void codon_sequence_XMLParser::implement_characters(const XMLCh* const chars, const XMLSize_t length) {
	string message = "";
	int i;
	for(i=0;i<length;i++) message += chars[i];
	istringstream oss(message);
	vector<string> val(0);
	string val_j;
	while(!(oss >> val_j).fail()) {
		val.push_back(val_j);
	}
	if(val.size()==0) error("codon_sequence_XMLParser: no values entered");
	// Are they entered as strings or ints?
	vector<int> ival(val.size());
	if(from_string<int>(ival[0],val[0])) {
		for(i=1;i<val.size();i++) {
			if(!from_string<int>(ival[i],val[i])) {
				error("codon_sequence_XMLParser: mixed integers and strings entered");
			}
		}
		new Codon61SequenceRV(val.size(),sattr[0],getDAG(),ival);
	}
	else {
		new Codon61SequenceRV(val.size(),sattr[0],getDAG(),vector<int>(0),val);
	}
	if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],"distribution",sattr[1]);
}

// TRANSFORMATIONS

ny98_pdrm_XMLParser::ny98_pdrm_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<ny98_pdrm_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 6;
	const char* attrNames[nattr] = {"id","theta","kappa","omega","pi","length"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	int int_length;
	if(!from_string<int>(int_length,sattr[5])) {
		RandomVariable* rv = getDAG()->get_random_variable(sattr[5]);
		if(rv==0) error("ny98_pdrm_XMLParser: could not convert length to int nor find named variable");
		LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
		if(lp==0) error("ny98_pdrm_XMLParser: named variable does not have length property");
		int_length = lp->length();
	}
	if(int_length<=0) error("ny98_pdrm_XMLParser: length must be a positive integer");
	new NY98_ParentDependentRateMatrix(int_length,sattr[0],getDAG());
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[4],sattr[4]);
}

omega_to_gamma_vector_XMLParser::omega_to_gamma_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omega_to_gamma_vector_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 2;
	const char* attrNames[nattr] = {"id","omega"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	// Automatically obtain length of the quantile mosaic
	int int_length;
	RandomVariable* rv = getDAG()->get_random_variable(sattr[1]);
	if(rv==0) error("omega_to_gamma_vector_XMLParser: could not find named variable for omega");
	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
	if(lp==0) error("omega_to_gamma_vector_XMLParser: named variable quantile does not have length property");
	int_length = lp->length();
	// Instantiate
	new Omega2GammaVectorTransform(int_length,sattr[0],getDAG());
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
}

parent_dependent_rate_matrix_XMLParser::parent_dependent_rate_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<parent_dependent_rate_matrix_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 4;
	const char* attrNames[nattr] = {"id","theta","kappa","pi"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	new ParentDependentRateMatrix(sattr[0],getDAG());
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
}

// INFERENCE

codon61_sequence_gibbs_sampler_XMLParser::codon61_sequence_gibbs_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon61_sequence_gibbs_sampler_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 2;
	const char* attrNames[nattr] = {"parameter","weight"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	vector<string> target(1,sattr[0]);
	double double_weight;
	if(!from_string<double>(double_weight,sattr[1])) error("codon61_sequence_gibbs_sampler_XMLParser: cannot convert parameter weight to double");
	// Get _mcmc from parent parser via dynamic type-checking
	MCMC_XMLParser* MCMC_XMLParser_parent_parser = dynamic_cast<MCMC_XMLParser*>(parent_parser);
	if(!MCMC_XMLParser_parent_parser) error("codon61_sequence_gibbs_sampler_XMLParser: parent parser must be of type MCMC_XMLParser");
	_mcmc = MCMC_XMLParser_parent_parser->get_mcmc();
	new Codon61SequenceGibbsSampler(_mcmc,target,double_weight);
}

int _GENOMEGAMAP_LIBRARY_IS_LOADED = 0;

xsd_string load_genomegaMap_library() {
	if(_GENOMEGAMAP_LIBRARY_IS_LOADED!=0) {
		throw std::runtime_error("load_genomegaMap_library(): library already loaded");
	} else {
		_GENOMEGAMAP_LIBRARY_IS_LOADED = 1;
	}
	// GSL is used, so must set this
	gsl_set_error_handler_off();
	// DISTRIBUTIONS
	distributions_XMLParser::add_child("genomegaMap",&genomegaMap_XMLParser::factory);
	// RANDOM VARIABLES
	data_XMLParser::add_child("codon_count",&codon_count_XMLParser::factory);
	data_XMLParser::add_child("codon_sequence",&codon_sequence_XMLParser::factory);
	parameters_XMLParser::add_child("codon_count",&codon_count_XMLParser::factory);
	parameters_XMLParser::add_child("codon_sequence",&codon_sequence_XMLParser::factory);
	// TRANSFORMATIONS
	transformations_XMLParser::add_child("ny98_pdrm",&ny98_pdrm_XMLParser::factory);
	transformations_XMLParser::add_child("omega_to_gamma_vector",&omega_to_gamma_vector_XMLParser::factory);
	transformations_XMLParser::add_child("parent_dependent_rate_matrix",&parent_dependent_rate_matrix_XMLParser::factory);
	// INFERENCE
	MCMC_XMLParser::add_child("codon61_sequence_gibbs_sampler",&codon61_sequence_gibbs_sampler_XMLParser::factory);
	// SCHEMA
	string s(genomegaMap1_0_xsd_len,' ');
	unsigned int i;
	for(i=0;i<genomegaMap1_0_xsd_len;i++) s[i] = omegaMap1_0_xsd[i];
	return s;
}

} // namespace genomegaMap
