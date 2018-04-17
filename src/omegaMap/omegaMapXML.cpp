/*  Copyright 2018 Daniel Wilson.
 *
 *  omegaMapXML.cpp
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
#include <gsl/gsl_errno.h>
#include <RandomVariables/Continuous.h>
#include <RandomVariables/ContinuousVector.h>
#include <Distributions/DistributionsXML.h>
#include <RandomVariables/RandomVariablesXML.h>
#include <Transformations/TransformationsXML.h>
#include <Inference/InferenceXML.h>
#include <omegaMap/omegaMapXML.h>
//#include <omegaMap/Distributions/BiallelicNCD2.h>
//#include <omegaMap/Distributions/Codon61SequenceStationaryDistribution.h>
//#include <omegaMap/Distributions/mkprf3.h>
//#include <omegaMap/Distributions/mkprf3HMM.h>
//#include <omegaMap/Distributions/mkprf3HMMmag.h>
#include <omegaMap/Distributions/omegaMapNCD.h>
//#include <omegaMap/Distributions/omegaMapNCD2.h>
//#include <omegaMap/Distributions/omegaMapNCD3.h>
//#include <omegaMap/Distributions/omegaMapNCD3P.h>
//#include <omegaMap/Distributions/omegaMapNCDA.h>
//#include <omegaMap/Distributions/omegaMapNCDHMM.h>
//#include <omegaMap/Distributions/omegaMapNCDHMMHybrid.h>
//#include <omegaMap/Distributions/omegaMapNCDOG.h>
#include <omegaMap/Distributions/omegaMapUnlinked.h>
//#include <omegaMap/Distributions/phylogeny.h>
//#include <omegaMap/RandomVariables/BiallelicCodingGroupCount.h>
//#include <omegaMap/RandomVariables/Codon61Alignment.h>
//#include <omegaMap/RandomVariables/Codon61AncestralCount.h>
#include <omegaMap/RandomVariables/Codon61Count.h>
//#include <omegaMap/RandomVariables/Codon61GroupCount.h>
#include <omegaMap/RandomVariables/Codon61Sequence.h>
//#include <omegaMap/Transformations/LogNormalQuantile.h>
//#include <omegaMap/Transformations/LogNormalQuantileVector.h>
//#include <omegaMap/Transformations/mkprf3HMMPathSampler.h>
//#include <omegaMap/Transformations/mkprf3HMMmagPathSampler.h>
//#include <omegaMap/Transformations/NormalQuantile.h>
//#include <omegaMap/Transformations/NormalQuantileMosaic.h>
//#include <omegaMap/Transformations/NormalQuantileVector.h>
#include <omegaMap/Transformations/NY98_PDRM.h>
//#include <omegaMap/Transformations/NY98_TransProb.h>
//#include <omegaMap/Transformations/NY98_TransProbMosaic.h>
#include <omegaMap/Transformations/Omega2GammaVector.h>
#include <omegaMap/Transformations/ParentDependentRateMatrix.h>
#include <omegaMap/Inference/MCMC/omegaMapMoves.h>
#include <omegaMap/omegaMap1.0.xsd.h>
#include <stdexcept>

using namespace gcat;

namespace gcat_omegaMap {

// DISTRIBUTIONS

//biallelicNCD2_XMLParser::biallelicNCD2_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<biallelicNCD2_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 9;
//	const char* attrNames[nattr] = {"id","thetaS0","thetaR0","gamma0","tau0",
//	"thetaS1","thetaR1","gamma1","tau1"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new BiallelicNCD2(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//codon61_sequence_stationary_distribution_XMLParser::codon61_sequence_stationary_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon61_sequence_stationary_distribution_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 2;
//	const char* attrNames[nattr] = {"id","pi"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new Codon61SequenceStationaryDistribution(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//}
//
//mkprf3_XMLParser::mkprf3_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mkprf3_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 13;
//	const char* attrNames[nattr] = {"id","thetaS0","thetaR0","gamma0","tau0",
//		"thetaS1","thetaR1","gamma1","tau1",
//	"thetaS2","thetaR2","gamma2","tau2"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new mkprf3(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//mkprf3_hmm_XMLParser::mkprf3_hmm_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mkprf3_hmm_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 16;
//	const char* attrNames[nattr] = {"id","thetaS0","thetaR0","gamma0","tau0","p0",
//		"thetaS1","thetaR1","gamma1","tau1","p1",
//	"thetaS2","thetaR2","gamma2","tau2","p2"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new mkprf3HMM(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//mkprf3_hmm_mag_XMLParser::mkprf3_hmm_mag_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mkprf3_hmm_mag_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 19;
//	const char* attrNames[nattr] = {"id","thetaS0","thetaR0","gamma0","tau0","p0","q0",
//		"thetaS1","thetaR1","gamma1","tau1","p1","q1",
//	"thetaS2","thetaR2","gamma2","tau2","p2","q2"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new mkprf3HMMmag(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}

omegaMapNCD_XMLParser::omegaMapNCD_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCD_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 3;
	const char* attrNames[nattr] = {"id","mut","sel"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	new omegaMapNCD(sattr[0],getDAG());
	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[2],sattr[2]);
}

//omegaMapNCD2_XMLParser::omegaMapNCD2_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCD2_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 7;
//	const char* attrNames[nattr] = {"id","mut0","sel0","phylo0","mut1","sel1","phylo1"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCD2(sattr[0],getDAG());
//	int i;
//	for(i=1;i<7;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//omegaMapNCD3_XMLParser::omegaMapNCD3_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCD3_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 10;
//	const char* attrNames[nattr] = {"id","mut0","sel0","phylo0","mut1","sel1","phylo1","mut2","sel2","phylo2"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCD3(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//omegaMapNCD3P_XMLParser::omegaMapNCD3P_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCD3P_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 10;
//	const char* attrNames[nattr] = {"id","mut0","sel0","phylo0","mut1","sel1","phylo1","mut2","sel2","phylo2"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCD3P(sattr[0],getDAG());
//	int i;
//	for(i=1;i<nattr;i++) {
//		getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[i],sattr[i]);
//	}
//}
//
//omegaMapNCDA_XMLParser::omegaMapNCDA_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDA_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 3;
//	const char* attrNames[nattr] = {"id","mut","sel"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCDA(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[2],sattr[2]);
//}
//
//omegaMapNCDHMM_XMLParser::omegaMapNCDHMM_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDHMM_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 9;
//	const char* attrNames[nattr] = {"id","anc","theta","kappa","gamma","T","p","pi","gamma_wt"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// theta, kappa, t and p can be specified as numeric, in which case they must be instantiated as Variables
//	double double_theta;
//	if(from_string<double>(double_theta,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_theta);
//		getDAG()->set_constant(sattr[2]);
//	}
//	double double_kappa;
//	if(from_string<double>(double_kappa,sattr[3])) {
//		// Internally-generated name
//		sattr[3] = "_" + sattr[0] + "." + attrNames[3];
//		new ContinuousRV(sattr[3],getDAG(),double_kappa);
//		getDAG()->set_constant(sattr[3]);
//	}
//	double double_T;
//	if(from_string<double>(double_T,sattr[5])) {
//		// Internally-generated name
//		sattr[5] = "_" + sattr[0] + "." + attrNames[5];
//		new ContinuousRV(sattr[5],getDAG(),double_T);
//		getDAG()->set_constant(sattr[5]);
//	}
//	double double_p;
//	if(from_string<double>(double_p,sattr[6])) {
//		// Internally-generated name
//		sattr[6] = "_" + sattr[0] + "." + attrNames[6];
//		new ContinuousRV(sattr[6],getDAG(),double_p);
//		getDAG()->set_constant(sattr[6]);
//	}
//	// Get lengths of anc and gamma
//	Parameter* rv = getDAG()->get_parameter(sattr[1]);
//	if(rv==0) error("omegaMapNCDHMM_XMLParser: could not find anc variable");
//	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("omegaMapNCDHMM_XMLParser: anc variable does not have length property");
//	int int_seqlen = lp->length();
//	rv = getDAG()->get_parameter(sattr[4]);
//	if(rv==0) error("omegaMapNCDHMM_XMLParser: could not find gamma variable");
//	lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("omegaMapNCDHMM_XMLParser: gamma variable does not have length property");
//	int int_ngamma = lp->length();
//	// gamma_wt can take the value 1, in which case it is a vector of 1s
//	if(sattr[8]=="1") {
//		sattr[8] = "_" + sattr[0] + "." + attrNames[8];
//		new ContinuousVectorRV(int_ngamma,sattr[8],getDAG(),vector<double>(int_ngamma,1.0));
//	}
//	else {
//		rv = getDAG()->get_parameter(sattr[8]);
//		if(rv==0) error("omegaMapNCDHMM_XMLParser: could not find gamma_wt variable");
//		lp = dynamic_cast<LengthProperty*>(rv);
//		if(lp==0) error("omegaMapNCDHMM_XMLParser: gamma_wt variable does not have length property");
//		if(lp->length()!=int_ngamma) error("omegaMapNCDHMM_XMLParser: gamma_wt variable has different length to gamma");
//	}
//	new omegaMapNCDHMM(int_seqlen,int_ngamma,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[3],sattr[3]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[4],sattr[4]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[5],sattr[5]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[6],sattr[6]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[7],sattr[7]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[8],sattr[8]);
//}
//
//omegaMapNCDHMMHybrid_XMLParser::omegaMapNCDHMMHybrid_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDHMMHybrid_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 11;
//	const char* attrNames[nattr] = {"id","anc","theta","kappa","gamma1","T","p","pi","gamma1_wt","gamma2","gamma2_wt"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// theta, kappa, t and p can be specified as numeric, in which case they must be instantiated as Variables
//	double double_theta;
//	if(from_string<double>(double_theta,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_theta);
//		getDAG()->set_constant(sattr[2]);
//	}
//	double double_kappa;
//	if(from_string<double>(double_kappa,sattr[3])) {
//		// Internally-generated name
//		sattr[3] = "_" + sattr[0] + "." + attrNames[3];
//		new ContinuousRV(sattr[3],getDAG(),double_kappa);
//		getDAG()->set_constant(sattr[3]);
//	}
//	double double_T;
//	if(from_string<double>(double_T,sattr[5])) {
//		// Internally-generated name
//		sattr[5] = "_" + sattr[0] + "." + attrNames[5];
//		new ContinuousRV(sattr[5],getDAG(),double_T);
//		getDAG()->set_constant(sattr[5]);
//	}
//	double double_p;
//	if(from_string<double>(double_p,sattr[6])) {
//		// Internally-generated name
//		sattr[6] = "_" + sattr[0] + "." + attrNames[6];
//		new ContinuousRV(sattr[6],getDAG(),double_p);
//		getDAG()->set_constant(sattr[6]);
//	}
//	// Get lengths of anc and gamma
//	Parameter* rv = getDAG()->get_parameter(sattr[1]);
//	if(rv==0) error("omegaMapNCDHMMHybrid_XMLParser: could not find anc variable");
//	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("omegaMapNCDHMMHybrid_XMLParser: anc variable does not have length property");
//	int int_seqlen = lp->length();
//	rv = getDAG()->get_parameter(sattr[4]);
//	if(rv==0) error("omegaMapNCDHMMHybrid_XMLParser: could not find gamma1 variable");
//	lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("omegaMapNCDHMMHybrid_XMLParser: gamma1 variable does not have length property");
//	int int_ngamma1 = lp->length();
//	rv = getDAG()->get_parameter(sattr[9]);
//	if(rv==0) error("omegaMapNCDHMMHybrid_XMLParser: could not find gamma2 variable");
//	lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("omegaMapNCDHMMHybrid_XMLParser: gamma2 variable does not have length property");
//	int int_ngamma2 = lp->length();
//	// gamma1_wt can take the value 1, in which case it is a vector of 1s
//	if(sattr[8]=="1") {
//		sattr[8] = "_" + sattr[0] + "." + attrNames[8];
//		new ContinuousVectorRV(int_ngamma1,sattr[8],getDAG(),vector<double>(int_ngamma1,1.0));
//	}
//	else {
//		rv = getDAG()->get_parameter(sattr[8]);
//		if(rv==0) error("omegaMapNCDHMMHybrid_XMLParser: could not find gamma_wt variable");
//		lp = dynamic_cast<LengthProperty*>(rv);
//		if(lp==0) error("omegaMapNCDHMMHybrid_XMLParser: gamma1_wt variable does not have length property");
//		if(lp->length()!=int_ngamma1) error("omegaMapNCDHMM_XMLParser: gamma1_wt variable has different length to gamma1");
//	}
//	// gamma2_wt can take the value 1, in which case it is a vector of 1s
//	if(sattr[10]=="1") {
//		sattr[10] = "_" + sattr[0] + "." + attrNames[10];
//		new ContinuousVectorRV(int_ngamma2,sattr[10],getDAG(),vector<double>(int_ngamma2,1.0));
//	}
//	else {
//		rv = getDAG()->get_parameter(sattr[10]);
//		if(rv==0) error("omegaMapNCDHMMHybrid_XMLParser: could not find gamma_wt variable");
//		lp = dynamic_cast<LengthProperty*>(rv);
//		if(lp==0) error("omegaMapNCDHMMHybrid_XMLParser: gamma2_wt variable does not have length property");
//		if(lp->length()!=int_ngamma2) error("omegaMapNCDHMM_XMLParser: gamma2_wt variable has different length to gamma2");
//	}
//	new omegaMapNCDHMMHybrid(int_seqlen,int_ngamma1,int_ngamma2,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[3],sattr[3]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[4],sattr[4]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[5],sattr[5]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[6],sattr[6]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[7],sattr[7]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[8],sattr[8]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[9],sattr[9]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[10],sattr[10]);
//}
//
//omegaMapNCDOG_XMLParser::omegaMapNCDOG_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDOG_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mut","sel","phylo"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCDOG(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[3],sattr[3]);
//}

omegaMapUnlinked_XMLParser::omegaMapUnlinked_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapUnlinked_XMLParser>(master_parser,parent_parser) {
	// Read in the attributes
	const int nattr = 2;
	const char* attrNames[nattr] = {"id","mut"};
	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
	new omegaMapUnlinked(sattr[0],getDAG());
	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
}

//phylogeny_XMLParser::phylogeny_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<phylogeny_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 2;
//	const char* attrNames[nattr] = {"id","phylo"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new Phylogeny(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_distribution(sattr[0],attrNames[1],sattr[1]);
//}

// RANDOM VARIABLES

//codon_alignment_XMLParser::codon_alignment_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_alignment_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 5;
//	const char* attrNames[nattr] = {"id","distribution","file","format","encoding"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	if(sattr[3]!="fasta") error("codon_alignment_XMLParser: only fasta format supported");
//	// Instantiate the variable
//	if(sattr[4]=="codon61") {
//		new Codon61Alignment(sattr[2],sattr[0],getDAG());
//	}
//	else error("codon_alignment_XMLParser only codon61 encoding supported");
//	if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->set_constant(sattr[0]);
//}
//
//codon_ancestral_count_XMLParser::codon_ancestral_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_ancestral_count_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 6;
//	const char* attrNames[nattr] = {"id","distribution","ancestor","file","format","encoding"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	if(sattr[4]!="fasta") error("codon_ancestral_count_XMLParser: only fasta format supported");
//	// Instantiate the variable
//	if(sattr[5]=="codon61") {
//		new Codon61AncestralCount(sattr[3],sattr[2],sattr[0],getDAG());
//	}
//	else error("codon_ancestral_count_XMLParser only codon61 encoding supported");
//	if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->set_constant(sattr[0]);
//}

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

//codon_group_count_group_XMLParser::codon_group_count_group_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_group_count_group_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 3;
//	const char* attrNames[nattr] = {"name","file","format"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	((codon_group_count_XMLParser*)parent_parser)->groupname.push_back(sattr[0]);
//	((codon_group_count_XMLParser*)parent_parser)->filename.push_back(sattr[1]);
//	if(sattr[2]!="fasta") error("codon_group_count_group_XMLParser: format must equal \"fasta\" currently");
//}
//
//codon_group_count_XMLParser::codon_group_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<codon_group_count_XMLParser>(master_parser,parent_parser), sattr(0), groupname(0), filename(0) {
//	// Read in the attributes
//	const int nattr = 3;
//	const char* attrNames[nattr] = {"id","distribution","encoding"};
//	sattr = attributesToStrings(nattr,attrNames,attrs);
//	if(sattr[2]=="codon61") encoding = codon61;
//	else if	(sattr[2]=="biallelic") encoding = biallelic;
//	else error("codon_group_count_XMLParser: currently encoding must equal \"codon61\" or \"biallelic\"");
//	// Don't instantiate until endElement!
//}
//
//void codon_group_count_XMLParser::implement_startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs) {
//	string element = elementName(localname);
//	if(element=="group") {
//		_child_parser = new codon_group_count_group_XMLParser(uri,localname,qname,attrs,_master_parser,this);
//	}
//	else DAGXMLParserTemplate<codon_group_count_XMLParser>::implement_startElement(uri,localname,qname,attrs);
//}
//
//void codon_group_count_XMLParser::implement_endElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname) {
//	// Instantiate!
//	if(encoding==codon61) new Codon61GroupCount(filename,groupname,sattr[0],getDAG());
//	else if(encoding==biallelic) new BiallelicCodingGroupCount(filename,groupname,sattr[0],getDAG());
//	else error("codon_group_count_XMLParser: expecting one of \"codon61\" or \"biallelic\"");
//	if(sattr[1]!="") getDAG()->assign_distribution_to_random_variable(sattr[0],"distribution",sattr[1]);
//	getDAG()->set_constant(sattr[0]);
//	// Return to default behaviour
//	DAGXMLParserTemplate<codon_group_count_XMLParser>::implement_endElement(uri,localname,qname);
//}

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

//mkprf3_hmm_path_sampler_XMLParser::mkprf3_hmm_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mkprf3_hmm_path_sampler_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","distribution","rv","species"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// species must be specified as an integer
//	int int_species;
//	if(!from_string<int>(int_species,sattr[3])) {
//		string errTxt = "mkprf3_hmm_path_sampler_XMLParser: could not convert species to an integer";
//	}
//	new mkprf3HMMPathSampler(int_species,sattr[2],sattr[1],sattr[0],getDAG());
//}
//
//mkprf3_hmm_mag_path_sampler_XMLParser::mkprf3_hmm_mag_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<mkprf3_hmm_mag_path_sampler_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","distribution","rv","species"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// species must be specified as an integer
//	int int_species;
//	if(!from_string<int>(int_species,sattr[3])) {
//		string errTxt = "mkprf3_hmm_mag_path_sampler_XMLParser: could not convert species to an integer";
//	}
//	new mkprf3HMMmagPathSampler(int_species,sattr[2],sattr[1],sattr[0],getDAG());
//}

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

//ny98_transprob_XMLParser::ny98_transprob_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<ny98_transprob_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 5;
//	const char* attrNames[nattr] = {"id","thetaT","kappa","gamma","pi"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new NY98_TransProb(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[4],sattr[4]);
//}
//
//ny98_transprob_mosaic_XMLParser::ny98_transprob_mosaic_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<ny98_transprob_mosaic_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 6;
//	const char* attrNames[nattr] = {"id","thetaT","kappa","gamma","pi","length"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	int int_length;
//	if(!from_string<int>(int_length,sattr[5])) {
//		RandomVariable* rv = getDAG()->get_random_variable(sattr[5]);
//		if(rv==0) error("ny98_transprob_mosaic_XMLParser: could not convert length to int nor find named variable");
//		LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//		if(lp==0) error("ny98_transprob_mosaic_XMLParser: named variable does not have length property");
//		int_length = lp->length();
//	}
//	if(int_length<=0) error("ny98_transprob_mosaic_XMLParser: length must be a positive integer");
//	new NY98_TransProbMosaic(int_length,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[4],sattr[4]);
//}
//
//omegaMapNCDHMM_path_sampler_XMLParser::omegaMapNCDHMM_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDHMM_path_sampler_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 3;
//	const char* attrNames[nattr] = {"id","distribution","rv"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCDHMMPathSampler(sattr[2],sattr[1],sattr[0],getDAG());
//}
//
//omegaMapNCDHMMHybrid_path_sampler_XMLParser::omegaMapNCDHMMHybrid_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<omegaMapNCDHMMHybrid_path_sampler_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 3;
//	const char* attrNames[nattr] = {"id","distribution","rv"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	new omegaMapNCDHMMHybridPathSampler(sattr[2],sattr[1],sattr[0],getDAG());
//}

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

//log_normal_quantile_function_XMLParser::log_normal_quantile_function_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<log_normal_quantile_function_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mean","sd","quantile"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// mean, sd and quantile can be specified as numeric, in which case they must be instantiated as Variables
//	double double_mean;
//	if(from_string<double>(double_mean,sattr[1])) {
//		// Internally-generated name
//		sattr[1] = "_" + sattr[0] + "." + attrNames[1];
//		new ContinuousRV(sattr[1],getDAG(),double_mean);
//		getDAG()->set_constant(sattr[1]);
//	}
//	double double_sd;
//	if(from_string<double>(double_sd,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[2]);
//	}
//	double double_quantile;
//	if(from_string<double>(double_quantile,sattr[3])) {
//		// Internally-generated name
//		sattr[3] = "_" + sattr[0] + "." + attrNames[3];
//		new ContinuousRV(sattr[3],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[3]);
//	}
//	new LogNormalQuantileTransform(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//}
//
//log_normal_quantile_function_vector_XMLParser::log_normal_quantile_function_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<log_normal_quantile_function_vector_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mean","sd","quantile"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// mean and sd can be specified as numeric, in which case they must be instantiated as Variables
//	double double_mean;
//	if(from_string<double>(double_mean,sattr[1])) {
//		// Internally-generated name
//		sattr[1] = "_" + sattr[0] + "." + attrNames[1];
//		new ContinuousRV(sattr[1],getDAG(),double_mean);
//		getDAG()->set_constant(sattr[1]);
//	}
//	double double_sd;
//	if(from_string<double>(double_sd,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[2]);
//	}
//	// Automatically obtain length of the quantile mosaic
//	int int_length;
//	RandomVariable* rv = getDAG()->get_random_variable(sattr[3]);
//	if(rv==0) error("log_normal_quantile_function_mosaic_XMLParser: could not find named variable for quantile");
//	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("log_normal_quantile_function_mosaic_XMLParser: named variable quantile does not have length property");
//	int_length = lp->length();
//	// Instantiate
//	new LogNormalQuantileVectorTransform(int_length,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//}
//
//normal_quantile_function_XMLParser::normal_quantile_function_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<normal_quantile_function_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mean","sd","quantile"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// mean, sd and quantile can be specified as numeric, in which case they must be instantiated as Variables
//	double double_mean;
//	if(from_string<double>(double_mean,sattr[1])) {
//		// Internally-generated name
//		sattr[1] = "_" + sattr[0] + "." + attrNames[1];
//		new ContinuousRV(sattr[1],getDAG(),double_mean);
//		getDAG()->set_constant(sattr[1]);
//	}
//	double double_sd;
//	if(from_string<double>(double_sd,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[2]);
//	}
//	double double_quantile;
//	if(from_string<double>(double_quantile,sattr[3])) {
//		// Internally-generated name
//		sattr[3] = "_" + sattr[0] + "." + attrNames[3];
//		new ContinuousRV(sattr[3],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[3]);
//	}
//	new NormalQuantileTransform(sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//}
//
//normal_quantile_function_mosaic_XMLParser::normal_quantile_function_mosaic_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<normal_quantile_function_mosaic_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mean","sd","quantile"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// mean and sd can be specified as numeric, in which case they must be instantiated as Variables
//	double double_mean;
//	if(from_string<double>(double_mean,sattr[1])) {
//		// Internally-generated name
//		sattr[1] = "_" + sattr[0] + "." + attrNames[1];
//		new ContinuousRV(sattr[1],getDAG(),double_mean);
//		getDAG()->set_constant(sattr[1]);
//	}
//	double double_sd;
//	if(from_string<double>(double_sd,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[2]);
//	}
//	// Automatically obtain length of the quantile mosaic
//	int int_length;
//	RandomVariable* rv = getDAG()->get_random_variable(sattr[3]);
//	if(rv==0) error("normal_quantile_function_mosaic_XMLParser: could not find named variable for quantile");
//	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("normal_quantile_function_mosaic_XMLParser: named variable quantile does not have length property");
//	int_length = lp->length();
//	// Instantiate
//	new NormalQuantileMosaicTransform(int_length,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//}
//
//normal_quantile_function_vector_XMLParser::normal_quantile_function_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser) : DAGXMLParserTemplate<normal_quantile_function_vector_XMLParser>(master_parser,parent_parser) {
//	// Read in the attributes
//	const int nattr = 4;
//	const char* attrNames[nattr] = {"id","mean","sd","quantile"};
//	vector<string> sattr = attributesToStrings(nattr,attrNames,attrs);
//	// mean and sd can be specified as numeric, in which case they must be instantiated as Variables
//	double double_mean;
//	if(from_string<double>(double_mean,sattr[1])) {
//		// Internally-generated name
//		sattr[1] = "_" + sattr[0] + "." + attrNames[1];
//		new ContinuousRV(sattr[1],getDAG(),double_mean);
//		getDAG()->set_constant(sattr[1]);
//	}
//	double double_sd;
//	if(from_string<double>(double_sd,sattr[2])) {
//		// Internally-generated name
//		sattr[2] = "_" + sattr[0] + "." + attrNames[2];
//		new ContinuousRV(sattr[2],getDAG(),double_sd);
//		getDAG()->set_constant(sattr[2]);
//	}
//	// Automatically obtain length of the quantile mosaic
//	int int_length;
//	RandomVariable* rv = getDAG()->get_random_variable(sattr[3]);
//	if(rv==0) error("normal_quantile_function_mosaic_XMLParser: could not find named variable for quantile");
//	LengthProperty* lp = dynamic_cast<LengthProperty*>(rv);
//	if(lp==0) error("normal_quantile_function_mosaic_XMLParser: named variable quantile does not have length property");
//	int_length = lp->length();
//	// Instantiate
//	new NormalQuantileVectorTransform(int_length,sattr[0],getDAG());
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[1],sattr[1]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[2],sattr[2]);
//	getDAG()->assign_parameter_to_transformation(sattr[0],attrNames[3],sattr[3]);
//}

int _OMEGAMAP_LIBRARY_IS_LOADED = 0;

xsd_string load_omegaMap_library() {
	if(_OMEGAMAP_LIBRARY_IS_LOADED!=0) {
		throw std::runtime_error("load_omegaMap_library(): library already loaded");
	} else {
		_OMEGAMAP_LIBRARY_IS_LOADED = 1;
	}
	// GSL is used, so must set this
	gsl_set_error_handler_off();
	// DISTRIBUTIONS
//	distributions_XMLParser::add_child("biallelicNCD2",&biallelicNCD2_XMLParser::factory);
//	distributions_XMLParser::add_child("codon61_sequence_stationary_distribution",&codon61_sequence_stationary_distribution_XMLParser::factory);
//	distributions_XMLParser::add_child("mkprf3",&mkprf3_XMLParser::factory);
//	distributions_XMLParser::add_child("mkprf3_hmm",&mkprf3_hmm_XMLParser::factory);
//	distributions_XMLParser::add_child("mkprf3_hmm_mag",&mkprf3_hmm_mag_XMLParser::factory);
	distributions_XMLParser::add_child("omegaMapNCD",&omegaMapNCD_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCD2",&omegaMapNCD2_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCD3",&omegaMapNCD3_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCD3P",&omegaMapNCD3P_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCDA",&omegaMapNCDA_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCDHMM",&omegaMapNCDHMM_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCDHMMHybrid",&omegaMapNCDHMMHybrid_XMLParser::factory);
//	distributions_XMLParser::add_child("omegaMapNCDOG",&omegaMapNCDOG_XMLParser::factory);
	distributions_XMLParser::add_child("omegaMapUnlinked",&omegaMapUnlinked_XMLParser::factory);
//	distributions_XMLParser::add_child("phylogeny",&phylogeny_XMLParser::factory);
	// RANDOM VARIABLES
//	data_XMLParser::add_child("codon_alignment",&codon_alignment_XMLParser::factory);
//	data_XMLParser::add_child("codon_ancestral_count",&codon_ancestral_count_XMLParser::factory);
	data_XMLParser::add_child("codon_count",&codon_count_XMLParser::factory);
//	data_XMLParser::add_child("codon_group_count",&codon_group_count_XMLParser::factory);
	data_XMLParser::add_child("codon_sequence",&codon_sequence_XMLParser::factory);
//	parameters_XMLParser::add_child("codon_alignment",&codon_alignment_XMLParser::factory);
//	parameters_XMLParser::add_child("codon_ancestral_count",&codon_ancestral_count_XMLParser::factory);
	parameters_XMLParser::add_child("codon_count",&codon_count_XMLParser::factory);
//	parameters_XMLParser::add_child("codon_group_count",&codon_group_count_XMLParser::factory);
	parameters_XMLParser::add_child("codon_sequence",&codon_sequence_XMLParser::factory);
	// TRANSFORMATIONS
//	transformations_XMLParser::add_child("log_normal_quantile_function",&log_normal_quantile_function_XMLParser::factory);
//	transformations_XMLParser::add_child("log_normal_quantile_function_vector",&log_normal_quantile_function_vector_XMLParser::factory);
//	transformations_XMLParser::add_child("mkprf3_hmm_path_sampler",&mkprf3_hmm_path_sampler_XMLParser::factory);
//	transformations_XMLParser::add_child("mkprf3_hmm_mag_path_sampler",&mkprf3_hmm_mag_path_sampler_XMLParser::factory);
//	transformations_XMLParser::add_child("normal_quantile_function",&normal_quantile_function_XMLParser::factory);
//	transformations_XMLParser::add_child("normal_quantile_function_mosaic",&normal_quantile_function_mosaic_XMLParser::factory);
//	transformations_XMLParser::add_child("normal_quantile_function_vector",&normal_quantile_function_vector_XMLParser::factory);
	transformations_XMLParser::add_child("ny98_pdrm",&ny98_pdrm_XMLParser::factory);
//	transformations_XMLParser::add_child("ny98_transprob",&ny98_transprob_XMLParser::factory);
//	transformations_XMLParser::add_child("ny98_transprob_mosaic",&ny98_transprob_mosaic_XMLParser::factory);
//	transformations_XMLParser::add_child("omegaMapNCDHMM_path_sampler",&omegaMapNCDHMM_path_sampler_XMLParser::factory);
//	transformations_XMLParser::add_child("omegaMapNCDHMMHybrid_path_sampler",&omegaMapNCDHMMHybrid_path_sampler_XMLParser::factory);
	transformations_XMLParser::add_child("omega_to_gamma_vector",&omega_to_gamma_vector_XMLParser::factory);
	transformations_XMLParser::add_child("parent_dependent_rate_matrix",&parent_dependent_rate_matrix_XMLParser::factory);
	// INFERENCE
	MCMC_XMLParser::add_child("codon61_sequence_gibbs_sampler",&codon61_sequence_gibbs_sampler_XMLParser::factory);
	// SCHEMA
	string s(omegaMap1_0_xsd_len,' ');
	unsigned int i;
	for(i=0;i<omegaMap1_0_xsd_len;i++) s[i] = omegaMap1_0_xsd[i];
	return s;
}

} // namespace gcat_omegaMap
