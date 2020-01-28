/*  Copyright 2018 Daniel Wilson.
 *
 *  omegaMapXML.h
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
#ifndef _OMEGAMAP_XML_H_
#define _OMEGAMAP_XML_H_
#include <DAG/DAGXMLParser.h>
#include <Inference/MCMC/MCMC.h>

using namespace gcat;

namespace gcat_omegaMap {
	
// DISTRIBUTIONS

/*	<xs:element name="biallelicNCD2">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaS0" type="xs:string" use="required"/>
			<xs:attribute name="thetaR0" type="xs:string" use="required"/>
			<xs:attribute name="gamma0" type="xs:string" use="required"/>
			<xs:attribute name="tau0" type="xs:string" use="required"/>
			<xs:attribute name="thetaS1" type="xs:string" use="required"/>
			<xs:attribute name="thetaR1" type="xs:string" use="required"/>
			<xs:attribute name="gamma1" type="xs:string" use="required"/>
			<xs:attribute name="tau1" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class biallelicNCD2_XMLParser : public DAGXMLParserTemplate<biallelicNCD2_XMLParser> {
public:
	biallelicNCD2_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="codon61_sequence_stationary_distribution">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
	</xs:element>
 */
class codon61_sequence_stationary_distribution_XMLParser : public DAGXMLParserTemplate<codon61_sequence_stationary_distribution_XMLParser> {
public:
	codon61_sequence_stationary_distribution_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="mkprf3">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaS0" type="xs:string" use="required"/>
			<xs:attribute name="thetaR0" type="xs:string" use="required"/>
			<xs:attribute name="gamma0" type="xs:string" use="required"/>
			<xs:attribute name="tau0" type="xs:string" use="required"/>
			<xs:attribute name="thetaS1" type="xs:string" use="required"/>
			<xs:attribute name="thetaR1" type="xs:string" use="required"/>
			<xs:attribute name="gamma1" type="xs:string" use="required"/>
			<xs:attribute name="tau1" type="xs:string" use="required"/>
			<xs:attribute name="thetaS2" type="xs:string" use="required"/>
			<xs:attribute name="thetaR2" type="xs:string" use="required"/>
			<xs:attribute name="gamma2" type="xs:string" use="required"/>
			<xs:attribute name="tau2" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class mkprf3_XMLParser : public DAGXMLParserTemplate<mkprf3_XMLParser> {
public:
	mkprf3_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="mkprf3_hmm">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaS0" type="xs:string" use="required"/>
			<xs:attribute name="thetaR0" type="xs:string" use="required"/>
			<xs:attribute name="gamma0" type="xs:string" use="required"/>
			<xs:attribute name="tau0" type="xs:string" use="required"/>
			<xs:attribute name="p0" type="xs:string" use="required"/>
			<xs:attribute name="thetaS1" type="xs:string" use="required"/>
			<xs:attribute name="thetaR1" type="xs:string" use="required"/>
			<xs:attribute name="gamma1" type="xs:string" use="required"/>
			<xs:attribute name="tau1" type="xs:string" use="required"/>
			<xs:attribute name="p1" type="xs:string" use="required"/>
			<xs:attribute name="thetaS2" type="xs:string" use="required"/>
			<xs:attribute name="thetaR2" type="xs:string" use="required"/>
			<xs:attribute name="gamma2" type="xs:string" use="required"/>
			<xs:attribute name="tau2" type="xs:string" use="required"/>
			<xs:attribute name="p2" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class mkprf3_hmm_XMLParser : public DAGXMLParserTemplate<mkprf3_hmm_XMLParser> {
public:
	mkprf3_hmm_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="mkprf3_hmm_mag">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaS0" type="xs:string" use="required"/>
			<xs:attribute name="thetaR0" type="xs:string" use="required"/>
			<xs:attribute name="gamma0" type="xs:string" use="required"/>
			<xs:attribute name="tau0" type="xs:string" use="required"/>
			<xs:attribute name="p0" type="xs:string" use="required"/>
			<xs:attribute name="q0" type="xs:string" use="required"/>
			<xs:attribute name="thetaS1" type="xs:string" use="required"/>
			<xs:attribute name="thetaR1" type="xs:string" use="required"/>
			<xs:attribute name="gamma1" type="xs:string" use="required"/>
			<xs:attribute name="tau1" type="xs:string" use="required"/>
			<xs:attribute name="p1" type="xs:string" use="required"/>
			<xs:attribute name="q1" type="xs:string" use="required"/>
			<xs:attribute name="thetaS2" type="xs:string" use="required"/>
			<xs:attribute name="thetaR2" type="xs:string" use="required"/>
			<xs:attribute name="gamma2" type="xs:string" use="required"/>
			<xs:attribute name="tau2" type="xs:string" use="required"/>
			<xs:attribute name="p2" type="xs:string" use="required"/>
			<xs:attribute name="q2" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class mkprf3_hmm_mag_XMLParser : public DAGXMLParserTemplate<mkprf3_hmm_mag_XMLParser> {
public:
	mkprf3_hmm_mag_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCD">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut" type="xs:string" use="required"/>
			<xs:attribute name="sel" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/
class omegaMapNCD_XMLParser : public DAGXMLParserTemplate<omegaMapNCD_XMLParser> {
public:
	omegaMapNCD_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCD2">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut0" type="xs:string" use="required"/>
			<xs:attribute name="sel0" type="xs:string" use="required"/>
			<xs:attribute name="phylo0" type="xs:string" use="required"/>
			<xs:attribute name="mut1" type="xs:string" use="required"/>
			<xs:attribute name="sel1" type="xs:string" use="required"/>
			<xs:attribute name="phylo1" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCD2_XMLParser : public DAGXMLParserTemplate<omegaMapNCD2_XMLParser> {
public:
	omegaMapNCD2_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCD3">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut0" type="xs:string" use="required"/>
			<xs:attribute name="sel0" type="xs:string" use="required"/>
			<xs:attribute name="phylo0" type="xs:string" use="required"/>
			<xs:attribute name="mut1" type="xs:string" use="required"/>
			<xs:attribute name="sel1" type="xs:string" use="required"/>
			<xs:attribute name="phylo1" type="xs:string" use="required"/>
			<xs:attribute name="mut2" type="xs:string" use="required"/>
			<xs:attribute name="sel2" type="xs:string" use="required"/>
			<xs:attribute name="phylo2" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCD3_XMLParser : public DAGXMLParserTemplate<omegaMapNCD3_XMLParser> {
public:
	omegaMapNCD3_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCD3P">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut0" type="xs:string" use="required"/>
			<xs:attribute name="sel0" type="xs:string" use="required"/>
			<xs:attribute name="phylo0" type="xs:string" use="required"/>
			<xs:attribute name="mut1" type="xs:string" use="required"/>
			<xs:attribute name="sel1" type="xs:string" use="required"/>
			<xs:attribute name="phylo1" type="xs:string" use="required"/>
			<xs:attribute name="mut2" type="xs:string" use="required"/>
			<xs:attribute name="sel2" type="xs:string" use="required"/>
			<xs:attribute name="phylo2" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCD3P_XMLParser : public DAGXMLParserTemplate<omegaMapNCD3P_XMLParser> {
public:
	omegaMapNCD3P_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDA">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut" type="xs:string" use="required"/>
			<xs:attribute name="sel" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDA_XMLParser : public DAGXMLParserTemplate<omegaMapNCDA_XMLParser> {
public:
	omegaMapNCDA_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDHMM">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="anc" type="xs:string" use="required"/>
			<xs:attribute name="theta" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="gamma" type="xs:string" use="required"/>
			<xs:attribute name="T" type="xs:string" use="required"/>
			<xs:attribute name="p" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDHMM_XMLParser : public DAGXMLParserTemplate<omegaMapNCDHMM_XMLParser> {
public:
	omegaMapNCDHMM_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDHMMHybrid">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="anc" type="xs:string" use="required"/>
			<xs:attribute name="theta" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="gamma1" type="xs:string" use="required"/>
			<xs:attribute name="gamma2" type="xs:string" use="required"/>
			<xs:attribute name="T" type="xs:string" use="required"/>
			<xs:attribute name="p" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
			<xs:attribute name="gamma1_wt" type="xs:string" default="1"/>
			<xs:attribute name="gamma2_wt" type="xs:string" default="1"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDHMMHybrid_XMLParser : public DAGXMLParserTemplate<omegaMapNCDHMMHybrid_XMLParser> {
public:
	omegaMapNCDHMMHybrid_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDOG">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut" type="xs:string" use="required"/>
			<xs:attribute name="sel" type="xs:string" use="required"/>
			<xs:attribute name="phylo" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDOG_XMLParser : public DAGXMLParserTemplate<omegaMapNCDOG_XMLParser> {
public:
	omegaMapNCDOG_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapUnlinked">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mut" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapUnlinked_XMLParser : public DAGXMLParserTemplate<omegaMapUnlinked_XMLParser> {
public:
	omegaMapUnlinked_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="phylogeny">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="phylo" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class phylogeny_XMLParser : public DAGXMLParserTemplate<phylogeny_XMLParser> {
public:
	phylogeny_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

// RANDOM VARIABLES

/*	<xs:element name="codon_alignment">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" default=""/>
			<xs:attribute name="file" type="xs:string" use="required"/>
			<xs:attribute name="format" type="xs:string" default="fasta"/>
			<xs:attribute name="encoding" type="xs:string" default="codon61"/>
		</xs:complexType>
	</xs:element>
 */
class codon_alignment_XMLParser : public DAGXMLParserTemplate<codon_alignment_XMLParser> {
public:
	codon_alignment_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="codon_ancestral_count">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" default=""/>
			<xs:attribute name="ancestor" type="xs:integer" use="required"/>
			<xs:attribute name="file" type="xs:string" use="required"/>
			<xs:attribute name="format" type="xs:string" default="fasta"/>
			<xs:attribute name="encoding" type="xs:string" default="codon61"/>
		</xs:complexType>
	</xs:element>
 */
class codon_ancestral_count_XMLParser : public DAGXMLParserTemplate<codon_ancestral_count_XMLParser> {
public:
	codon_ancestral_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="codon_count">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" default=""/>
			<xs:attribute name="file" type="xs:string" use="required"/>
			<xs:attribute name="format" type="xs:string" default="fasta"/>
			<xs:attribute name="encoding" type="xs:string" default="codon61"/>
		</xs:complexType>
	</xs:element>
 */
class codon_count_XMLParser : public DAGXMLParserTemplate<codon_count_XMLParser> {
public:
	codon_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:complexType name="codon_group_count_type">
		<xs:choice minOccurs="0" maxOccurs="unbounded">
			<xs:element name="group">
				<xs:complexType>
					<xs:attribute name="name" type="xs:string" use="required"/>
					<xs:attribute name="file" type="xs:string" use="required"/>
					<xs:attribute name="format" type="xs:string" default="fasta"/>
				</xs:complexType>
			</xs:element>
		</xs:choice>
	</xs:complexType>

	<xs:element name="codon_group_count">
		<xs:complexType>
			<xs:complexContent>
				<xs:extension base="codon_group_count_type">
					<xs:attribute name="id" type="xs:string" use="required"/>
					<xs:attribute name="distributon" type="xs:string" default=""/>
					<xs:attribute name="encoding" type="xs:string" default="codon61"/>
				</xs:extension>
			</xs:complexContent>
		</xs:complexType>
	</xs:element>
*/
class codon_group_count_group_XMLParser : public DAGXMLParserTemplate<codon_group_count_group_XMLParser> {
public:
	codon_group_count_group_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

class codon_group_count_XMLParser : public DAGXMLParserTemplate<codon_group_count_XMLParser> {
private:
	friend class codon_group_count_group_XMLParser;
	vector<string> sattr, groupname, filename;
	enum {codon61=0, biallelic} encoding;
public:
	codon_group_count_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	void implement_startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs);
	void implement_endElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname);
};

/*	<xs:simpleType name="string_list">
		<xs:list itemType="xs:string"/>
	</xs:simpleType>

	<xs:element name="codon_sequence">
		<xs:complexType>
			<xs:simpleContent>
				<xs:extension base="string_list">
					<xs:attribute name="id" type="xs:string" use="required"/>
					<xs:attribute name="distribution" type="xs:string" default=""/>
					<xs:attribute name="encoding" type="xs:string" default="codon61"/>
				</xs:extension>
			</xs:simpleContent>
		</xs:complexType>
	</xs:element>
 */
class codon_sequence_XMLParser : public DAGXMLParserTemplate<codon_sequence_XMLParser> {
	vector<string> sattr;
public:
	codon_sequence_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
	void implement_characters(const XMLCh* const chars, const XMLSize_t length);
};

// TRANSFORMATIONS

/*	<xs:element name="mkprf3_hmm_path_sampler">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" use="required"/>
			<xs:attribute name="rv" type="xs:string" use="required"/>
			<xs:attribute name="species" type="xs:integer" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class mkprf3_hmm_path_sampler_XMLParser : public DAGXMLParserTemplate<mkprf3_hmm_path_sampler_XMLParser> {
public:
	mkprf3_hmm_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="mkprf3_hmm_mag_path_sampler">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" use="required"/>
			<xs:attribute name="rv" type="xs:string" use="required"/>
			<xs:attribute name="species" type="xs:integer" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class mkprf3_hmm_mag_path_sampler_XMLParser : public DAGXMLParserTemplate<mkprf3_hmm_mag_path_sampler_XMLParser> {
public:
	mkprf3_hmm_mag_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="ny98_pdrm">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="theta" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="omega" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
			<xs:attribute name="length" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class ny98_pdrm_XMLParser : public DAGXMLParserTemplate<ny98_pdrm_XMLParser> {
public:
	ny98_pdrm_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="ny98_transprob">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaT" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="gamma" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class ny98_transprob_XMLParser : public DAGXMLParserTemplate<ny98_transprob_XMLParser> {
public:
	ny98_transprob_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="ny98_transprob_mosaic">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="thetaT" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="gamma" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
			<xs:attribute name="length" type="xs:length" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class ny98_transprob_mosaic_XMLParser : public DAGXMLParserTemplate<ny98_transprob_mosaic_XMLParser> {
public:
	ny98_transprob_mosaic_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDHMM_path_sampler_XMLParser">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" use="required"/>
			<xs:attribute name="rv" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDHMM_path_sampler_XMLParser : public DAGXMLParserTemplate<omegaMapNCDHMM_path_sampler_XMLParser> {
public:
	omegaMapNCDHMM_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omegaMapNCDHMMHybrid_path_sampler">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="distribution" type="xs:string" use="required"/>
			<xs:attribute name="rv" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omegaMapNCDHMMHybrid_path_sampler_XMLParser : public DAGXMLParserTemplate<omegaMapNCDHMMHybrid_path_sampler_XMLParser> {
public:
	omegaMapNCDHMMHybrid_path_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="omega_to_gamma_vector">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="omega" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class omega_to_gamma_vector_XMLParser : public DAGXMLParserTemplate<omega_to_gamma_vector_XMLParser> {
public:
	omega_to_gamma_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="parent_dependent_rate_matrix">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="theta" type="xs:string" use="required"/>
			<xs:attribute name="kappa" type="xs:string" use="required"/>
			<xs:attribute name="pi" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/
class parent_dependent_rate_matrix_XMLParser : public DAGXMLParserTemplate<parent_dependent_rate_matrix_XMLParser> {
public:
	parent_dependent_rate_matrix_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="log_normal_quantile_function">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" default="0"/>
			<xs:attribute name="sd" type="xs:string" default="1"/>
			<xs:attribute name="quantile" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class log_normal_quantile_function_XMLParser : public DAGXMLParserTemplate<log_normal_quantile_function_XMLParser> {
public:
	log_normal_quantile_function_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="log_normal_quantile_function_vector">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" default="0"/>
			<xs:attribute name="sd" type="xs:string" default="1"/>
			<xs:attribute name="quantile" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class log_normal_quantile_function_vector_XMLParser : public DAGXMLParserTemplate<log_normal_quantile_function_vector_XMLParser> {
public:
	log_normal_quantile_function_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="normal_quantile_function">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" default="0"/>
			<xs:attribute name="sd" type="xs:string" default="1"/>
			<xs:attribute name="quantile" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class normal_quantile_function_XMLParser : public DAGXMLParserTemplate<normal_quantile_function_XMLParser> {
public:
	normal_quantile_function_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*	<xs:element name="normal_quantile_function_mosaic">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" default="0"/>
			<xs:attribute name="sd" type="xs:string" default="1"/>
			<xs:attribute name="quantile" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
*/
class normal_quantile_function_mosaic_XMLParser : public DAGXMLParserTemplate<normal_quantile_function_mosaic_XMLParser> {
public:
	normal_quantile_function_mosaic_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

/*<xs:element name="normal_quantile_function_vector">
		<xs:complexType>
			<xs:attribute name="id" type="xs:string" use="required"/>
			<xs:attribute name="mean" type="xs:string" default="0"/>
			<xs:attribute name="sd" type="xs:string" default="1"/>
			<xs:attribute name="quantile" type="xs:string" use="required"/>
		</xs:complexType>
	</xs:element>
 */
class normal_quantile_function_vector_XMLParser : public DAGXMLParserTemplate<normal_quantile_function_vector_XMLParser> {
public:
	normal_quantile_function_vector_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

// INFERENCE

/*	<xs:element name="codon61_sequence_gibbs_sampler">
		<xs:complexType>
			<xs:attribute name="parameter" type="xs:string" use="required"/>
			<xs:attribute name="weight" type="xs:decimal" default="1"/>
		</xs:complexType>
	</xs:element>
 */
class codon61_sequence_gibbs_sampler_XMLParser : public DAGXMLParserTemplate<codon61_sequence_gibbs_sampler_XMLParser> {
protected:
	MCMC* _mcmc;
public:
	codon61_sequence_gibbs_sampler_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

// Load the library
xsd_string load_omegaMap_library();
	
} // namespace gcat_omegaMap

#endif//_OMEGAMAP_XML_H_
