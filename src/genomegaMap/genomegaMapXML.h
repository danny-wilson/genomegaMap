/*  Copyright 2018 Daniel Wilson.
 *
 *  genomegaMapXML.h
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
#ifndef _GENOMEGAMAP_XML_H_
#define _GENOMEGAMAP_XML_H_
#include <DAG/DAGXMLParser.h>
#include <Inference/MCMC/MCMC.h>

using namespace gcat;

namespace genomegaMap {
	
// DISTRIBUTIONS

/*	<xs:element name="genomegaMap" substitutionGroup="abstract_distribution">
			<xs:complexType>
				<xs:attribute name="id" type="xs:string" use="required"/>
				<xs:attribute name="mut" type="xs:string" use="required"/>
			</xs:complexType>
		</xs:element>
 */
class genomegaMap_XMLParser : public DAGXMLParserTemplate<genomegaMap_XMLParser> {
public:
	genomegaMap_XMLParser(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const Attributes& attrs, DAGXMLMasterParser* const master_parser, DAGXMLParser* const parent_parser);
};

// RANDOM VARIABLES

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
xsd_string load_genomegaMap_library();
	
} // namespace genomegaMap

#endif//_GENOMEGAMAP_XML_H_
